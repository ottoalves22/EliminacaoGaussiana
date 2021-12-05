import java.util.*;

// classe que representa uma matriz de valores do tipo double.

class Matriz {

	// constante para ser usada na comparacao de valores double.
	// Se a diferenca absoluta entre dois valores double for menor
	// do que o valor definido por esta constante, eles devem ser
	// considerados iguais.
	public static final double SMALL = 0.000001;
	
	private int lin, col;	
	private double [][] m;

	// metodo estatico que cria uma matriz identidade de tamanho n x n.

	public static Matriz identidade(int n){

		Matriz mat = new Matriz(n, n);
		for(int i = 0; i < mat.lin; i++) mat.m[i][i] = 1;
		return mat;
	}	

	// construtor que cria uma matriz de n linhas por m colunas com todas as entradas iguais a zero.

	public Matriz(int n, int m){

		this.lin = n;
		this.col = m;
		this.m = new double[lin][col];
	}

	public void set(int i, int j, double valor){

		m[i][j] = valor;
	}

	public double get(int i, int j){

		return m[i][j];
	}

	// metodo que imprime as entradas da matriz.

	public void imprime(){

		for(int i = 0; i < lin; i++){

			for(int j = 0; j < col; j++){
	
				System.out.printf("%7.2f ", m[i][j]);
			}

			System.out.println();
		}
	}

	public void imprimeDeterminante(){

		for(int i = 0; i < lin; i++){
			System.out.print("|");
			for(int j = 0; j < col; j++){

				System.out.printf("%7.2f ", m[i][j]);
			}
			System.out.print("|");
			System.out.println();
		}
	}

	// metodo que imprime a matriz expandida formada pela combinacao da matriz que 
	// chama o metodo com a matriz "agregada" recebida como parametro. Ou seja, cada 
	// linha da matriz impressa possui as entradas da linha correspondente da matriz 
	// que chama o metodo, seguida das entradas da linha correspondente em "agregada".

	public void imprime(Matriz agregada){

		for(int i = 0; i < lin; i++){

			for(int j = 0; j < col; j++){
	
				System.out.printf("%7.2f ", m[i][j]);
			}

			System.out.print(" |");

			for(int j = 0; j < agregada.col; j++){
	
				System.out.printf("%7.2f ", agregada.m[i][j]);
			}

			System.out.println();
		}
	}
	
	// metodo que troca as linhas i1 e i2 de lugar.

	private void trocaLinha(int i1, int i2){
        int i = 0;
        double aux;
        while(i<this.col){
			aux = this.m[i1-1][i];
			set(i1-1, i, this.m[i2-1][i]);
			set(i2-1, i, aux);
        }
	}

	// metodo que multiplica as entradas da linha i pelo escalar k

	private void multiplicaLinha(int i, double k){
        int j = 0;
        while(j<this.col){
			set(i, j, (get(i, j)*k));
            j++;
        }
	}

	// metodo que faz a seguinte combinacao de duas linhas da matriz:
	//	
	// 	(linha i1) = (linha i1) + (linha i2 * k)
	//

	private void combinaLinhas(int i1, int i2, double k){
        int i = 0;
        while(i<this.col){
            set(i1, i, (get(i1, i)+get(i2, i)*k));
			i++;
        }
	}

	// metodo que procura, a partir da linha ini, a linha com uma entrada nao nula que
	// esteja o mais a esquerda possivel dentre todas as linhas. Os indices da linha e da 
	// coluna referentes a entrada nao nula encontrada sao devolvidos como retorno do metodo.
	// Este metodo ja esta pronto para voces usarem na implementacao da eliminacao gaussiana
	// e eleminacao de Gauss-Jordan.

	private int [] encontraLinhaPivo(int ini){

		int pivo_col, pivo_lin;

		pivo_lin = lin;
		pivo_col = col;

		for(int i = ini; i < lin; i++){
		
			int j;
			
			for(j = 0; j < col; j++) if(Math.abs(m[i][j]) > 0) break;

			if(j < pivo_col) {

				pivo_lin = i;
				pivo_col = j;
			}
		}

		return new int [] { pivo_lin, pivo_col };
	}


	// metodo que implementa a eliminacao gaussiana, que coloca a matriz (que chama o metodo)
	// na forma escalonada. As operacoes realizadas para colocar a matriz na forma escalonada
	// tambem devem ser aplicadas na matriz "agregada" caso esta seja nao nula. Este metodo
	// tambem deve calcular e devolver o determinante da matriz que invoca o metodo. Assumimos
	// que a matriz que invoca este metodo eh uma matriz quadrada.

	public double[] formaEscalonada(Matriz ma, Matriz ag){

		double[][] m = ma.m;
		double[][] agregada = ag.m;

		int n = agregada.length;

		for(int i=0; i<n; i++){
			int maximo = i;

			for(int j=i+1; j<n; j++){
				if(m[j][i] > m[maximo][i]){
					maximo = j;
				}
			}

			double[] temp = m[i];
			m[i] = m[maximo];
			m[maximo] = temp;

			double temp2 = agregada[i][0]; // acho que ta pegando da coluna correta
			agregada[i][0] = agregada[maximo][0];
			agregada[maximo][0] = temp2;

			for(int j=i+1; j<n ; j++){
				double mL = m[j][i]/m[i][i];
				agregada[j][0] = mL * agregada[i][0];
				for(int k=0; k<n; k++){
					m[j][k] -= mL * m[i][k];
				}
			}
		}

		//vetor de solucoes
		double[] x = new double[n];
		for (int i = n - 1; i >= 0; i--) {
			double sum = 0.0;
			for (int j = i + 1; j < n; j++) {
				sum += m[i][j] * x[j];
			}
			x[i] = (agregada[i][0] - sum) / m[i][i];
		}

		return x; // ARRANJO COM VALORES DE X E Y
	}

	// fase de eliminacao:
		//etapa_eliminacao(agregada); // falta fazer com a porra da MAtriz Agregada (sera q ela é o b do [A | b]?)

	// metodo que implementa a eliminacao de Gauss-Jordan, que coloca a matriz (que chama o metodo)
	// na forma escalonada reduzida. As operacoes realizadas para colocar a matriz na forma escalonada 
	// reduzida tambem devem ser aplicadas na matriz "agregada" caso esta seja nao nula. Assumimos que
	// a matriz que invoca esta metodo eh uma matriz quadrada. Não se pode assumir, contudo, que esta
	// matriz ja esteja na forma escalonada (mas voce pode usar o metodo acima para isso).

	public void formaEscalonadaReduzida(Matriz agregada){

		// TODO: implementar este metodo.		
	}
}

// Classe "executavel".

public class EP1 {

	// metodo principal.

	public static void main(String [] args){

		Scanner in = new Scanner(System.in);	// Scanner para facilitar a leitura de dados a partir da entrada padrao.
		String operacao = in.next();		// le, usando o scanner, a string que determina qual operacao deve ser realizada.
		int n = in.nextInt();			// le a dimensão da matriz a ser manipulada pela operacao escolhida.

		// TODO: completar este metodo.
		Matriz m = new Matriz(n, n);
		Matriz aux = new Matriz(n,1);

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				int nTemp = in.nextInt();
				m.set(i, j, nTemp);
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < 1; j++){
				int nTemp = in.nextInt();
				aux.set(i, j, nTemp);
			}
		}
		/*
		for(int i = 0; i < n; i++){
			for(int j = 0; j <= n; j++){
				if (j < n) {
					m.set(i, j, in.nextInt());
				} else {
					int nTemp = in.nextInt();
					aux.set(1, n-j, nTemp);
				}
			}
		}*/





		if("resolve".equals(operacao)){


		}
		else if("inverte".equals(operacao)){

		}
		else if("determinante".equals(operacao)){
			System.out.println("determinante:");
			m.formaEscalonada(m, aux);
			m.imprime();
		}
		else {
			System.out.println("Operação desconhecida!");
			System.exit(1);
		}
	}
}
