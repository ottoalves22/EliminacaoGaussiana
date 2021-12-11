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
		//fazendo substituicoes
		double[] temp = m[i1];
		m[i1] = m[i2];
		m[i2] = temp;
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




	public double formaEscalonada(Matriz agregada){
		//m e a matriz e agregada_aux e o vetor de resultados de m

		double[] agregada_aux = new double[this.lin];

		int n = agregada.m.length;

		for(int i=0; i<n; i++){
			agregada_aux[i] = this.m[i][n];
		}

		int i=0;
		while(i<n){
			//pivotando
			int maximo = this.encontraLinhaPivo(i)[0];
			this.trocaLinha(i, maximo);
			//troca de linha no vetor resposta
			double temp2 = agregada_aux[i]; // ta pegando da coluna correta?
			agregada_aux[i] = agregada_aux[maximo];
			agregada_aux[maximo] = temp2;

			//faz a eliminacao
			for(int j=i+1; j<n ; j++){
				double mL = this.m[j][i]/this.m[i][i];
				if(agregada_aux!=null){
					agregada_aux[j] -= mL * agregada_aux[i];
				}
				for(int k=0; k<n; k++){
					this.m[j][k] -= mL * this.m[i][k];
				}
			}
			i++;
		}

		//calculo das solucoes x y z...
		double[] x = new double[n];
		for (int z = n - 1; z >= 0; z--) {
			double sum = 0.0;
			for (int j = z + 1; j < n; j++) {
				sum += this.m[z][j] * x[j];
			}
			x[z] = (agregada_aux[z]-sum)/this.m[z][z];
		}

		//exibe os resultados
		for(i=0; i<n; i++){
			System.out.println(x[i]);
		}

		//atribui as arranjos nas Matriz
		for(i=0; i<n; i++){
			agregada.m[i][0] = agregada_aux[i];
		}

		return determinante(this);
	}


	public double determinante(Matriz ma){
		Matriz aux_m = ma;
		double[][] matriz_entrada = aux_m.m;
		int order = matriz_entrada.length;
		double [][] m = new double[order-1][order-1];
		double somatorio = 0;
		int i,j,k;

		if(order==1){
			somatorio=matriz_entrada[0][0];
		}else{
			for(int z=0; z<order; z++){
				int s = 0;
				k=0;
				for(i=1; i<order; i++){
					for(j=0; j<order; j++){
						if(j!=z) {
							m[s][k] = matriz_entrada[i][j];
							k++;
							if(k==order-1) {
								s++;
								k = 0;
							}
						}
					}
				}
				aux_m.m = m;
				aux_m.lin = order-1;
				aux_m.col = order-1;
				somatorio = somatorio + matriz_entrada[0][z] * Math.pow((-1), z) * determinante(aux_m);
			}
		}
		return somatorio;
	}

	public void triangularSuperior(){

	}

	public void formaEscalonadaReduzida(Matriz agregada){
		Matriz aux = new Matriz(this.lin, this.col);

		for(int i=0; i<this.m.length; i++){
			for(int j=0; j<this.m.length; j++) {
				aux.set(i, j, this.get(i, j));
			}
		}

		double res = this.determinante(aux); // a determinante ta zoando a matriz...
		if(res==0.0){
			System.out.println("matriz singular");
			return;
		} else {

			int n = this.m.length;
			double x[][] = new double[n][n];
			double b[][] = new double[n][n];
			for(int i=0; i<n; ++i){
				b[i][i] = 1;
			}

			//GAUSS (IGUAL AO DO METODO DE CIMA) TRIANGULAR INFERIOR
			double[] agregada_aux = new double[this.lin];


			int i=0;
			while(i<n){
				//pivotando
				int maximo = this.encontraLinhaPivo(i)[0];
				this.trocaLinha(i, maximo);
				//troca de linha no vetor resposta

				//faz a eliminacao
				for(int j=i+1; j<n ; j++){
					double mL = this.m[j][i]/this.m[i][i];
					for(int k=0; k<n; k++){
						this.m[j][k] -= mL * this.m[i][k];
					}
				}
				i++;
			}

			//hora de fazer a triangular infeior
			i=0;
			while(i<n){
				//pivotando
				int maximo = this.encontraLinhaPivo(i)[0];
				this.trocaLinha(i, maximo);
				//troca de linha no vetor resposta

				//faz a eliminacao
				for(int j=i+1; j<n ; j++){
					double mL = this.m[j][i]/this.m[i][i];
					this.m[i][j] = mL;
					for(int k=i+1; k<n; k++){
						this.m[j][k] -= mL * this.m[i][k];
					}
				}
				i++;
			}
			//matriz triangulada certinho

			for(int v=0; v<n-1; ++v){
				for(int j=v+1; j<n; ++j){
					for(int k=0; k<n; ++k){
						b[j][k] -= this.m[j][v]*b[v][k];
					}
				}
			}

			for (int z=0; z<n; ++z)
			{
				x[n-1][z] = b[n-1][z]/this.m[n-1][n-1];
				for (int j=n-2; j>=0; --j)
				{
					x[j][z] = b[j][z];
					for (int k=j+1; k<n; ++k)
					{
						x[j][z] -= this.m[j][k]*x[k][z];
					}
					x[j][z] /= this.m[j][j];
				}
			}

			this.m = x;
			this.imprime();
		}
	}
}

// Classe "executavel".

public class EP1 {

	// metodo principal.

	public static void main(String [] args){

		Scanner in = new Scanner(System.in);	// Scanner para facilitar a leitura de dados a partir da entrada padrao.
		String operacao = in.next();		// le, usando o scanner, a string que determina qual operacao deve ser realizada.
		int n = in.nextInt();			// le a dimensão da matriz a ser manipulada pela operacao escolhida.

		if("resolve".equals(operacao)){
			Matriz m1 = new Matriz(n, n+1);

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n+1; j++){
					int nTemp = in.nextInt();
					m1.set(i, j, nTemp);
				}
			}

			Matriz agregada = new Matriz(n, 1);

			double resultado = m1.formaEscalonada(agregada); //TODO ARRUMA ESSA PORRA
		}
		else if("inverte".equals(operacao)){
			Matriz m = new Matriz(n, n);

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					int nTemp = in.nextInt();
					m.set(i, j, nTemp);
				}
			}
			Matriz agregada = new Matriz(n, 1);

			m.formaEscalonadaReduzida(agregada);
		}
		else if("determinante".equals(operacao)){
			Matriz m = new Matriz(n, n);

			for(int i = 0; i < n; i++){
				for(int j = 0; j < n; j++){
					int nTemp = in.nextInt();
					m.set(i, j, nTemp);
				}
			}
			System.out.println("determinante:");
			double resultado = m.determinante(m);
			System.out.println(resultado);
		}
		else {
			System.out.println("Operação desconhecida!");
			System.exit(1);
		}
	}
}
