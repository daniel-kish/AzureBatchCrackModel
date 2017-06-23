
class Crack;

namespace CrackModel{
	class PreruptureArea
	{
	public:
		PreruptureArea()
		{
		}

		~PreruptureArea()
		{
		}
		double getLen()
		{
			return len;
		}
		double calcLen(double crackLen, double maxPermissibleLen){}
	private:
		double len;
		double initLen;
		double maxToInitRatio;
		double alpha, beta;

	};
	PreruptureArea
}