#include "stdafx.h"
#include "CppUnitTest.h"

#include "definitions.h"
#include "excess_power_convergence.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace OpenFFD_Tests
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			// TODO: Your test code here
            epc::checkTrapz(0, bmp_real(1.0));
		}

	};
}