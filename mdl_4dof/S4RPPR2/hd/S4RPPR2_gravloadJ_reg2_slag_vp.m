% Calculate inertial parameters regressor of gravitation load for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t15 = sin(pkin(6));
t17 = sin(qJ(1));
t26 = t17 * t15;
t18 = cos(qJ(1));
t25 = t18 * t15;
t24 = t18 * pkin(1) + t17 * qJ(2);
t23 = pkin(6) + qJ(4);
t22 = cos(t23);
t21 = sin(t23);
t1 = -t17 * t21 - t18 * t22;
t2 = -t17 * t22 + t18 * t21;
t20 = g(1) * t2 - g(2) * t1;
t19 = g(1) * t1 + g(2) * t2;
t16 = cos(pkin(6));
t12 = t18 * qJ(2);
t10 = t16 * pkin(3) + pkin(2);
t6 = g(1) * t18 + g(2) * t17;
t5 = g(1) * t17 - g(2) * t18;
t4 = t18 * t16 + t26;
t3 = -t17 * t16 + t25;
t7 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (-t17 * pkin(1) + t12) - g(2) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t4, -g(1) * t4 + g(2) * t3, 0, -g(1) * (t12 + (-pkin(1) - pkin(2)) * t17) - g(2) * (t18 * pkin(2) + t24) 0, 0, 0, 0, 0, 0, -t20, t19, 0, -g(1) * (pkin(3) * t25 + t12 + (-pkin(1) - t10) * t17) - g(2) * (pkin(3) * t26 + t18 * t10 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, 0;];
taug_reg  = t7;
