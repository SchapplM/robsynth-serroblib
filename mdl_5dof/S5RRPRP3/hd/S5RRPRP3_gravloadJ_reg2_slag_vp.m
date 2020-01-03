% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = pkin(8) + qJ(4);
t18 = sin(t23);
t19 = cos(t23);
t33 = t19 * pkin(4) + t18 * qJ(5);
t24 = qJ(1) + qJ(2);
t20 = sin(t24);
t21 = cos(t24);
t7 = g(1) * t20 - g(2) * t21;
t8 = g(1) * t21 + g(2) * t20;
t28 = sin(qJ(1));
t41 = t28 * pkin(1);
t27 = -pkin(7) - qJ(3);
t40 = t21 * t27;
t39 = t21 * pkin(2) + t20 * qJ(3);
t26 = cos(pkin(8));
t15 = t26 * pkin(3) + pkin(2);
t10 = t21 * t15;
t37 = t33 * t21 + t10;
t36 = -t20 * pkin(2) + t21 * qJ(3);
t35 = -t20 * t27 + t10;
t29 = cos(qJ(1));
t34 = g(1) * t28 - g(2) * t29;
t31 = -t20 * t15 - t40;
t30 = (-g(1) * (-t15 - t33) + g(2) * t27) * t20;
t22 = t29 * pkin(1);
t6 = t7 * t26;
t5 = t7 * sin(pkin(8));
t4 = t7 * t19;
t3 = t7 * t18;
t2 = g(3) * t18 + t8 * t19;
t1 = -g(3) * t19 + t8 * t18;
t9 = [0, 0, 0, 0, 0, 0, t34, g(1) * t29 + g(2) * t28, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t34 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t36 - t41) - g(2) * (t22 + t39), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t31 - t41) - g(2) * (t22 + t35), 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(1) * (-t40 - t41) - g(2) * (t22 + t37) + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t36 - g(2) * t39, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t31 - g(2) * t35, 0, 0, 0, 0, 0, 0, t4, -t8, t3, g(1) * t40 - g(2) * t37 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t33 + t8 * (pkin(4) * t18 - qJ(5) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
