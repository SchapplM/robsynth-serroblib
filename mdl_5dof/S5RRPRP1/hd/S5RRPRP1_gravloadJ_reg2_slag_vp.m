% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP1
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t26 = qJ(1) + qJ(2);
t21 = pkin(8) + t26;
t16 = sin(t21);
t17 = cos(t21);
t22 = sin(t26);
t18 = pkin(2) * t22;
t30 = cos(qJ(4));
t20 = t30 * pkin(4) + pkin(3);
t27 = -qJ(5) - pkin(7);
t36 = t16 * t20 + t17 * t27 + t18;
t23 = cos(t26);
t19 = pkin(2) * t23;
t35 = t17 * pkin(3) + t16 * pkin(7) + t19;
t34 = t16 * pkin(3) - t17 * pkin(7) + t18;
t33 = -t16 * t27 + t17 * t20 + t19;
t6 = g(2) * t17 + g(3) * t16;
t5 = g(2) * t16 - g(3) * t17;
t8 = -g(2) * t23 - g(3) * t22;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t32 = -g(2) * t31 - g(3) * t29;
t28 = sin(qJ(4));
t2 = -g(1) * t30 + t5 * t28;
t25 = t31 * pkin(1);
t24 = t29 * pkin(1);
t7 = g(2) * t22 - g(3) * t23;
t4 = t6 * t30;
t3 = t6 * t28;
t1 = g(1) * t28 + t5 * t30;
t9 = [0, 0, 0, 0, 0, 0, t32, g(2) * t29 - g(3) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, -t6, t5, 0, -g(2) * (t19 + t25) - g(3) * (t18 + t24), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t25 + t35) - g(3) * (t24 + t34), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t25 + t33) - g(3) * (t24 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t35 - g(3) * t34, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t33 - g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6;];
taug_reg = t9;
