% Calculate inertial parameters regressor of gravitation load for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t11 = g(1) * t28 + g(2) * t25;
t24 = sin(qJ(2));
t45 = t11 * t24;
t17 = t24 * qJ(3);
t27 = cos(qJ(2));
t37 = t27 * pkin(2) + t17;
t43 = g(1) * t25;
t40 = t27 * pkin(3);
t39 = t24 * t28;
t38 = t27 * t28;
t36 = t28 * pkin(1) + t25 * pkin(5);
t35 = qJ(3) * t27;
t34 = pkin(2) * t38 + t28 * t17 + t36;
t33 = -g(2) * t28 + t43;
t23 = sin(qJ(4));
t26 = cos(qJ(4));
t32 = t27 * t23 - t24 * t26;
t9 = t24 * t23 + t27 * t26;
t31 = -pkin(1) - t37;
t3 = t32 * t25;
t5 = t23 * t38 - t26 * t39;
t30 = g(1) * t5 + g(2) * t3 + g(3) * t9;
t4 = t9 * t25;
t6 = t9 * t28;
t29 = g(1) * t6 + g(2) * t4 - g(3) * t32;
t20 = t28 * pkin(5);
t15 = t28 * t35;
t13 = t25 * t35;
t8 = t33 * t27;
t7 = t33 * t24;
t2 = g(3) * t24 + t11 * t27;
t1 = -g(3) * t27 + t45;
t10 = [0, 0, 0, 0, 0, 0, t33, t11, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * (-t25 * pkin(1) + t20) - g(2) * t36, 0, 0, 0, 0, 0, 0, t8, -t11, t7, -g(1) * t20 - g(2) * t34 - t31 * t43, 0, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5, t11, -g(1) * (-t28 * pkin(6) + t20) - g(2) * (pkin(3) * t38 + t34) + (-g(1) * (t31 - t40) + g(2) * pkin(6)) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(2) * t39 + t15) - g(2) * (-t25 * t24 * pkin(2) + t13) - g(3) * t37, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, -g(1) * t15 - g(2) * t13 - g(3) * (t37 + t40) + (pkin(2) + pkin(3)) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t29, 0, 0;];
taug_reg = t10;
