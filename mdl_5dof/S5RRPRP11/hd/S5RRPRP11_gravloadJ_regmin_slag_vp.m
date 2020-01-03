% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(2));
t19 = t26 * qJ(3);
t29 = cos(qJ(2));
t51 = t29 * pkin(2) + t19;
t27 = sin(qJ(1));
t30 = cos(qJ(1));
t13 = g(1) * t30 + g(2) * t27;
t6 = g(3) * t26 + t13 * t29;
t49 = pkin(2) * t26;
t48 = g(1) * t27;
t44 = g(3) * t29;
t25 = sin(qJ(4));
t43 = t27 * t25;
t28 = cos(qJ(4));
t42 = t27 * t28;
t41 = t29 * t30;
t40 = t30 * t25;
t39 = t30 * t28;
t38 = qJ(3) * t29;
t37 = g(3) * t51;
t36 = pkin(2) * t41 + t27 * pkin(6) + (pkin(1) + t19) * t30;
t7 = -t26 * t39 + t43;
t9 = t26 * t42 + t40;
t35 = g(1) * t9 + g(2) * t7;
t34 = -g(2) * t30 + t48;
t33 = pkin(4) * t25 - qJ(5) * t28;
t32 = -pkin(1) - t51;
t1 = g(1) * t7 - g(2) * t9 + t28 * t44;
t10 = -t26 * t43 + t39;
t8 = t26 * t40 + t42;
t31 = -g(1) * t8 + g(2) * t10 + t25 * t44;
t22 = t30 * pkin(6);
t17 = t30 * t38;
t15 = t27 * t38;
t12 = t34 * t29;
t11 = t34 * t26;
t5 = t13 * t26 - t44;
t4 = t6 * t28;
t3 = t6 * t25;
t2 = -g(1) * t10 - g(2) * t8;
t14 = [0, t34, t13, 0, 0, 0, 0, 0, t12, -t11, -t13, -t12, t11, -g(1) * t22 - g(2) * t36 - t32 * t48, 0, 0, 0, 0, 0, t2, t35, t2, t12, -t35, -g(1) * (t30 * pkin(3) + t10 * pkin(4) + t9 * qJ(5) + t22) - g(2) * (t8 * pkin(4) + pkin(7) * t41 + t7 * qJ(5) + t36) + (-g(1) * (-t29 * pkin(7) + t32) - g(2) * pkin(3)) * t27; 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, -t5, -t6, -g(1) * (-t30 * t49 + t17) - g(2) * (-t27 * t49 + t15) - t37, 0, 0, 0, 0, 0, -t3, -t4, -t3, t5, t4, -g(1) * t17 - g(2) * t15 - t37 + (-g(3) * pkin(7) - t13 * t33) * t29 + (-g(3) * t33 + t13 * (pkin(2) + pkin(7))) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t31, t1, 0, t31, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (t9 * pkin(4) - t10 * qJ(5)) - (-pkin(4) * t28 - qJ(5) * t25) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
