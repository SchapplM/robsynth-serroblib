% Calculate minimal parameter regressor of gravitation load for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(3));
t61 = qJ(4) * t33 + pkin(2);
t34 = sin(qJ(2));
t30 = sin(pkin(7));
t32 = cos(pkin(7));
t46 = g(1) * t32 + g(2) * t30;
t60 = t46 * t34;
t36 = cos(qJ(2));
t59 = pkin(6) * t36;
t56 = t33 * t34;
t55 = t33 * t36;
t31 = cos(pkin(8));
t54 = t34 * t31;
t35 = cos(qJ(3));
t53 = t34 * t35;
t52 = t35 * t36;
t50 = g(3) * t56;
t13 = t30 * t55 + t32 * t35;
t14 = t30 * t52 - t32 * t33;
t49 = -t13 * pkin(3) + qJ(4) * t14;
t15 = -t30 * t35 + t32 * t55;
t16 = t30 * t33 + t32 * t52;
t48 = -t15 * pkin(3) + qJ(4) * t16;
t47 = pkin(3) * t52 + t34 * pkin(6) + t61 * t36;
t29 = sin(pkin(8));
t45 = -pkin(4) * t31 - qJ(5) * t29;
t43 = t29 * t36 - t31 * t53;
t42 = t29 * t53 + t31 * t36;
t17 = t29 * t52 - t54;
t7 = t42 * t30;
t9 = t42 * t32;
t41 = g(1) * t9 + g(2) * t7 - g(3) * t17;
t40 = g(1) * t15 + g(2) * t13 + t50;
t39 = g(1) * t16 + g(2) * t14 + g(3) * t53;
t38 = -g(3) * t36 + t60;
t37 = (pkin(3) * t35 + t61) * t60;
t21 = t32 * t59;
t20 = qJ(4) * t53;
t19 = t30 * t59;
t18 = t34 * t29 + t31 * t52;
t10 = t43 * t32;
t8 = t43 * t30;
t6 = t38 * t33;
t3 = t40 * t31;
t2 = t40 * t29;
t1 = -g(1) * t10 - g(2) * t8 - g(3) * t18;
t4 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, t38, g(3) * t34 + t46 * t36, 0, 0, 0, 0, 0, t38 * t35, -t6, t1, -t41, t6, -g(1) * t21 - g(2) * t19 - g(3) * t47 + t37, t1, t6, t41, -g(1) * (pkin(4) * t10 - qJ(5) * t9 + t21) - g(2) * (pkin(4) * t8 - qJ(5) * t7 + t19) - g(3) * (pkin(4) * t18 + qJ(5) * t17 + t47) + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, t3, -t2, -t39, -g(1) * t48 - g(2) * t49 - g(3) * (-pkin(3) * t56 + t20), t3, -t39, t2, -g(1) * (t45 * t15 + t48) - g(2) * (t45 * t13 + t49) - g(3) * t20 - (-pkin(3) + t45) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t29 - t32 * t54) - g(2) * (t14 * t29 - t30 * t54) - g(3) * t42;];
taug_reg = t4;
