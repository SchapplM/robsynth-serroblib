% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:46:33
% EndTime: 2019-05-05 02:46:35
% DurationCPUTime: 0.38s
% Computational Cost: add. (276->81), mult. (520->139), div. (0->0), fcn. (623->12), ass. (0->53)
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t50 = cos(pkin(6));
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t56 = t28 * t32;
t65 = -t31 * t56 + t50 * t34;
t35 = cos(qJ(2));
t27 = sin(pkin(10));
t46 = t27 * t50;
t49 = cos(pkin(10));
t17 = -t32 * t46 + t49 * t35;
t55 = t28 * t34;
t64 = -t17 * t31 + t27 * t55;
t63 = g(3) * t28;
t40 = t50 * t49;
t14 = t27 * t32 - t35 * t40;
t15 = t27 * t35 + t32 * t40;
t23 = t34 * pkin(3) + pkin(2);
t29 = -qJ(4) - pkin(8);
t62 = -t14 * t23 - t15 * t29;
t16 = t49 * t32 + t35 * t46;
t61 = -t16 * t23 - t17 * t29;
t26 = qJ(3) + pkin(11);
t24 = sin(t26);
t30 = sin(qJ(6));
t59 = t24 * t30;
t33 = cos(qJ(6));
t58 = t24 * t33;
t57 = t27 * t28;
t54 = t28 * t35;
t53 = t29 * t32;
t52 = t30 * t35;
t51 = t33 * t35;
t45 = t28 * t49;
t43 = t64 * pkin(3);
t25 = cos(t26);
t42 = pkin(4) * t25 + qJ(5) * t24;
t41 = t65 * pkin(3);
t11 = t50 * t24 + t25 * t56;
t5 = t15 * t25 - t24 * t45;
t7 = t17 * t25 + t24 * t57;
t39 = -g(1) * t7 - g(2) * t5 - g(3) * t11;
t38 = -t15 * t31 - t34 * t45;
t2 = -g(1) * t16 - g(2) * t14 + g(3) * t54;
t37 = g(1) * t17 + g(2) * t15 + g(3) * t56;
t36 = t38 * pkin(3);
t18 = t23 * t54;
t10 = t24 * t56 - t50 * t25;
t6 = t17 * t24 - t25 * t57;
t4 = t15 * t24 + t25 * t45;
t1 = -g(1) * t6 - g(2) * t4 - g(3) * t10;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t2, t37, 0, 0, 0, 0, 0, -t2 * t34, t2 * t31, -t37, -g(1) * t61 - g(2) * t62 - g(3) * (-t28 * t53 + t18) -t37, t2 * t25, -t2 * t24, -g(1) * (-t42 * t16 + t61) - g(2) * (-t42 * t14 + t62) - g(3) * t18 - (t42 * t35 - t53) * t63, 0, 0, 0, 0, 0, -g(1) * (-t16 * t59 + t17 * t33) - g(2) * (-t14 * t59 + t15 * t33) - (t24 * t52 + t32 * t33) * t63, -g(1) * (-t16 * t58 - t17 * t30) - g(2) * (-t14 * t58 - t15 * t30) - (t24 * t51 - t30 * t32) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t64 - g(2) * t38 - g(3) * t65, -g(1) * (-t17 * t34 - t31 * t57) - g(2) * (-t15 * t34 + t31 * t45) - g(3) * (-t50 * t31 - t32 * t55) 0, -g(1) * t43 - g(2) * t36 - g(3) * t41, 0, t1, t39, -g(1) * (-t6 * pkin(4) + t7 * qJ(5) + t43) - g(2) * (-t4 * pkin(4) + t5 * qJ(5) + t36) - g(3) * (-t10 * pkin(4) + t11 * qJ(5) + t41) 0, 0, 0, 0, 0, t39 * t30, t39 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t30 + t6 * t33) - g(2) * (-t14 * t30 + t4 * t33) - g(3) * (t10 * t33 + t28 * t52) -g(1) * (-t16 * t33 - t6 * t30) - g(2) * (-t14 * t33 - t4 * t30) - g(3) * (-t10 * t30 + t28 * t51);];
taug_reg  = t3;
