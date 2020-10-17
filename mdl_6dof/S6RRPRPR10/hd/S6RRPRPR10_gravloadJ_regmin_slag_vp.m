% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:43:25
% EndTime: 2019-05-06 15:43:27
% DurationCPUTime: 0.48s
% Computational Cost: add. (412->101), mult. (750->158), div. (0->0), fcn. (916->12), ass. (0->54)
t38 = sin(qJ(2));
t39 = sin(qJ(1));
t41 = cos(qJ(2));
t55 = cos(pkin(6));
t65 = cos(qJ(1));
t46 = t55 * t65;
t17 = t39 * t38 - t41 * t46;
t37 = sin(qJ(6));
t40 = cos(qJ(6));
t18 = t38 * t46 + t39 * t41;
t32 = pkin(11) + qJ(4);
t29 = sin(t32);
t30 = cos(t32);
t34 = sin(pkin(6));
t53 = t34 * t65;
t9 = t18 * t29 + t30 * t53;
t70 = t17 * t40 + t9 * t37;
t69 = -t17 * t37 + t9 * t40;
t51 = t39 * t55;
t19 = t65 * t38 + t41 * t51;
t66 = g(3) * t34;
t6 = -g(1) * t19 - g(2) * t17 + t41 * t66;
t62 = t29 * t37;
t61 = t29 * t40;
t60 = t34 * t38;
t59 = t34 * t39;
t58 = t37 * t41;
t57 = t40 * t41;
t56 = t65 * pkin(1) + pkin(8) * t59;
t33 = sin(pkin(11));
t54 = t33 * t59;
t52 = -t39 * pkin(1) + pkin(8) * t53;
t10 = t18 * t30 - t29 * t53;
t50 = t33 * t53;
t20 = -t38 * t51 + t65 * t41;
t13 = t20 * t29 - t30 * t59;
t49 = -g(1) * t9 + g(2) * t13;
t14 = t20 * t30 + t29 * t59;
t48 = -g(1) * t10 + g(2) * t14;
t8 = g(1) * t17 - g(2) * t19;
t47 = g(1) * t20 + g(2) * t18;
t15 = t29 * t60 - t55 * t30;
t44 = g(1) * t13 + g(2) * t9 + g(3) * t15;
t16 = t55 * t29 + t30 * t60;
t43 = g(1) * t14 + g(2) * t10 + g(3) * t16;
t42 = g(3) * t60 + t47;
t36 = -pkin(9) - qJ(3);
t35 = cos(pkin(11));
t28 = t35 * pkin(3) + pkin(2);
t5 = t13 * t37 + t19 * t40;
t4 = t13 * t40 - t19 * t37;
t3 = t6 * t30;
t2 = t6 * t29;
t1 = [0, g(1) * t39 - g(2) * t65, g(1) * t65 + g(2) * t39, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t20, -t8, -g(1) * (-t18 * t35 + t50) - g(2) * (t20 * t35 + t54) -g(1) * (t18 * t33 + t35 * t53) - g(2) * (-t20 * t33 + t35 * t59) t8, -g(1) * (-t18 * pkin(2) - t17 * qJ(3) + t52) - g(2) * (t20 * pkin(2) + t19 * qJ(3) + t56) 0, 0, 0, 0, 0, -t48, t49, t8, t48, -t49, -g(1) * (pkin(3) * t50 - pkin(4) * t10 - qJ(5) * t9 + t17 * t36 - t18 * t28 + t52) - g(2) * (pkin(3) * t54 + t14 * pkin(4) + t13 * qJ(5) - t19 * t36 + t20 * t28 + t56) 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t5, g(1) * t69 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, -t6, t42, -t6 * t35, t6 * t33, -t42, -g(1) * (-t19 * pkin(2) + t20 * qJ(3)) - g(2) * (-t17 * pkin(2) + t18 * qJ(3)) - (pkin(2) * t41 + qJ(3) * t38) * t66, 0, 0, 0, 0, 0, -t3, t2, -t42, t3, -t2 (t38 * t66 + t47) * t36 - t6 * (pkin(4) * t30 + qJ(5) * t29 + t28) 0, 0, 0, 0, 0, -g(1) * (-t19 * t62 + t20 * t40) - g(2) * (-t17 * t62 + t18 * t40) - (t29 * t58 + t38 * t40) * t66, -g(1) * (-t19 * t61 - t20 * t37) - g(2) * (-t17 * t61 - t18 * t37) - (t29 * t57 - t37 * t38) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43, 0, -t44, -t43, -g(1) * (-t13 * pkin(4) + t14 * qJ(5)) - g(2) * (-t9 * pkin(4) + t10 * qJ(5)) - g(3) * (-t15 * pkin(4) + t16 * qJ(5)) 0, 0, 0, 0, 0, -t43 * t37, -t43 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t69 - g(3) * (t15 * t40 + t34 * t58) g(1) * t5 + g(2) * t70 - g(3) * (-t15 * t37 + t34 * t57);];
taug_reg  = t1;
