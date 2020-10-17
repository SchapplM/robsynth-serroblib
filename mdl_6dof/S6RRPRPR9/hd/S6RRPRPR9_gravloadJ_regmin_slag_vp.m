% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:18:28
% EndTime: 2019-05-06 15:18:30
% DurationCPUTime: 0.56s
% Computational Cost: add. (490->116), mult. (829->193), div. (0->0), fcn. (1020->14), ass. (0->58)
t38 = sin(qJ(2));
t39 = sin(qJ(1));
t40 = cos(qJ(2));
t54 = cos(pkin(6));
t66 = cos(qJ(1));
t45 = t54 * t66;
t13 = t39 * t38 - t40 * t45;
t30 = pkin(12) + qJ(6);
t25 = sin(t30);
t27 = cos(t30);
t14 = t38 * t45 + t39 * t40;
t31 = pkin(11) + qJ(4);
t26 = sin(t31);
t28 = cos(t31);
t34 = sin(pkin(6));
t52 = t34 * t66;
t6 = t14 * t28 - t26 * t52;
t72 = -t13 * t27 + t6 * t25;
t71 = t13 * t25 + t6 * t27;
t50 = t39 * t54;
t15 = t66 * t38 + t40 * t50;
t70 = -g(1) * t15 - g(2) * t13;
t67 = g(3) * t34;
t63 = t25 * t28;
t62 = t27 * t28;
t32 = sin(pkin(12));
t61 = t28 * t32;
t35 = cos(pkin(12));
t60 = t28 * t35;
t59 = t28 * t40;
t58 = t34 * t38;
t57 = t34 * t39;
t56 = t34 * t40;
t55 = t66 * pkin(1) + pkin(8) * t57;
t33 = sin(pkin(11));
t53 = t33 * t57;
t51 = -t39 * pkin(1) + pkin(8) * t52;
t49 = t33 * t52;
t5 = t14 * t26 + t28 * t52;
t16 = -t38 * t50 + t66 * t40;
t9 = t16 * t26 - t28 * t57;
t48 = -g(1) * t5 + g(2) * t9;
t47 = g(1) * t13 - g(2) * t15;
t46 = g(1) * t16 + g(2) * t14;
t11 = t26 * t58 - t54 * t28;
t43 = g(1) * t9 + g(2) * t5 + g(3) * t11;
t10 = t16 * t28 + t26 * t57;
t12 = t54 * t26 + t28 * t58;
t42 = g(1) * t10 + g(2) * t6 + g(3) * t12;
t4 = g(3) * t56 + t70;
t41 = g(3) * t58 + t46;
t37 = -pkin(9) - qJ(3);
t36 = cos(pkin(11));
t24 = t36 * pkin(3) + pkin(2);
t3 = t4 * t26;
t2 = t10 * t27 + t15 * t25;
t1 = -t10 * t25 + t15 * t27;
t7 = [0, g(1) * t39 - g(2) * t66, g(1) * t66 + g(2) * t39, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t16, -t47, -g(1) * (-t14 * t36 + t49) - g(2) * (t16 * t36 + t53) -g(1) * (t14 * t33 + t36 * t52) - g(2) * (-t16 * t33 + t36 * t57) t47, -g(1) * (-t14 * pkin(2) - t13 * qJ(3) + t51) - g(2) * (t16 * pkin(2) + t15 * qJ(3) + t55) 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t10, t48, -g(1) * (-t13 * t32 - t35 * t6) - g(2) * (t10 * t35 + t15 * t32) -g(1) * (-t13 * t35 + t32 * t6) - g(2) * (-t10 * t32 + t15 * t35) -t48, -g(1) * (pkin(3) * t49 - pkin(4) * t6 - qJ(5) * t5 + t13 * t37 - t14 * t24 + t51) - g(2) * (pkin(3) * t53 + t10 * pkin(4) + t9 * qJ(5) - t15 * t37 + t16 * t24 + t55) 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t2, -g(1) * t72 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t4, t41, -t4 * t36, t4 * t33, -t41, -g(1) * (-t15 * pkin(2) + t16 * qJ(3)) - g(2) * (-t13 * pkin(2) + t14 * qJ(3)) - (pkin(2) * t40 + qJ(3) * t38) * t67, 0, 0, 0, 0, 0, -t4 * t28, t3, -g(1) * (-t15 * t60 + t16 * t32) - g(2) * (-t13 * t60 + t14 * t32) - (t32 * t38 + t35 * t59) * t67, -g(1) * (t15 * t61 + t16 * t35) - g(2) * (t13 * t61 + t14 * t35) - (-t32 * t59 + t35 * t38) * t67, -t3 (t38 * t67 + t46) * t37 + (-t40 * t67 - t70) * (pkin(4) * t28 + qJ(5) * t26 + t24) 0, 0, 0, 0, 0, -g(1) * (-t15 * t62 + t16 * t25) - g(2) * (-t13 * t62 + t14 * t25) - (t25 * t38 + t27 * t59) * t67, -g(1) * (t15 * t63 + t16 * t27) - g(2) * (t13 * t63 + t14 * t27) - (-t25 * t59 + t27 * t38) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t42, t43 * t35, -t43 * t32, -t42, -g(1) * (-t9 * pkin(4) + t10 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - g(3) * (-t11 * pkin(4) + t12 * qJ(5)) 0, 0, 0, 0, 0, t43 * t27, -t43 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t72 - g(3) * (-t12 * t25 - t27 * t56) g(1) * t2 + g(2) * t71 - g(3) * (-t12 * t27 + t25 * t56);];
taug_reg  = t7;
