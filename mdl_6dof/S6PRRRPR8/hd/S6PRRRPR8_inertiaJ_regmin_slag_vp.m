% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:10:38
% EndTime: 2019-05-05 09:10:41
% DurationCPUTime: 1.03s
% Computational Cost: add. (759->142), mult. (1986->275), div. (0->0), fcn. (2335->12), ass. (0->102)
t56 = cos(pkin(7));
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t54 = sin(pkin(7));
t60 = sin(qJ(3));
t99 = t54 * t60;
t34 = t59 * t56 + t63 * t99;
t111 = -0.2e1 * t34;
t110 = 0.2e1 * t34;
t109 = -0.2e1 * t59;
t108 = 0.2e1 * t63;
t107 = 2 * qJ(5);
t106 = pkin(4) + pkin(11);
t105 = pkin(2) * t60;
t64 = cos(qJ(3));
t104 = pkin(2) * t64;
t55 = sin(pkin(6));
t57 = cos(pkin(6));
t65 = cos(qJ(2));
t94 = t56 * t65;
t61 = sin(qJ(2));
t97 = t55 * t61;
t98 = t54 * t64;
t18 = -t55 * t64 * t94 - t57 * t98 + t60 * t97;
t103 = t18 * t59;
t102 = t18 * t63;
t33 = -t63 * t56 + t59 * t99;
t58 = sin(qJ(6));
t62 = cos(qJ(6));
t21 = t58 * t33 - t62 * t98;
t101 = t21 * t62;
t30 = t34 * t59;
t49 = t54 ^ 2;
t100 = t49 * t64;
t96 = t55 * t65;
t95 = t56 * t60;
t93 = t58 * t34;
t92 = t58 * t59;
t91 = t58 * t63;
t90 = t58 * t106;
t89 = t59 * t63;
t88 = t62 * t58;
t87 = t62 * t63;
t86 = t62 * t106;
t81 = pkin(9) * t98;
t27 = t81 + (pkin(10) + t105) * t56;
t28 = (-pkin(3) * t64 - pkin(10) * t60 - pkin(2)) * t54;
t15 = t63 * t27 + t59 * t28;
t51 = t59 ^ 2;
t53 = t63 ^ 2;
t85 = t51 + t53;
t84 = qJ(5) * t63;
t83 = 0.2e1 * t98;
t82 = -0.2e1 * t89;
t80 = t59 * t98;
t79 = t63 * t98;
t78 = qJ(5) * t98;
t77 = -t59 * qJ(5) - pkin(3);
t14 = -t59 * t27 + t63 * t28;
t76 = pkin(10) * t80;
t75 = pkin(10) * t79;
t44 = pkin(4) * t98;
t13 = -t14 + t44;
t19 = t57 * t99 + (t60 * t94 + t61 * t64) * t55;
t31 = -t54 * t96 + t57 * t56;
t10 = t19 * t63 + t31 * t59;
t9 = t19 * t59 - t31 * t63;
t74 = t10 * t63 + t9 * t59;
t73 = -pkin(4) * t59 + t84;
t12 = t78 - t15;
t72 = -t12 * t63 + t13 * t59;
t71 = -t106 * t59 + t84;
t70 = t18 * t33 + t9 * t98;
t69 = t10 * t98 + t18 * t34;
t43 = pkin(9) * t99;
t26 = t43 + (-pkin(3) - t104) * t56;
t68 = -t34 * qJ(5) + t26;
t52 = t62 ^ 2;
t50 = t58 ^ 2;
t48 = t63 * pkin(10);
t47 = t59 * pkin(10);
t46 = t62 * t59;
t42 = t63 * pkin(5) + t48;
t41 = t59 * pkin(5) + t47;
t39 = -t63 * pkin(4) + t77;
t37 = -t106 * t63 + t77;
t36 = pkin(2) * t95 + t81;
t35 = t56 * t104 - t43;
t32 = t34 ^ 2;
t29 = t62 * t34;
t20 = t62 * t33 + t58 * t98;
t17 = t62 * t37 + t58 * t41;
t16 = -t58 * t37 + t62 * t41;
t8 = t33 * pkin(4) + t68;
t7 = -t33 * pkin(5) - t12;
t6 = t106 * t33 + t68;
t5 = t34 * pkin(5) + pkin(11) * t98 + t13;
t4 = t18 * t62 + t9 * t58;
t3 = -t18 * t58 + t9 * t62;
t2 = t58 * t5 + t62 * t6;
t1 = t62 * t5 - t58 * t6;
t11 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t18 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t96, -t97, 0, 0, 0, 0, 0, -t18 * t56 - t31 * t98, -t19 * t56 + t31 * t99, 0, 0, 0, 0, 0, t70, t69, -t10 * t33 + t9 * t34, -t70, -t69, -t10 * t12 + t9 * t13 + t18 * t8, 0, 0, 0, 0, 0, -t10 * t20 + t3 * t34, t10 * t21 - t4 * t34; 0, 1, 0, 0, t49 * t60 ^ 2, 0.2e1 * t60 * t100, 0.2e1 * t54 * t95, t56 * t83, t56 ^ 2, 0.2e1 * pkin(2) * t100 + 0.2e1 * t35 * t56, -0.2e1 * t49 * t105 - 0.2e1 * t36 * t56, t32, t33 * t111, t98 * t111, t33 * t83, t49 * t64 ^ 2, -0.2e1 * t14 * t98 + 0.2e1 * t26 * t33, 0.2e1 * t15 * t98 + 0.2e1 * t26 * t34, 0.2e1 * t12 * t33 + 0.2e1 * t13 * t34, -0.2e1 * t13 * t98 - 0.2e1 * t8 * t33, 0.2e1 * t12 * t98 - 0.2e1 * t8 * t34, t12 ^ 2 + t13 ^ 2 + t8 ^ 2, t21 ^ 2, 0.2e1 * t21 * t20, t21 * t110, t20 * t110, t32, 0.2e1 * t1 * t34 - 0.2e1 * t7 * t20, -0.2e1 * t2 * t34 + 0.2e1 * t7 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, 0, 0, 0, 0, -t102, t103, t74, t102, -t103, t74 * pkin(10) + t18 * t39, 0, 0, 0, 0, 0, t10 * t87 + t3 * t59, -t10 * t91 - t4 * t59; 0, 0, 0, 0, 0, 0, t99, t98, t56, t35, -t36, t30, -t59 * t33 + t34 * t63, -t80, -t79, 0, -pkin(3) * t33 - t26 * t63 + t76, -pkin(3) * t34 + t26 * t59 + t75 (-t33 * t63 + t30) * pkin(10) + t72, -t39 * t33 + t8 * t63 - t76, -t39 * t34 - t8 * t59 - t75, t72 * pkin(10) + t8 * t39, -t21 * t91 (-t20 * t58 - t101) * t63, t21 * t59 - t34 * t91, t20 * t59 - t34 * t87, t30, t1 * t59 + t16 * t34 - t42 * t20 + t7 * t87, -t17 * t34 - t2 * t59 + t42 * t21 - t7 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t51, 0.2e1 * t89, 0, 0, 0, pkin(3) * t108, pkin(3) * t109, 0.2e1 * t85 * pkin(10), t39 * t108, t39 * t109, t85 * pkin(10) ^ 2 + t39 ^ 2, t50 * t53, 0.2e1 * t53 * t88, t58 * t82, t62 * t82, t51, 0.2e1 * t16 * t59 + 0.2e1 * t42 * t87, -0.2e1 * t17 * t59 - 0.2e1 * t42 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, t9, t10, -t9 * pkin(4) + t10 * qJ(5), 0, 0, 0, 0, 0, t10 * t58, t10 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, -t98, t14, -t15, -pkin(4) * t34 - qJ(5) * t33, -t14 + 0.2e1 * t44, -0.2e1 * t78 + t15, -t13 * pkin(4) - t12 * qJ(5), t101, t62 * t20 - t21 * t58, t29, -t93, 0, -qJ(5) * t20 - t34 * t86 + t7 * t58, qJ(5) * t21 + t34 * t90 + t7 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t63, 0, -t47, -t48, t73, t47, t48, t73 * pkin(10), -t58 * t87 (t50 - t52) * t63, t46, -t92, 0, t42 * t58 + t71 * t62, t42 * t62 - t71 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t107, pkin(4) ^ 2 + (qJ(5) ^ 2) t52, -0.2e1 * t88, 0, 0, 0, t58 * t107, t62 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t98, 0, t13, 0, 0, 0, 0, 0, t29, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, t47, 0, 0, 0, 0, 0, t46, -t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20, t34, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t87, t59, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t58, 0, -t86, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
