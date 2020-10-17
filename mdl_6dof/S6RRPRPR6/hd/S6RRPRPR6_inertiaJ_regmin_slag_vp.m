% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:26:12
% EndTime: 2019-05-06 14:26:16
% DurationCPUTime: 1.06s
% Computational Cost: add. (1180->130), mult. (2921->252), div. (0->0), fcn. (3310->10), ass. (0->99)
t60 = sin(pkin(11));
t61 = sin(pkin(6));
t62 = cos(pkin(11));
t66 = sin(qJ(2));
t69 = cos(qJ(2));
t38 = (t60 * t69 + t62 * t66) * t61;
t63 = cos(pkin(6));
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t30 = t38 * t68 + t63 * t65;
t111 = -0.2e1 * t30;
t110 = 0.2e1 * t30;
t109 = 2 * qJ(5);
t108 = pkin(4) + pkin(10);
t107 = pkin(1) * t66;
t100 = t61 * t66;
t99 = t61 * t69;
t37 = t100 * t60 - t62 * t99;
t106 = t37 * pkin(4);
t105 = t68 * pkin(4);
t29 = t38 * t65 - t63 * t68;
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t20 = t29 * t64 + t37 * t67;
t104 = t20 * t67;
t103 = t30 * t64;
t25 = t30 * t65;
t102 = t30 * t68;
t55 = t61 ^ 2;
t101 = t55 * t69;
t98 = t64 * t65;
t97 = t64 * t68;
t96 = t64 * t108;
t19 = -t29 * t67 + t37 * t64;
t95 = t65 * t19;
t94 = t65 * t29;
t93 = t65 * t37;
t52 = pkin(2) * t60 + pkin(9);
t46 = t65 * t52;
t92 = t65 * t68;
t91 = t67 * t64;
t90 = t67 * t68;
t89 = t67 * t108;
t47 = t68 * t52;
t88 = pkin(8) + qJ(3);
t49 = t63 * t69 * pkin(1);
t31 = t63 * pkin(2) - t100 * t88 + t49;
t82 = t63 * t107;
t34 = t88 * t99 + t82;
t17 = t31 * t60 + t34 * t62;
t14 = pkin(9) * t63 + t17;
t45 = (-pkin(2) * t69 - pkin(1)) * t61;
t21 = t37 * pkin(3) - t38 * pkin(9) + t45;
t10 = t14 * t68 + t21 * t65;
t57 = t65 ^ 2;
t59 = t68 ^ 2;
t87 = t57 + t59;
t86 = qJ(5) * t68;
t36 = t37 * qJ(5);
t85 = t65 * qJ(5);
t84 = 0.2e1 * t61 * t63;
t83 = -0.2e1 * t92;
t81 = t30 * t97;
t80 = t37 * t46;
t79 = t30 * t90;
t78 = t37 * t47;
t77 = t36 + t10;
t53 = -t62 * pkin(2) - pkin(3);
t9 = -t14 * t65 + t68 * t21;
t16 = t62 * t31 - t34 * t60;
t7 = -t9 - t106;
t76 = t7 * t65 + t68 * t77;
t75 = -pkin(4) * t65 + t86;
t74 = -t108 * t65 + t86;
t73 = t53 - t85;
t13 = -t63 * pkin(3) - t16;
t72 = -t30 * qJ(5) + t13;
t58 = t67 ^ 2;
t56 = t64 ^ 2;
t54 = t67 * t65;
t44 = pkin(5) * t68 + t47;
t43 = pkin(5) * t65 + t46;
t42 = t73 - t105;
t41 = pkin(8) * t99 + t82;
t40 = -pkin(8) * t100 + t49;
t39 = -t108 * t68 + t73;
t35 = t68 * t37;
t28 = t30 ^ 2;
t26 = t30 * t67;
t23 = t39 * t67 + t43 * t64;
t22 = -t39 * t64 + t43 * t67;
t15 = t20 * t65;
t8 = pkin(4) * t29 + t72;
t5 = t108 * t29 + t72;
t4 = -pkin(5) * t29 + t77;
t3 = t30 * pkin(5) - t108 * t37 - t9;
t2 = t3 * t64 + t5 * t67;
t1 = t3 * t67 - t5 * t64;
t6 = [1, 0, 0, t55 * t66 ^ 2, 0.2e1 * t66 * t101, t66 * t84, t69 * t84, t63 ^ 2, 0.2e1 * pkin(1) * t101 + 0.2e1 * t40 * t63, -0.2e1 * t107 * t55 - 0.2e1 * t41 * t63, -0.2e1 * t16 * t38 - 0.2e1 * t17 * t37, t16 ^ 2 + t17 ^ 2 + t45 ^ 2, t28, t29 * t111, t37 * t110, -0.2e1 * t29 * t37, t37 ^ 2, 0.2e1 * t13 * t29 + 0.2e1 * t37 * t9, -0.2e1 * t10 * t37 + 0.2e1 * t13 * t30, -0.2e1 * t29 * t77 + 0.2e1 * t30 * t7, -0.2e1 * t29 * t8 + 0.2e1 * t37 * t7, -0.2e1 * t30 * t8 + 0.2e1 * t37 * t77, t7 ^ 2 + t77 ^ 2 + t8 ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, t20 * t110, t19 * t111, t28, 0.2e1 * t1 * t30 + 0.2e1 * t19 * t4, -0.2e1 * t2 * t30 + 0.2e1 * t20 * t4; 0, 0, 0, 0, 0, t100, t99, t63, t40, -t41 (-t37 * t60 - t38 * t62) * pkin(2) (t16 * t62 + t17 * t60) * pkin(2), t25, -t94 + t102, t93, t35, 0, -t13 * t68 + t29 * t53 - t80, t13 * t65 + t30 * t53 - t78 (-t29 * t68 + t25) * t52 + t76, -t29 * t42 + t68 * t8 + t80, -t30 * t42 - t65 * t8 + t78, t8 * t42 + t52 * t76, -t20 * t97 (t19 * t64 - t104) * t68, t15 - t81, -t79 - t95, t25, t1 * t65 + t19 * t44 + t22 * t30 + t4 * t90, -t2 * t65 + t20 * t44 - t23 * t30 - t4 * t97; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t60 ^ 2 + t62 ^ 2) * pkin(2) ^ 2, t57, 0.2e1 * t92, 0, 0, 0, -0.2e1 * t53 * t68, 0.2e1 * t53 * t65, 0.2e1 * t87 * t52, 0.2e1 * t42 * t68, -0.2e1 * t42 * t65, t52 ^ 2 * t87 + t42 ^ 2, t56 * t59, 0.2e1 * t59 * t91, t64 * t83, t67 * t83, t57, 0.2e1 * t22 * t65 + 0.2e1 * t44 * t90, -0.2e1 * t23 * t65 - 0.2e1 * t44 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, t35, -t93, -t94 - t102, -t35, t93, t65 * t77 - t68 * t7, 0, 0, 0, 0, 0, -t79 + t95, t15 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, t37, t9, -t10, -pkin(4) * t30 - qJ(5) * t29, -t9 - 0.2e1 * t106, t77 + t36, -pkin(4) * t7 + qJ(5) * t77, t104, -t19 * t67 - t20 * t64, t26, -t103, 0, qJ(5) * t19 - t30 * t89 + t4 * t64, qJ(5) * t20 + t30 * t96 + t4 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t68, 0, -t46, -t47, t75, t46, t47, t75 * t52, -t64 * t90 (t56 - t58) * t68, t54, -t98, 0, t44 * t64 + t67 * t74, t44 * t67 - t64 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t65, 0, -t68, t65, t85 + t105, 0, 0, 0, 0, 0, t98, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t109, pkin(4) ^ 2 + (qJ(5) ^ 2) t58, -0.2e1 * t91, 0, 0, 0, t64 * t109, t67 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t37, 0, t7, 0, 0, 0, 0, 0, t26, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, t46, 0, 0, 0, 0, 0, t54, -t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t30, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t90, t65, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t64, 0, -t89, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t6;
