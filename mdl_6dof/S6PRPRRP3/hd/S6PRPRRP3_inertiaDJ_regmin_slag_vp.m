% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:45
% EndTime: 2019-03-08 20:07:48
% DurationCPUTime: 0.95s
% Computational Cost: add. (1307->164), mult. (3331->300), div. (0->0), fcn. (3365->10), ass. (0->99)
t57 = cos(pkin(11));
t63 = cos(qJ(4));
t106 = t63 * t57;
t55 = sin(pkin(11));
t60 = sin(qJ(4));
t40 = t60 * t55 - t106;
t41 = t63 * t55 + t60 * t57;
t49 = -t57 * pkin(3) - pkin(2);
t27 = t40 * pkin(4) - t41 * pkin(9) + t49;
t105 = pkin(8) + qJ(3);
t43 = t105 * t55;
t44 = t105 * t57;
t30 = -t60 * t43 + t63 * t44;
t62 = cos(qJ(5));
t28 = t62 * t30;
t59 = sin(qJ(5));
t103 = t59 * t27 + t28;
t56 = sin(pkin(6));
t61 = sin(qJ(2));
t112 = t56 * t61;
t58 = cos(pkin(6));
t33 = -t55 * t112 + t58 * t57;
t34 = t57 * t112 + t58 * t55;
t72 = t63 * t33 - t60 * t34;
t65 = t72 * qJD(4);
t53 = t59 ^ 2;
t54 = t62 ^ 2;
t101 = t53 - t54;
t82 = t101 * qJD(5);
t64 = cos(qJ(2));
t111 = t56 * t64;
t19 = t60 * t33 + t63 * t34;
t13 = -t59 * t111 + t62 * t19;
t66 = t40 * t64;
t70 = t62 * t111 + t59 * t19;
t98 = qJD(2) * t61;
t90 = t56 * t98;
t99 = qJD(2) * t56;
t5 = -t62 * (-t66 * t99 + t65) - t59 * t90 + t70 * qJD(5);
t96 = qJD(5) * t62;
t97 = qJD(5) * t59;
t6 = -t19 * t96 - t59 * t65 + (t64 * t97 + (t59 * t66 + t61 * t62) * qJD(2)) * t56;
t118 = (-t13 * t62 - t59 * t70) * qJD(5) + t5 * t59 - t6 * t62;
t100 = qJ(6) * t41;
t36 = t41 * qJD(4);
t35 = t40 * qJD(4);
t71 = qJ(6) * t35 - qJD(6) * t41;
t107 = t63 * t43;
t15 = qJD(4) * t107 - qJD(3) * t106 + (qJD(3) * t55 + qJD(4) * t44) * t60;
t26 = t36 * pkin(4) + t35 * pkin(9);
t85 = t59 * t15 + t62 * t26;
t1 = t36 * pkin(5) + t71 * t62 + (-t28 + (-t27 + t100) * t59) * qJD(5) + t85;
t92 = t41 * t96;
t94 = -t62 * t15 + t59 * t26 + t27 * t96;
t2 = -qJ(6) * t92 + (-qJD(5) * t30 + t71) * t59 + t94;
t84 = t62 * t27 - t59 * t30;
t7 = t40 * pkin(5) - t62 * t100 + t84;
t8 = -t59 * t100 + t103;
t117 = -t1 * t62 - t2 * t59 + (t59 * t7 - t62 * t8) * qJD(5);
t116 = 0.2e1 * t49;
t115 = t41 * t35;
t114 = t41 * t59;
t113 = t41 * t62;
t110 = t59 * t36;
t109 = t62 * t35;
t108 = t62 * t36;
t104 = -qJ(6) - pkin(9);
t102 = t55 ^ 2 + t57 ^ 2;
t95 = -0.2e1 * pkin(4) * qJD(5);
t93 = pkin(5) * t97;
t91 = t72 * t97;
t89 = t64 * t99;
t88 = t59 * t96;
t87 = t102 * t64;
t86 = -0.4e1 * t59 * t113;
t29 = t60 * t44 + t107;
t83 = qJD(5) * t104;
t81 = 0.2e1 * t102 * qJD(3);
t80 = pkin(4) * t35 - pkin(9) * t36;
t79 = pkin(4) * t41 + pkin(9) * t40;
t75 = -t13 * t59 + t62 * t70;
t73 = -t33 * t55 + t34 * t57;
t69 = -t59 * t35 + t92;
t68 = t41 * t97 + t109;
t67 = t40 * t96 + t110;
t16 = t41 * qJD(3) + t30 * qJD(4);
t50 = -t62 * pkin(5) - pkin(4);
t46 = t104 * t62;
t45 = t104 * t59;
t38 = t41 ^ 2;
t32 = -t59 * qJD(6) + t62 * t83;
t31 = t62 * qJD(6) + t59 * t83;
t25 = -t40 * t97 + t108;
t17 = pkin(5) * t114 + t29;
t11 = t19 * qJD(4) + t41 * t89;
t9 = t69 * pkin(5) + t16;
t4 = -t103 * qJD(5) + t85;
t3 = t30 * t97 - t94;
t10 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t73 - t112) * t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t11 * t72 - 0.2e1 * t13 * t5 - 0.2e1 * t6 * t70; 0, 0, -t90, -t89, -t57 * t90, t55 * t90, t87 * t99, t73 * qJD(3) + (-pkin(2) * t61 + qJ(3) * t87) * t99, 0, 0, 0, 0, 0 (-t36 * t64 + t40 * t98) * t56 (t35 * t64 + t41 * t98) * t56, 0, 0, 0, 0, 0, t11 * t114 - t36 * t70 + t6 * t40 - t69 * t72, t11 * t113 - t13 * t36 + t5 * t40 + t68 * t72, t118 * t41 - t75 * t35, -t1 * t70 + t11 * t17 + t13 * t2 - t5 * t8 + t6 * t7 - t72 * t9; 0, 0, 0, 0, 0, 0, t81, qJ(3) * t81, -0.2e1 * t115, 0.2e1 * t35 * t40 - 0.2e1 * t41 * t36, 0, 0, 0, t36 * t116, -t35 * t116, -0.2e1 * t54 * t115 - 0.2e1 * t38 * t88, -t35 * t86 + 0.2e1 * t38 * t82, 0.2e1 * t41 * t108 - 0.2e1 * t68 * t40, -0.2e1 * t41 * t110 - 0.2e1 * t69 * t40, 0.2e1 * t40 * t36, 0.2e1 * t16 * t114 + 0.2e1 * t69 * t29 + 0.2e1 * t84 * t36 + 0.2e1 * t4 * t40, -0.2e1 * t103 * t36 + 0.2e1 * t16 * t113 - 0.2e1 * t68 * t29 + 0.2e1 * t3 * t40, -0.2e1 * (-t59 * t8 - t62 * t7) * t35 + 0.2e1 * t117 * t41, 0.2e1 * t7 * t1 + 0.2e1 * t17 * t9 + 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, 0, 0, 0, 0, t25, -t67 -(-t53 - t54) * t35, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t40 * t89 - t65, 0, 0, 0, 0, 0, -t11 * t62 - t91, t11 * t59 - t72 * t96, t75 * qJD(5) - t5 * t62 - t6 * t59, -pkin(5) * t91 + t11 * t50 + t13 * t31 - t32 * t70 + t6 * t45 + t5 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t16, t15, -t59 * t109 - t41 * t82, qJD(5) * t86 + t101 * t35, t67, t25, 0, -t16 * t62 + t80 * t59 + (t29 * t59 - t79 * t62) * qJD(5), t16 * t59 + t80 * t62 + (t29 * t62 + t79 * t59) * qJD(5) (-t32 * t41 + t35 * t45 + t2 + (t41 * t46 - t7) * qJD(5)) * t62 + (-t31 * t41 - t35 * t46 - t1 + (t41 * t45 - t8) * qJD(5)) * t59, t1 * t45 + t17 * t93 - t2 * t46 + t8 * t31 + t7 * t32 + t9 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t31 + t62 * t32 + (-t45 * t59 - t46 * t62) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t88, -0.2e1 * t82, 0, 0, 0, t59 * t95, t62 * t95, 0.2e1 * t31 * t62 - 0.2e1 * t32 * t59 + 0.2e1 * (-t45 * t62 + t46 * t59) * qJD(5), -0.2e1 * t46 * t31 + 0.2e1 * t45 * t32 + 0.2e1 * t50 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t6 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t69, t36, t4, t3, t68 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t97, 0, -pkin(9) * t96, pkin(9) * t97, -pkin(5) * t96, t32 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
