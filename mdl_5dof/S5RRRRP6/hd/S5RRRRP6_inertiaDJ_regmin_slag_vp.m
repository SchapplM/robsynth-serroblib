% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:37
% DurationCPUTime: 0.97s
% Computational Cost: add. (1411->172), mult. (3254->296), div. (0->0), fcn. (2831->6), ass. (0->105)
t71 = sin(qJ(3));
t72 = sin(qJ(2));
t74 = cos(qJ(3));
t75 = cos(qJ(2));
t46 = t71 * t72 - t74 * t75;
t47 = t71 * t75 + t74 * t72;
t63 = -t75 * pkin(2) - pkin(1);
t24 = t46 * pkin(3) - t47 * pkin(8) + t63;
t122 = -pkin(7) - pkin(6);
t55 = t122 * t72;
t56 = t122 * t75;
t34 = t71 * t55 - t74 * t56;
t73 = cos(qJ(4));
t32 = t73 * t34;
t70 = sin(qJ(4));
t112 = t70 * t24 + t32;
t69 = t73 ^ 2;
t110 = t70 ^ 2 - t69;
t87 = t110 * qJD(4);
t123 = qJD(2) + qJD(3);
t121 = t74 * pkin(2);
t30 = t123 * t46;
t120 = t47 * t30;
t119 = t47 * t70;
t118 = t47 * t73;
t31 = t123 * t47;
t117 = t70 * t31;
t116 = t73 * t30;
t115 = t73 * t31;
t114 = -qJ(5) - pkin(8);
t92 = qJD(2) * t122;
t83 = t75 * t92;
t84 = t72 * t92;
t17 = t34 * qJD(3) + t71 * t84 - t74 * t83;
t33 = -t74 * t55 - t71 * t56;
t66 = qJD(4) * t73;
t113 = t17 * t70 + t33 * t66;
t61 = -pkin(3) - t121;
t107 = qJD(3) * t71;
t97 = pkin(2) * t107;
t111 = t61 * t66 + t70 * t97;
t109 = qJ(5) * t47;
t67 = t73 * qJ(5);
t60 = t71 * pkin(2) + pkin(8);
t108 = -qJ(5) - t60;
t106 = qJD(3) * t74;
t105 = qJD(4) * t70;
t104 = t72 * qJD(2);
t103 = t75 * qJD(2);
t102 = -0.2e1 * pkin(1) * qJD(2);
t98 = pkin(2) * t104;
t13 = t31 * pkin(3) + t30 * pkin(8) + t98;
t16 = -t55 * t106 - t56 * t107 - t71 * t83 - t74 * t84;
t101 = t70 * t13 - t73 * t16 + t24 * t66;
t100 = pkin(3) * t105;
t99 = pkin(3) * t66;
t96 = pkin(2) * t106;
t64 = pkin(4) * t105;
t95 = pkin(4) * t66;
t94 = t47 * t66;
t93 = t70 * t66;
t62 = -t73 * pkin(4) - pkin(3);
t91 = -0.4e1 * t70 * t118;
t90 = t73 * t13 + t70 * t16;
t89 = t73 * t24 - t70 * t34;
t88 = qJD(4) * t114;
t86 = qJD(4) * t108;
t85 = t73 * t96;
t82 = t46 * t60 - t47 * t61;
t81 = qJ(5) * t30 - qJD(5) * t47;
t80 = t61 * t105 - t73 * t97;
t79 = -t70 * t30 + t94;
t78 = t47 * t105 + t116;
t77 = t46 * t105 - t115;
t76 = -t30 * t61 - t31 * t60 + (-t46 * t74 + t47 * t71) * qJD(3) * pkin(2);
t65 = t73 * qJD(5);
t58 = 0.2e1 * t93;
t54 = t73 * pkin(8) + t67;
t53 = t114 * t70;
t52 = t62 - t121;
t49 = t64 + t97;
t45 = -0.2e1 * t87;
t44 = t47 ^ 2;
t42 = t73 * t60 + t67;
t41 = t108 * t70;
t37 = -t70 * qJD(5) + t73 * t88;
t36 = t70 * t88 + t65;
t35 = t36 * t73;
t29 = (-qJD(5) - t96) * t70 + t73 * t86;
t28 = t70 * t86 + t65 + t85;
t26 = t33 * t105;
t25 = t28 * t73;
t20 = pkin(4) * t119 + t33;
t19 = t46 * t66 + t117;
t12 = -t70 * t116 - t47 * t87;
t9 = -t70 * t109 + t112;
t8 = qJD(4) * t91 + t110 * t30;
t7 = t46 * pkin(4) - t47 * t67 + t89;
t6 = t79 * pkin(4) + t17;
t5 = -t112 * qJD(4) + t90;
t4 = t34 * t105 - t101;
t3 = -qJ(5) * t94 + (-qJD(4) * t34 + t81) * t70 + t101;
t2 = t3 * t73;
t1 = t31 * pkin(4) + t81 * t73 + (-t32 + (-t24 + t109) * t70) * qJD(4) + t90;
t10 = [0, 0, 0, 0.2e1 * t72 * t103, 0.2e1 * (-t72 ^ 2 + t75 ^ 2) * qJD(2), 0, 0, 0, t72 * t102, t75 * t102, -0.2e1 * t120, 0.2e1 * t30 * t46 - 0.2e1 * t47 * t31, 0, 0, 0, 0.2e1 * t63 * t31 + 0.2e1 * t46 * t98, -0.2e1 * t63 * t30 + 0.2e1 * t47 * t98, -0.2e1 * t69 * t120 - 0.2e1 * t44 * t93, -t30 * t91 + 0.2e1 * t44 * t87, 0.2e1 * t47 * t115 - 0.2e1 * t78 * t46, -0.2e1 * t47 * t117 - 0.2e1 * t79 * t46, 0.2e1 * t46 * t31, 0.2e1 * t17 * t119 + 0.2e1 * t89 * t31 + 0.2e1 * t79 * t33 + 0.2e1 * t5 * t46, -0.2e1 * t112 * t31 + 0.2e1 * t17 * t118 - 0.2e1 * t78 * t33 + 0.2e1 * t4 * t46, -0.2e1 * (-t7 * t73 - t70 * t9) * t30 + 0.2e1 * (-t1 * t73 - t3 * t70 + (t7 * t70 - t73 * t9) * qJD(4)) * t47, 0.2e1 * t7 * t1 + 0.2e1 * t20 * t6 + 0.2e1 * t9 * t3; 0, 0, 0, 0, 0, t103, -t104, 0, -pkin(6) * t103, pkin(6) * t104, 0, 0, -t30, -t31, 0, -t17, t16, t12, t8, t19, -t77, 0, t26 + (-t82 * qJD(4) - t17) * t73 + t76 * t70, t82 * t105 + t76 * t73 + t113, t2 + (-t29 * t47 + t30 * t41 + (-t42 * t47 - t7) * qJD(4)) * t73 + (-t28 * t47 + t30 * t42 - t1 + (t41 * t47 - t9) * qJD(4)) * t70, t1 * t41 + t20 * t49 + t9 * t28 + t7 * t29 + t3 * t42 + t6 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t97, -0.2e1 * t96, t58, t45, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t111, -0.2e1 * t29 * t70 + 0.2e1 * t25 + 0.2e1 * (-t41 * t73 - t42 * t70) * qJD(4), 0.2e1 * t42 * t28 + 0.2e1 * t41 * t29 + 0.2e1 * t52 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31, 0, -t17, t16, t12, t8, t19, -t77, 0, t26 + (pkin(3) * t30 - pkin(8) * t31) * t70 + (-t17 + (-pkin(3) * t47 - pkin(8) * t46) * qJD(4)) * t73, t78 * pkin(3) + t77 * pkin(8) + t113, t2 + (t30 * t53 - t37 * t47 + (-t47 * t54 - t7) * qJD(4)) * t73 + (t30 * t54 - t36 * t47 - t1 + (t47 * t53 - t9) * qJD(4)) * t70, t1 * t53 + t20 * t64 + t3 * t54 + t9 * t36 + t7 * t37 + t6 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, t58, t45, 0, 0, 0, t80 - t100, -t99 + t111, t25 + t35 + (-t29 - t37) * t70 + ((-t41 - t53) * t73 + (-t42 - t54) * t70) * qJD(4), t28 * t54 + t29 * t53 + t42 * t36 + t41 * t37 + t49 * t62 + t52 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t45, 0, 0, 0, -0.2e1 * t100, -0.2e1 * t99, -0.2e1 * t37 * t70 + 0.2e1 * t35 + 0.2e1 * (-t53 * t73 - t54 * t70) * qJD(4), 0.2e1 * t54 * t36 + 0.2e1 * t53 * t37 + 0.2e1 * t62 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t79, t31, t5, t4, t78 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t105, 0, -t60 * t66 - t70 * t96, t60 * t105 - t85, -t95, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t105, 0, -pkin(8) * t66, pkin(8) * t105, -t95, t37 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
