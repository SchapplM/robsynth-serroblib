% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:49
% EndTime: 2021-01-16 02:54:56
% DurationCPUTime: 1.49s
% Computational Cost: add. (1890->190), mult. (4727->356), div. (0->0), fcn. (4602->10), ass. (0->113)
t123 = cos(pkin(11));
t69 = cos(qJ(3));
t101 = t123 * t69;
t67 = sin(qJ(3));
t122 = qJD(3) * t67;
t64 = sin(pkin(11));
t43 = qJD(3) * t101 - t64 * t122;
t102 = t123 * t67;
t48 = t64 * t69 + t102;
t68 = cos(qJ(5));
t61 = qJD(5) * t68;
t66 = sin(qJ(5));
t84 = t43 * t66 + t48 * t61;
t131 = t64 * t67;
t47 = -t101 + t131;
t58 = -pkin(3) * t69 - pkin(2);
t33 = pkin(4) * t47 - pkin(9) * t48 + t58;
t128 = qJ(4) + pkin(8);
t53 = t128 * t69;
t37 = t123 * t53 - t128 * t131;
t147 = t66 * t33 + t68 * t37;
t132 = t48 * t68;
t65 = sin(pkin(6));
t70 = cos(qJ(2));
t130 = t65 * t70;
t108 = qJD(2) * t130;
t138 = sin(qJ(2));
t106 = t65 * t138;
t124 = cos(pkin(6));
t78 = t67 * t106 - t124 * t69;
t145 = t78 * qJD(3) - t69 * t108;
t44 = t69 * t106 + t124 * t67;
t73 = t44 * qJD(3) + t67 * t108;
t14 = t123 * t73 - t145 * t64;
t114 = t66 * t130;
t28 = t123 * t44 - t64 * t78;
t19 = t68 * t28 - t114;
t27 = t123 * t78 + t44 * t64;
t42 = t48 * qJD(3);
t18 = t68 * t130 + t66 * t28;
t71 = -t123 * t145 - t64 * t73;
t103 = qJD(2) * t138;
t97 = t65 * t103;
t6 = t18 * qJD(5) - t66 * t97 - t68 * t71;
t119 = qJD(5) * t66;
t83 = t48 * t119 - t43 * t68;
t146 = -t14 * t132 + t19 * t42 + t83 * t27 - t6 * t47;
t7 = -qJD(5) * t114 + t28 * t61 + t66 * t71 - t68 * t97;
t144 = -t6 * t66 - t7 * t68 + (t18 * t66 + t19 * t68) * qJD(5);
t118 = qJD(6) * t47;
t125 = qJ(6) * t42;
t100 = qJD(3) * t128;
t40 = t69 * qJD(4) - t67 * t100;
t81 = -t67 * qJD(4) - t69 * t100;
t24 = t123 * t40 + t64 * t81;
t60 = pkin(3) * t122;
t25 = pkin(4) * t42 - pkin(9) * t43 + t60;
t3 = t37 * t119 - t68 * t24 - t66 * t25 - t33 * t61;
t1 = t118 - t3 + t125;
t10 = qJ(6) * t47 + t147;
t88 = t33 * t68 - t37 * t66;
t11 = -pkin(5) * t47 - t88;
t139 = t42 * pkin(5);
t4 = -qJD(5) * t147 - t66 * t24 + t68 * t25;
t2 = -t139 - t4;
t143 = t1 * t66 - t2 * t68 + (t10 * t68 + t11 * t66) * qJD(5);
t94 = t68 * pkin(5) + t66 * qJ(6);
t142 = t94 * qJD(5) - t68 * qJD(6);
t141 = 0.2e1 * qJD(5);
t140 = 0.2e1 * qJD(6);
t137 = t27 * t14;
t56 = pkin(3) * t64 + pkin(9);
t136 = t42 * t56;
t134 = t47 * t56;
t133 = t48 * t66;
t129 = t66 * t68;
t62 = t66 ^ 2;
t63 = t68 ^ 2;
t126 = t62 - t63;
t121 = qJD(3) * t69;
t120 = qJD(3) * t70;
t117 = t66 * qJD(6);
t115 = -0.2e1 * pkin(2) * qJD(3);
t57 = -t123 * pkin(3) - pkin(4);
t113 = t57 * t141;
t112 = t67 * t120;
t110 = t56 * t119;
t109 = t56 * t61;
t107 = t66 * t61;
t105 = -0.4e1 * t48 * t129;
t23 = -t123 * t81 + t40 * t64;
t36 = t128 * t102 + t64 * t53;
t99 = t126 * qJD(5);
t93 = pkin(5) * t66 - qJ(6) * t68;
t91 = -t10 * t66 + t11 * t68;
t90 = t18 * t68 - t19 * t66;
t86 = t43 * t57 - t136;
t85 = -t48 * t57 + t134;
t32 = t42 * t66 + t47 * t61;
t30 = -t47 * t119 + t42 * t68;
t82 = t14 * t133 - t18 * t42 + t84 * t27 - t47 * t7;
t45 = -t94 + t57;
t5 = t142 * t48 + t93 * t43 + t23;
t80 = -t5 + (t45 * t48 - t134) * qJD(5);
t15 = t93 * t48 + t36;
t41 = -pkin(5) * t119 + qJ(6) * t61 + t117;
t76 = qJD(5) * t15 - t41 * t48 + t43 * t45 - t136;
t75 = t91 * qJD(5) + t1 * t68 + t2 * t66;
t74 = t90 * qJD(5) - t6 * t68 + t7 * t66;
t46 = t48 ^ 2;
t9 = t27 * t119 - t14 * t68;
t8 = t14 * t66 + t27 * t61;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t65 ^ 2 * t70 * t103 + 0.2e1 * t28 * t71 + 0.2e1 * t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t18 * t7 - 0.2e1 * t19 * t6 + 0.2e1 * t137; 0, 0, -t97, -t108, 0, 0, 0, 0, 0, (-t69 * t103 - t112) * t65, (t67 * t103 - t69 * t120) * t65, (t47 * t103 - t42 * t70) * t65, (t48 * t103 - t43 * t70) * t65, t14 * t48 + t27 * t43 - t28 * t42 - t71 * t47, -t65 * pkin(3) * t112 + t14 * t36 + t27 * t23 + t28 * t24 + t71 * t37 + t58 * t97, 0, 0, 0, 0, 0, t82, -t146, t82, -t144 * t48 + t90 * t43, t146, t1 * t19 - t10 * t6 + t11 * t7 + t14 * t15 + t18 * t2 + t27 * t5; 0, 0, 0, 0, 0.2e1 * t67 * t121, 0.2e1 * (-t67 ^ 2 + t69 ^ 2) * qJD(3), 0, 0, 0, t67 * t115, t69 * t115, 0.2e1 * t42 * t58 + 0.2e1 * t47 * t60, 0.2e1 * t43 * t58 + 0.2e1 * t48 * t60, 0.2e1 * t23 * t48 - 0.2e1 * t24 * t47 + 0.2e1 * t36 * t43 - 0.2e1 * t37 * t42, 0.2e1 * t23 * t36 + 0.2e1 * t24 * t37 + 0.2e1 * t58 * t60, 0.2e1 * t43 * t48 * t63 - 0.2e1 * t46 * t107, t126 * t46 * t141 + t43 * t105, 0.2e1 * t42 * t132 - 0.2e1 * t83 * t47, -0.2e1 * t42 * t133 - 0.2e1 * t84 * t47, 0.2e1 * t47 * t42, 0.2e1 * t23 * t133 + 0.2e1 * t84 * t36 + 0.2e1 * t4 * t47 + 0.2e1 * t88 * t42, 0.2e1 * t23 * t132 - 0.2e1 * t147 * t42 + 0.2e1 * t3 * t47 - 0.2e1 * t83 * t36, -0.2e1 * t11 * t42 + 0.2e1 * t5 * t133 + 0.2e1 * t15 * t84 - 0.2e1 * t2 * t47, -0.2e1 * t143 * t48 + 0.2e1 * t91 * t43, 0.2e1 * t1 * t47 + 0.2e1 * t10 * t42 - 0.2e1 * t5 * t132 + 0.2e1 * t15 * t83, 0.2e1 * t1 * t10 + 0.2e1 * t11 * t2 + 0.2e1 * t15 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t145, -t14, -t71, 0, (-t14 * t123 + t71 * t64) * pkin(3), 0, 0, 0, 0, 0, t9, t8, t9, t74, -t8, t14 * t45 - t27 * t41 + t56 * t74; 0, 0, 0, 0, 0, 0, t121, -t122, 0, -pkin(8) * t121, pkin(8) * t122, -t23, -t24, (-t123 * t43 - t42 * t64) * pkin(3), (-t123 * t23 + t24 * t64) * pkin(3), t43 * t129 - t48 * t99, qJD(5) * t105 - t126 * t43, t32, t30, 0, -t23 * t68 + t86 * t66 + (t36 * t66 - t85 * t68) * qJD(5), t23 * t66 + t86 * t68 + (t36 * t68 + t66 * t85) * qJD(5), t66 * t76 + t68 * t80, t75, t66 * t80 - t68 * t76, -t15 * t41 + t45 * t5 + t56 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t107, -0.2e1 * t99, 0, 0, 0, t66 * t113, t68 * t113, 0.2e1 * t45 * t119 + 0.2e1 * t41 * t68, 0, 0.2e1 * t41 * t66 - 0.2e1 * t45 * t61, -0.2e1 * t45 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t43, 0, t60, 0, 0, 0, 0, 0, t30, -t32, t30, (-t62 - t63) * t43, t32, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t7, 0, -t6, -pkin(5) * t7 - qJ(6) * t6 + qJD(6) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t84, t42, t4, t3, t4 + 0.2e1 * t139, -t94 * t43 + (qJD(5) * t93 - t117) * t48, 0.2e1 * t118 - t3 + 0.2e1 * t125, -pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t119, 0, -t109, t110, -t109, -t142, -t110, -t142 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, -t61, -t119, 0, t61, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, qJ(6) * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t83, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
