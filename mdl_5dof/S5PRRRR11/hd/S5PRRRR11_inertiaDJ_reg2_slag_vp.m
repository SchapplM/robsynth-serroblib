% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:59
% DurationCPUTime: 1.38s
% Computational Cost: add. (2190->164), mult. (6703->305), div. (0->0), fcn. (6194->8), ass. (0->100)
t120 = cos(qJ(4));
t59 = sin(pkin(5));
t64 = cos(qJ(3));
t113 = t59 * t64;
t60 = cos(pkin(5));
t124 = pkin(2) * t60;
t63 = sin(qJ(3));
t56 = t63 * t124;
t46 = pkin(7) * t113 + t56;
t38 = pkin(8) * t113 + t46;
t36 = t120 * t38;
t127 = qJD(4) + qJD(5);
t58 = t59 ^ 2;
t126 = -0.2e1 * t60;
t125 = pkin(7) + pkin(8);
t123 = pkin(3) * t60;
t122 = pkin(4) * t60;
t62 = sin(qJ(4));
t121 = t62 * pkin(3);
t119 = cos(qJ(5));
t80 = t120 * t64 - t62 * t63;
t74 = t80 * qJD(3);
t27 = (-t80 * qJD(4) - t74) * t59;
t95 = t120 * t63;
t40 = (t62 * t64 + t95) * t59;
t61 = sin(qJ(5));
t76 = t59 * t80;
t73 = t61 * t76;
t109 = qJD(4) * t62;
t110 = qJD(3) * t64;
t96 = t59 * t110;
t90 = t109 * t113 + t62 * t96 + (qJD(3) + qJD(4)) * t59 * t95;
t92 = t119 * qJD(5);
t11 = qJD(5) * t73 + t119 * t90 - t61 * t27 + t40 * t92;
t71 = t119 * t76;
t25 = t40 * t61 - t71;
t118 = t25 * t11;
t108 = qJD(5) * t61;
t10 = -qJD(5) * t71 + t40 * t108 + t119 * t27 + t61 * t90;
t26 = t119 * t40 + t73;
t117 = t26 * t10;
t116 = t38 * t62;
t115 = t40 * t27;
t114 = t59 * t63;
t57 = t64 * t124;
t100 = t59 * t125;
t88 = t63 * t100;
t75 = t57 - t88 + t123;
t20 = t62 * t75 + t36;
t16 = pkin(9) * t76 + t20;
t112 = t61 * t16;
t111 = qJD(3) * t63;
t33 = t120 * t75;
t19 = t33 - t116;
t15 = -pkin(9) * t40 + t122 + t19;
t106 = pkin(7) * t114;
t53 = qJD(3) * t57;
t66 = (-t62 * (-pkin(8) * t114 - t106) - t36) * qJD(3) - t62 * t53;
t13 = -t20 * qJD(4) + t66;
t65 = t27 * pkin(9) + t13;
t103 = -qJD(4) * t33 - t62 * (-t64 * t100 - t56) * qJD(3) - t120 * (-qJD(3) * t88 + t53);
t12 = t38 * t109 + t103;
t69 = -t90 * pkin(9) - t12;
t107 = -t119 * t69 - t15 * t92 - t61 * t65;
t104 = t61 * t121;
t102 = pkin(3) * t109;
t101 = pkin(4) * t108;
t99 = t120 * pkin(3);
t98 = t58 * t110;
t97 = t59 * t111;
t94 = t119 * t16;
t93 = t119 * t62;
t91 = -0.2e1 * t101;
t89 = pkin(3) * t97;
t87 = t63 * t98;
t86 = t60 * t97;
t84 = qJD(4) * t99;
t83 = pkin(4) * t92;
t82 = t99 + pkin(4);
t81 = pkin(3) * t90;
t1 = t16 * t108 + t107;
t9 = t61 * t15 + t94;
t78 = t119 * t82;
t28 = -qJD(5) * t78 + t127 * t104 - t119 * t84;
t72 = pkin(3) * qJD(4) * t76;
t70 = t90 * t76;
t68 = t119 * t65 - t61 * t69;
t67 = t127 * (-t61 * t120 - t93) * pkin(3);
t2 = -t9 * qJD(5) + t68;
t47 = (-t64 * pkin(3) - pkin(2)) * t59;
t45 = t57 - t106;
t44 = pkin(3) * t93 + t61 * t82;
t43 = t78 - t104;
t42 = t46 * qJD(3);
t41 = pkin(7) * t97 - t53;
t30 = -t59 * pkin(2) - pkin(3) * t113 - pkin(4) * t76;
t29 = t67 - t101;
t21 = t90 * pkin(4) + t89;
t8 = t119 * t15 - t112;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t70 - 0.2e1 * t115, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t117 + 0.2e1 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t41 * t63 - t42 * t64 + (-t45 * t63 + t46 * t64) * qJD(3)) * t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, pkin(3) * t86 - t40 * t12 + t13 * t76 - t90 * t19 - t27 * t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t26 - t10 * t9 - t11 * t8 - t2 * t25 + t21 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t87, 0.2e1 * (-t63 ^ 2 + t64 ^ 2) * t58 * qJD(3), 0.2e1 * t60 * t96, -0.2e1 * t87, -0.2e1 * t86, 0, -0.2e1 * pkin(2) * t58 * t111 - 0.2e1 * t42 * t60, -0.2e1 * pkin(2) * t98 + 0.2e1 * t41 * t60, 0.2e1 * (-t41 * t64 + t42 * t63 + (-t45 * t64 - t46 * t63) * qJD(3)) * t59, -0.2e1 * t41 * t46 - 0.2e1 * t42 * t45, -0.2e1 * t115, -0.2e1 * t27 * t76 - 0.2e1 * t40 * t90, t27 * t126, -0.2e1 * t70, t90 * t126, 0, -0.2e1 * t58 * t63 * pkin(3) * t74 + 0.2e1 * t13 * t60 + 0.2e1 * t47 * t90, 0.2e1 * t12 * t60 - 0.2e1 * t27 * t47 + 0.2e1 * t40 * t89, -0.2e1 * t12 * t76 - 0.2e1 * t13 * t40 + 0.2e1 * t19 * t27 - 0.2e1 * t20 * t90, -0.2e1 * t12 * t20 + 0.2e1 * t13 * t19 + 0.2e1 * t47 * t89, -0.2e1 * t117, 0.2e1 * t25 * t10 - 0.2e1 * t26 * t11, t10 * t126, 0.2e1 * t118, t11 * t126, 0, 0.2e1 * t11 * t30 + 0.2e1 * t2 * t60 + 0.2e1 * t21 * t25, 0.2e1 * t1 * t60 - 0.2e1 * t10 * t30 + 0.2e1 * t21 * t26, 0.2e1 * t1 * t25 + 0.2e1 * t10 * t8 - 0.2e1 * t11 * t9 - 0.2e1 * t2 * t26, -0.2e1 * t1 * t9 + 0.2e1 * t2 * t8 + 0.2e1 * t21 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t27, 0, -t120 * t81 - t27 * t121 + t40 * t84 - t62 * t72, 0, 0, 0, 0, 0, 0, -t11, t10, 0, -t10 * t44 - t11 * t43 - t25 * t29 - t26 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, -t97, 0, -t42, t41, 0, 0, 0, 0, -t27, 0, -t90, 0, (-t36 + (t125 * t114 - 0.2e1 * t123 - t57) * t62) * qJD(4) + t66, (-t60 * t99 + t116) * qJD(4) + t103, t40 * t102 + t120 * t72 + t27 * t99 - t62 * t81, (t120 * t13 - t12 * t62 + (t120 * t20 - t19 * t62) * qJD(4)) * pkin(3), 0, 0, -t10, 0, -t11, 0, t29 * t60 + t2, t28 * t60 + t1, t10 * t43 - t11 * t44 + t25 * t28 - t26 * t29, -t1 * t44 + t2 * t43 - t28 * t9 + t29 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t102, -0.2e1 * t84, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t29, 0.2e1 * t28, 0, -0.2e1 * t28 * t44 + 0.2e1 * t29 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, (-t119 * t11 - t10 * t61 + (t119 * t26 + t25 * t61) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, -t90, 0, t13, t12, 0, 0, 0, 0, -t10, 0, -t11, 0, (-t94 + (-t15 - t122) * t61) * qJD(5) + t68, (-t119 * t122 + t112) * qJD(5) + t107, (t119 * t10 - t11 * t61 + (-t119 * t25 + t26 * t61) * qJD(5)) * pkin(4), (t119 * t2 - t1 * t61 + (t119 * t9 - t61 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t84, 0, 0, 0, 0, 0, 0, 0, 0, t67 + t91, -t83 + t28, 0, (t119 * t29 - t28 * t61 + (t119 * t44 - t43 * t61) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -0.2e1 * t83, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, -t11, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t83, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
