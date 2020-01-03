% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:48
% DurationCPUTime: 1.05s
% Computational Cost: add. (1787->212), mult. (4735->309), div. (0->0), fcn. (3234->6), ass. (0->121)
t106 = (qJD(1) * qJD(2));
t142 = -2 * t106;
t75 = sin(pkin(7));
t79 = cos(qJ(2));
t113 = cos(pkin(7));
t77 = sin(qJ(2));
t99 = t113 * t77;
t58 = t75 * t79 + t99;
t112 = qJD(1) * t58;
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t36 = t76 * qJD(2) + t112 * t78;
t110 = qJD(1) * t77;
t98 = t113 * t79;
t65 = qJD(1) * t98;
t48 = t75 * t110 - t65;
t45 = qJD(4) + t48;
t96 = t45 * t76;
t141 = t36 * t96;
t117 = -qJ(3) - pkin(5);
t97 = qJD(2) * t117;
t46 = t79 * qJD(3) + t77 * t97;
t40 = t46 * qJD(1);
t100 = t113 * t40;
t84 = t77 * qJD(3) - t79 * t97;
t83 = t75 * t84;
t14 = -qJD(1) * t83 + t100;
t50 = t58 * qJD(2);
t42 = qJD(1) * t50;
t101 = t77 * t106;
t64 = t75 * t101;
t43 = qJD(2) * t65 - t64;
t66 = pkin(2) * t101;
t17 = t42 * pkin(3) - t43 * pkin(6) + t66;
t71 = -t79 * pkin(2) - pkin(1);
t111 = qJD(1) * t71;
t62 = qJD(3) + t111;
t18 = t48 * pkin(3) - pkin(6) * t112 + t62;
t63 = t117 * t79;
t61 = qJD(1) * t63;
t54 = t113 * t61;
t114 = qJD(2) * pkin(2);
t102 = t117 * t77;
t60 = qJD(1) * t102;
t56 = t60 + t114;
t28 = t75 * t56 - t54;
t24 = qJD(2) * pkin(6) + t28;
t5 = t78 * t18 - t76 * t24;
t1 = qJD(4) * t5 + t78 * t14 + t76 * t17;
t140 = -t5 * t45 + t1;
t6 = t76 * t18 + t78 * t24;
t2 = -qJD(4) * t6 - t76 * t14 + t78 * t17;
t139 = t6 * t45 + t2;
t138 = t112 ^ 2;
t137 = pkin(2) * t77;
t82 = t113 * t84;
t13 = qJD(1) * t82 + t75 * t40;
t32 = -t117 * t99 - t75 * t63;
t134 = t13 * t32;
t133 = t13 * t76;
t107 = t78 * qJD(2);
t109 = qJD(4) * t76;
t15 = -qJD(4) * t107 + t109 * t112 - t78 * t43;
t132 = t15 * t76;
t121 = t76 * t43;
t16 = t36 * qJD(4) + t121;
t131 = t16 * t78;
t34 = t112 * t76 - t107;
t130 = t34 * t48;
t129 = t36 * t34;
t128 = t36 * t112;
t86 = -t75 * t77 + t98;
t127 = t42 * t86;
t126 = t112 * t34;
t125 = t112 * t48;
t124 = t58 * t78;
t123 = t75 * t61;
t11 = t76 * t16;
t122 = t76 * t42;
t38 = t78 * t42;
t81 = qJD(1) ^ 2;
t120 = t79 * t81;
t80 = qJD(2) ^ 2;
t119 = t80 * t77;
t118 = t80 * t79;
t108 = qJD(4) * t78;
t116 = -t34 * t108 - t11;
t115 = t77 ^ 2 - t79 ^ 2;
t105 = t77 * t120;
t104 = t77 * t114;
t103 = pkin(2) * t110;
t95 = t45 * t78;
t94 = 0.2e1 * t112;
t93 = pkin(1) * t142;
t92 = t79 * t101;
t91 = -t5 * t78 - t6 * t76;
t26 = -pkin(3) * t86 - t58 * pkin(6) + t71;
t33 = t75 * t102 - t113 * t63;
t9 = t78 * t26 - t76 * t33;
t10 = t76 * t26 + t78 * t33;
t89 = -t45 * t109 - t48 * t96 + t38;
t53 = t86 * qJD(2);
t88 = t58 * t108 + t53 * t76;
t87 = -t58 * t109 + t53 * t78;
t27 = t113 * t56 + t123;
t23 = -qJD(2) * pkin(3) - t27;
t68 = t75 * pkin(2) + pkin(6);
t85 = t45 * t23 - t42 * t68;
t69 = -t113 * pkin(2) - pkin(3);
t47 = t48 ^ 2;
t30 = t113 * t60 + t123;
t29 = t75 * t60 - t54;
t22 = t50 * pkin(3) - t53 * pkin(6) + t104;
t21 = pkin(3) * t112 + t48 * pkin(6) + t103;
t20 = t113 * t46 - t83;
t19 = t75 * t46 + t82;
t8 = t76 * t21 + t78 * t30;
t7 = t78 * t21 - t76 * t30;
t4 = -t10 * qJD(4) - t76 * t20 + t78 * t22;
t3 = t9 * qJD(4) + t78 * t20 + t76 * t22;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t92, t115 * t142, t118, -0.2e1 * t92, -t119, 0, -pkin(5) * t118 + t77 * t93, pkin(5) * t119 + t79 * t93, 0, 0, t112 * t53 + t43 * t58, -t112 * t50 - t58 * t42 + t43 * t86 - t53 * t48, t53 * qJD(2), t48 * t50 - t127, -t50 * qJD(2), 0, t71 * t42 + t62 * t50 + (-t19 + (-qJD(1) * t86 + t48) * t137) * qJD(2), t71 * t43 + t62 * t53 + (t94 * t137 - t20) * qJD(2), t112 * t19 + t13 * t58 + t14 * t86 - t20 * t48 - t27 * t53 - t28 * t50 + t32 * t43 - t33 * t42, t134 + t14 * t33 - t27 * t19 + t28 * t20 + (t62 + t111) * t104, -t15 * t124 + t87 * t36, (-t34 * t78 - t36 * t76) * t53 + (t132 - t131 + (t34 * t76 - t36 * t78) * qJD(4)) * t58, t15 * t86 + t36 * t50 + t58 * t38 + t87 * t45, t58 * t11 + t88 * t34, -t58 * t122 + t16 * t86 - t34 * t50 - t88 * t45, t45 * t50 - t127, t58 * t133 + t32 * t16 + t19 * t34 - t2 * t86 + t88 * t23 + t4 * t45 + t9 * t42 + t5 * t50, t1 * t86 - t10 * t42 + t13 * t124 - t32 * t15 + t19 * t36 + t23 * t87 - t3 * t45 - t6 * t50, -t10 * t16 + t9 * t15 - t3 * t34 - t4 * t36 + t91 * t53 + (-t1 * t76 - t2 * t78 + (t5 * t76 - t6 * t78) * qJD(4)) * t58, t1 * t10 + t23 * t19 + t2 * t9 + t6 * t3 + t5 * t4 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t115 * t81, 0, t105, 0, 0, t81 * pkin(1) * t77, pkin(1) * t120, 0, 0, t125, -t47 + t138, -t64 + (t65 + t48) * qJD(2), -t125, 0, 0, t29 * qJD(2) - t48 * t103 - t112 * t62 - t13, -t100 + t30 * qJD(2) + t62 * t48 + (-t112 * t137 + t83) * qJD(1), (t28 - t29) * t112 + (-t27 + t30) * t48 + (-t113 * t43 - t42 * t75) * pkin(2), t27 * t29 - t28 * t30 + (-t62 * t110 - t113 * t13 + t14 * t75) * pkin(2), t36 * t95 - t132, (-t15 - t130) * t78 - t141 + t116, t45 * t95 + t122 - t128, t34 * t96 - t131, t89 + t126, -t45 * t112, -t13 * t78 + t69 * t16 - t29 * t34 - t5 * t112 + (-t68 * t108 - t7) * t45 + t85 * t76, t133 - t69 * t15 - t29 * t36 + t6 * t112 + (t109 * t68 + t8) * t45 + t85 * t78, t8 * t34 + t7 * t36 + (-t16 * t68 - t48 * t5 + t1 + (t36 * t68 - t5) * qJD(4)) * t78 + (-t15 * t68 - t48 * t6 - t2 + (t34 * t68 - t6) * qJD(4)) * t76, t13 * t69 - t23 * t29 - t5 * t7 - t6 * t8 + (qJD(4) * t91 + t1 * t78 - t2 * t76) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * qJD(2), -t64 + (t65 - t48) * qJD(2), -t47 - t138, t112 * t27 + t28 * t48 + t66, 0, 0, 0, 0, 0, 0, t89 - t126, -t45 ^ 2 * t78 - t122 - t128, (t15 - t130) * t78 + t141 + t116, -t112 * t23 + t139 * t78 + t140 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t34 ^ 2 + t36 ^ 2, t34 * t45 - t15, -t129, -t121 + (-qJD(4) + t45) * t36, t42, -t23 * t36 + t139, t23 * t34 - t140, 0, 0;];
tauc_reg = t12;
