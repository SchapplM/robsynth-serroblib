% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:39
% DurationCPUTime: 1.41s
% Computational Cost: add. (1108->227), mult. (2771->319), div. (0->0), fcn. (1766->6), ass. (0->132)
t83 = cos(qJ(3));
t132 = qJD(3) * t83;
t82 = sin(qJ(2));
t113 = t82 * t132;
t147 = t82 * t83;
t80 = cos(pkin(8));
t121 = t80 * t147;
t81 = sin(qJ(3));
t133 = qJD(3) * t81;
t114 = t82 * t133;
t84 = cos(qJ(2));
t146 = t84 * t80;
t79 = sin(pkin(8));
t153 = t79 * t82;
t92 = t83 * t146 + t153;
t26 = t92 * qJD(2) - t80 * t114;
t49 = -t84 * t79 + t121;
t129 = t79 * qJD(3);
t135 = qJD(2) * t81;
t58 = t80 * t135 + t129;
t167 = (((-t49 + t121) * qJD(3) + t58 * t84) * t81 + t26 * t83) * qJD(2) + t58 * t113;
t91 = t84 * t135 + t113;
t122 = qJD(2) * qJD(3);
t165 = -0.2e1 * t122;
t95 = pkin(3) * t81 - qJ(4) * t83;
t46 = t95 * qJD(3) - t81 * qJD(4);
t164 = -t92 * qJD(1) + t79 * t46;
t54 = t58 ^ 2;
t128 = t80 * qJD(3);
t56 = t79 * t135 - t128;
t163 = -t56 ^ 2 - t54;
t125 = t84 * qJD(1);
t127 = t82 * qJD(1);
t152 = t79 * t83;
t162 = t125 * t152 + (-t127 + t46) * t80;
t160 = pkin(6) * t80;
t109 = t83 * t122;
t100 = t80 * t109;
t110 = qJD(2) * t125;
t70 = t81 * t110;
t71 = qJD(2) * pkin(6) + t127;
t33 = t71 * t132 + t70;
t69 = t79 * t109;
t89 = pkin(4) * t69 - qJ(5) * t100 + t33;
t5 = -t58 * qJD(5) + t89;
t159 = t5 * t79;
t158 = t5 * t80;
t156 = t33 * t79;
t155 = t33 * t80;
t61 = t95 * qJD(2);
t150 = t80 * t61;
t136 = t81 * qJ(4);
t96 = -t83 * pkin(3) - t136;
t65 = -pkin(2) + t96;
t149 = t80 * t65;
t62 = t81 * t71;
t148 = t81 * t82;
t63 = t83 * t71;
t85 = qJD(3) ^ 2;
t145 = t85 * t81;
t144 = t85 * t83;
t126 = t83 * qJD(5);
t143 = -t126 + (qJ(5) - t160) * t133 + t164;
t112 = pkin(6) * t79 + pkin(4);
t142 = -t112 * t133 - t162;
t120 = pkin(6) * t133;
t105 = t79 * t120;
t141 = t105 + t162;
t104 = t80 * t120;
t140 = -t104 + t164;
t28 = (t46 + t127) * qJD(2);
t31 = t83 * t110 + (qJD(4) - t62) * qJD(3);
t4 = t79 * t28 + t80 * t31;
t39 = t65 * qJD(2) - t125;
t123 = qJD(3) * qJ(4);
t51 = t63 + t123;
t10 = t79 * t39 + t80 * t51;
t35 = t83 * t160 + t79 * t65;
t78 = t83 ^ 2;
t139 = t81 ^ 2 - t78;
t86 = qJD(2) ^ 2;
t138 = t85 + t86;
t137 = qJD(2) * pkin(2);
t134 = qJD(2) * t83;
t131 = qJD(4) * t58;
t130 = qJD(4) * t80;
t124 = qJ(5) * qJD(2);
t119 = t56 * t125;
t118 = t58 * t125;
t117 = t81 * t125;
t116 = t79 * t134;
t111 = t81 * t123;
t108 = t81 * t122;
t3 = t80 * t28 - t79 * t31;
t72 = -t125 - t137;
t107 = -t72 - t125;
t94 = pkin(4) * t79 - qJ(5) * t80;
t20 = t94 * t134 + t63;
t106 = qJD(5) * t79 + t20;
t103 = -qJD(3) * pkin(3) + qJD(4);
t101 = t84 * t165;
t99 = t83 * t108;
t97 = -t58 * t134 + t69;
t47 = t103 + t62;
t9 = t80 * t39 - t79 * t51;
t93 = pkin(6) + t94;
t90 = qJD(3) * (-t107 - t137);
t25 = -qJD(2) * t82 * t80 - t79 * t114 + t84 * t116;
t48 = t79 * t147 + t146;
t88 = t99 * t153 + (-t48 * t133 + t25 * t83) * qJD(2) + t91 * t56;
t87 = t25 * t58 - t26 * t56 + (t48 * t80 - t49 * t79) * t109;
t68 = qJD(4) * t116;
t64 = -t80 * pkin(4) - t79 * qJ(5) - pkin(3);
t50 = t79 * t61;
t45 = t56 * t134;
t44 = t56 * t130;
t40 = t93 * t81;
t34 = -pkin(6) * t152 + t149;
t30 = t112 * t83 - t149;
t29 = -t83 * qJ(5) + t35;
t27 = t45 + t100;
t18 = -t80 * t62 + t50;
t17 = t79 * t62 + t150;
t16 = -t80 * t81 * qJD(5) + t93 * t132;
t14 = -t150 + (-pkin(4) * qJD(2) - t71 * t79) * t81;
t13 = t50 + (-t71 * t80 + t124) * t81;
t8 = t56 * pkin(4) - t58 * qJ(5) + t47;
t7 = -t83 * t124 + t10;
t6 = pkin(4) * t134 + qJD(5) - t9;
t2 = -pkin(4) * t108 - t3;
t1 = (qJ(5) * t133 - t126) * qJD(2) + t4;
t11 = [0, 0, -t86 * t82, -t86 * t84, 0, 0, 0, 0, 0, t81 * t101 - t138 * t147, t83 * t101 + t138 * t148, t88, t167, t87, t10 * t26 + t33 * t148 - t9 * t25 - t3 * t48 + t4 * t49 + t91 * t47, t88, t87, -t167, t1 * t49 + t5 * t148 + t2 * t48 + t6 * t25 + t7 * t26 + t8 * t91; 0, 0, 0, 0, 0.2e1 * t99, t139 * t165, t144, -t145, 0, -pkin(6) * t144 + t81 * t90, pkin(6) * t145 + t83 * t90, (-t119 + t156 + (qJD(2) * t34 + t9) * qJD(3)) * t81 + (-t3 + (pkin(6) * t56 + t47 * t79) * qJD(3) + (t105 - t141) * qJD(2)) * t83, (-t118 + t155 + (-qJD(2) * t35 - t10) * qJD(3)) * t81 + (t4 + (pkin(6) * t58 + t47 * t80) * qJD(3) + (t104 + t140) * qJD(2)) * t83, (-t3 * t80 - t4 * t79) * t81 - t141 * t58 - t140 * t56 + (-t10 * t79 - t80 * t9 + (-t34 * t80 - t35 * t79) * qJD(2)) * t132, -t47 * t117 + t3 * t34 + t4 * t35 + t141 * t9 + t140 * t10 + (t47 * t132 + t33 * t81) * pkin(6), t16 * t56 + (-t119 + t159 + (-qJD(2) * t30 - t6) * qJD(3)) * t81 + (t8 * t129 + t2 + (t40 * t129 + t142) * qJD(2)) * t83, (-t1 * t79 + t2 * t80) * t81 + t142 * t58 - t143 * t56 + (t6 * t80 - t7 * t79 + (-t29 * t79 + t30 * t80) * qJD(2)) * t132, -t16 * t58 + (t118 - t158 + (qJD(2) * t29 + t7) * qJD(3)) * t81 + (-t8 * t128 - t1 + (-t128 * t40 - t143) * qJD(2)) * t83, t1 * t29 + t2 * t30 + t5 * t40 + (t16 - t117) * t8 + t143 * t7 + t142 * t6; 0, 0, 0, 0, -t81 * t86 * t83, t139 * t86, 0, 0, 0, -t72 * t135 - t70, t107 * t134, -t56 * t63 - t155 + t68 + (t17 * t83 - t81 * t9 + (t96 * qJD(3) - t47 * t83) * t79) * qJD(2), -t58 * t63 + t156 + (t10 * t81 - t18 * t83 + (-t111 + (t103 - t47) * t83) * t80) * qJD(2), t17 * t58 + t18 * t56 - t44 + (t9 * t134 + t4) * t80 + (t10 * t134 + t131 - t3) * t79, -t47 * t63 - t33 * pkin(3) - t10 * t18 - t9 * t17 + (t10 * t80 - t79 * t9) * qJD(4) + (-t3 * t79 + t4 * t80) * qJ(4), -t158 + t68 - t106 * t56 + (-t14 * t83 + t6 * t81 + (-t8 * t83 + (t64 * t83 - t136) * qJD(3)) * t79) * qJD(2), t13 * t56 - t14 * t58 - t44 + (-t6 * t134 + t1) * t80 + (t7 * t134 + t131 + t2) * t79, -t159 + t106 * t58 + (t13 * t83 - t7 * t81 + (t111 + (-qJD(3) * t64 - qJD(4) + t8) * t83) * t80) * qJD(2), t1 * t80 * qJ(4) - t6 * t14 - t8 * t20 + t5 * t64 + (-t13 + t130) * t7 + (qJ(4) * t2 + qJD(4) * t6 - qJD(5) * t8) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t27, t163, t10 * t56 + t9 * t58 + t33, t97, t163, -t27, t7 * t56 + (-qJD(5) - t6) * t58 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56 - t108, -t45 + t100, -t78 * t86 - t54, t8 * t58 + (-pkin(4) * t133 + t7 * t83) * qJD(2) - t3;];
tauc_reg = t11;
