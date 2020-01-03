% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:05
% EndTime: 2019-12-31 20:11:09
% DurationCPUTime: 1.20s
% Computational Cost: add. (1215->218), mult. (2807->316), div. (0->0), fcn. (1533->4), ass. (0->133)
t155 = pkin(3) + pkin(6);
t87 = cos(qJ(4));
t122 = t87 * qJD(2);
t88 = cos(qJ(2));
t132 = qJD(1) * t88;
t85 = sin(qJ(4));
t49 = -t85 * t132 + t122;
t120 = qJD(1) * qJD(2);
t111 = t88 * t120;
t86 = sin(qJ(2));
t112 = t86 * t120;
t125 = t85 * qJD(2);
t47 = t87 * t132 + t125;
t19 = t47 * qJD(4) - t85 * t112;
t110 = -t86 * qJ(3) - pkin(1);
t89 = -pkin(2) - pkin(7);
t44 = t89 * t88 + t110;
t26 = t44 * qJD(1);
t124 = t86 * qJD(1);
t74 = pkin(6) * t124;
t100 = pkin(3) * t124 + qJD(3) + t74;
t28 = t89 * qJD(2) + t100;
t10 = t87 * t26 + t85 * t28;
t72 = pkin(2) * t112;
t104 = pkin(7) * t86 - qJ(3) * t88;
t123 = t86 * qJD(3);
t94 = t104 * qJD(2) - t123;
t17 = t94 * qJD(1) + t72;
t71 = pkin(6) * t111;
t43 = pkin(3) * t111 + t71;
t109 = -t85 * t17 + t87 * t43;
t93 = -t10 * qJD(4) + t109;
t1 = pkin(4) * t111 + t19 * qJ(5) - t49 * qJD(5) + t93;
t20 = t49 * qJD(4) - t87 * t112;
t128 = qJD(4) * t87;
t117 = -t28 * t128 - t87 * t17 - t85 * t43;
t129 = qJD(4) * t85;
t97 = -t26 * t129 - t117;
t2 = -t20 * qJ(5) - t47 * qJD(5) + t97;
t9 = -t85 * t26 + t87 * t28;
t6 = -t49 * qJ(5) + t9;
t73 = qJD(4) + t124;
t5 = t73 * pkin(4) + t6;
t7 = -t47 * qJ(5) + t10;
t162 = -(t73 * t5 - t2) * t85 + (t73 * t7 + t1) * t87;
t161 = -0.2e1 * t120;
t64 = t155 * t86;
t137 = t87 * t44 + t85 * t64;
t149 = t47 * t73;
t160 = t19 - t149;
t148 = t49 * t73;
t159 = -t20 + t148;
t156 = t49 ^ 2;
t154 = t5 - t6;
t153 = pkin(4) * t88;
t78 = pkin(2) * t124;
t34 = t104 * qJD(1) + t78;
t75 = pkin(6) * t132;
t56 = pkin(3) * t132 + t75;
t108 = -t85 * t34 + t87 * t56;
t121 = t87 * qJD(5);
t133 = qJ(5) - t89;
t144 = t85 * t86;
t152 = t133 * t129 - t121 - (-qJ(5) * t144 + t153) * qJD(1) - t108;
t131 = qJD(2) * t86;
t55 = t155 * t131;
t81 = qJD(2) * qJD(3);
t31 = -qJD(1) * t55 + t81;
t151 = t31 * t85;
t150 = t31 * t87;
t147 = t49 * t88;
t146 = t73 * t87;
t145 = t73 * t89;
t143 = t87 * t19;
t91 = qJD(1) ^ 2;
t142 = t88 * t91;
t90 = qJD(2) ^ 2;
t141 = t90 * t86;
t140 = t90 * t88;
t139 = t87 * t34 + t85 * t56;
t134 = t87 * qJ(5);
t60 = t133 * t87;
t138 = -qJD(4) * t60 - t85 * qJD(5) - t124 * t134 - t139;
t65 = t155 * t88;
t83 = t86 ^ 2;
t84 = t88 ^ 2;
t136 = t83 - t84;
t135 = qJD(2) * pkin(2);
t62 = -t88 * pkin(2) + t110;
t38 = qJD(1) * t62;
t130 = qJD(2) * t88;
t127 = qJD(4) * t88;
t82 = qJD(2) * qJ(3);
t37 = t82 + t56;
t16 = t47 * pkin(4) + qJD(5) + t37;
t126 = t16 * qJD(2);
t119 = t86 * t146;
t118 = t86 * t142;
t116 = t85 * t127;
t115 = t73 * t128;
t114 = t87 * t127;
t107 = qJ(5) * t88 - t44;
t106 = pkin(1) * t161;
t105 = qJD(3) - t135;
t103 = -qJD(1) * t84 + t73 * t86;
t102 = -0.2e1 * qJD(2) * t38;
t101 = t73 * t85;
t95 = -t88 * t82 - t123;
t24 = t95 * qJD(1) + t72;
t77 = pkin(2) * t131;
t36 = t77 + t95;
t99 = pkin(6) * t90 + qJD(1) * t36 + t24;
t98 = t89 * t130 + t37 * t86;
t22 = t77 + t94;
t57 = t155 * t130;
t96 = t64 * t128 - t44 * t129 + t87 * t22 + t85 * t57;
t12 = t20 * pkin(4) + t31;
t58 = pkin(6) * t112 - t81;
t61 = t105 + t74;
t63 = -t75 - t82;
t92 = -t58 * t88 + (t61 * t88 + (t63 + t75) * t86) * qJD(2);
t67 = t87 * t111;
t59 = t133 * t85;
t53 = -qJ(3) * t132 + t78;
t52 = t87 * t64;
t46 = t47 ^ 2;
t42 = t87 * t57;
t27 = t38 * t124;
t14 = -t88 * t134 + t137;
t13 = t86 * pkin(4) + t107 * t85 + t52;
t4 = -t88 * t121 + (t86 * t122 + t116) * qJ(5) + t96;
t3 = pkin(4) * t130 + t42 + t107 * t128 + (-qJ(5) * t131 - qJD(4) * t64 + qJD(5) * t88 - t22) * t85;
t8 = [0, 0, 0, 0.2e1 * t86 * t111, t136 * t161, t140, -t141, 0, -pkin(6) * t140 + t86 * t106, pkin(6) * t141 + t88 * t106, t92, t86 * t102 + t99 * t88, t88 * t102 - t99 * t86, t92 * pkin(6) + t24 * t62 + t38 * t36, t19 * t85 * t88 + (t86 * t125 - t114) * t49, (-t47 * t85 + t49 * t87) * t131 + (t143 + t85 * t20 + (t47 * t87 + t49 * t85) * qJD(4)) * t88, -t73 * t114 - t19 * t86 + (t103 * t85 + t147) * qJD(2), t73 * t116 - t20 * t86 + (t103 * t87 - t47 * t88) * qJD(2), (t73 + t124) * t130, (-t85 * t22 + t42) * t73 - t55 * t47 + t65 * t20 + (-t37 * t122 + t109) * t86 + (-t10 * t86 - t137 * t73) * qJD(4) + (-t37 * t129 + t150 + ((-t85 * t44 + t52) * qJD(1) + t9) * qJD(2)) * t88, -t96 * t73 - t55 * t49 - t65 * t19 + ((qJD(2) * t37 + qJD(4) * t26) * t85 + t117) * t86 + (-t37 * t128 - t151 + (-t137 * qJD(1) - t10) * qJD(2)) * t88, t13 * t19 - t14 * t20 - t3 * t49 - t4 * t47 + (-t5 * t85 + t7 * t87) * t131 + (t1 * t85 - t2 * t87 + (t5 * t87 + t7 * t85) * qJD(4)) * t88, t2 * t14 + t7 * t4 + t1 * t13 + t5 * t3 + t12 * (t87 * t153 + t65) - t16 * pkin(4) * t116 + (-pkin(4) * t87 - t155) * t86 * t126; 0, 0, 0, -t118, t136 * t91, 0, 0, 0, t91 * pkin(1) * t86, pkin(1) * t142, ((-t63 - t82) * t86 + (t105 - t61) * t88) * qJD(1), -t53 * t132 + t27, 0.2e1 * t81 + (t38 * t88 + t53 * t86) * qJD(1), -t58 * qJ(3) - t63 * qJD(3) - t38 * t53 + (-t63 * t86 + (-t61 - t135) * t88) * qJD(1) * pkin(6), -t49 * t101 - t143, (-t20 - t148) * t87 + (t19 + t149) * t85, -t73 * t129 + t67 + (-t73 * t144 - t147) * qJD(1), -t115 + (-t119 + (t47 - t125) * t88) * qJD(1), -t73 * t132, qJ(3) * t20 + t151 - t108 * t73 + t100 * t47 + (-t85 * t145 + t37 * t87) * qJD(4) + (t98 * t87 - t9 * t88) * qJD(1), -qJ(3) * t19 + t150 + t139 * t73 + t100 * t49 + (-t87 * t145 - t37 * t85) * qJD(4) + (t10 * t88 - t98 * t85) * qJD(1), -t138 * t47 - t152 * t49 - t60 * t19 + t59 * t20 - t162, -t2 * t59 - t1 * t60 + t12 * (t85 * pkin(4) + qJ(3)) + t138 * t7 + t152 * t5 + (pkin(4) * t146 + t100) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t83 * t91 - t90, t63 * qJD(2) + t27 + t71, 0, 0, 0, 0, 0, -qJD(2) * t47 - t73 * t101 + t67, -t115 - qJD(2) * t49 + (-t88 * t125 - t119) * qJD(1), t159 * t85 + t160 * t87, -t126 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t47, -t46 + t156, -t160, t159, t111, t10 * t73 - t37 * t49 + t93, t37 * t47 + t9 * t73 - t97, pkin(4) * t19 - t154 * t47, t154 * t7 + (-t16 * t49 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 - t156, t7 * t47 + t5 * t49 + t12;];
tauc_reg = t8;
