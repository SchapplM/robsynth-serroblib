% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:41
% DurationCPUTime: 1.19s
% Computational Cost: add. (1468->238), mult. (3395->320), div. (0->0), fcn. (1986->6), ass. (0->130)
t79 = cos(qJ(3));
t125 = t79 * qJD(1);
t64 = -qJD(4) + t125;
t78 = cos(qJ(4));
t129 = qJD(4) * t78;
t76 = sin(qJ(4));
t130 = qJD(4) * t76;
t66 = sin(pkin(8)) * pkin(1) + pkin(6);
t60 = t66 * qJD(1);
t77 = sin(qJ(3));
t70 = t77 * qJD(2);
t37 = t79 * t60 + t70;
t30 = qJD(3) * pkin(7) + t37;
t165 = t79 * qJD(2) - t77 * t60;
t31 = t165 * qJD(3);
t67 = -cos(pkin(8)) * pkin(1) - pkin(2);
t47 = -t79 * pkin(3) - t77 * pkin(7) + t67;
t33 = t47 * qJD(1);
t97 = pkin(3) * t77 - pkin(7) * t79;
t58 = t97 * qJD(3);
t46 = qJD(1) * t58;
t103 = t30 * t129 + t33 * t130 + t76 * t31 - t78 * t46;
t122 = qJD(1) * qJD(3);
t105 = t77 * t122;
t99 = pkin(4) * t105;
t2 = -t99 + t103;
t11 = t78 * t30 + t76 * t33;
t7 = -t64 * qJ(5) + t11;
t166 = t7 * t64 + t2;
t145 = t78 * t79;
t136 = t66 * t145 + t76 * t47;
t111 = t77 * t129;
t121 = qJD(3) * qJD(4);
t127 = t76 * qJD(3);
t26 = (t79 * t127 + t111) * qJD(1) + t76 * t121;
t134 = qJD(1) * t77;
t54 = t78 * t134 + t127;
t164 = t54 ^ 2;
t107 = t78 * t122;
t110 = t76 * t134;
t25 = qJD(4) * t110 - t79 * t107 - t78 * t121;
t132 = qJD(3) * t79;
t32 = qJD(3) * t70 + t60 * t132;
t3 = t26 * pkin(4) + t25 * qJ(5) - t54 * qJD(5) + t32;
t163 = t3 * t76;
t162 = t3 * t78;
t29 = -qJD(3) * pkin(3) - t165;
t126 = t78 * qJD(3);
t52 = t110 - t126;
t9 = t52 * pkin(4) - t54 * qJ(5) + t29;
t160 = t9 * t54;
t159 = t29 * t76;
t158 = t29 * t78;
t157 = t32 * t76;
t156 = t32 * t78;
t155 = t52 * t64;
t154 = t54 * t52;
t153 = t54 * t64;
t152 = t64 * t78;
t151 = t66 * t76;
t150 = t66 * t78;
t149 = t76 * t79;
t148 = t78 * t47;
t147 = t78 * t58;
t146 = t78 * t77;
t144 = t79 * t25;
t143 = t79 * t26;
t80 = qJD(3) ^ 2;
t142 = t80 * t77;
t141 = t80 * t79;
t93 = pkin(4) * t76 - qJ(5) * t78;
t140 = t76 * qJD(5) + t64 * t93 + t37;
t109 = t79 * t126;
t139 = -t52 * t109 - t26 * t146;
t57 = t97 * qJD(1);
t138 = t165 * t78 + t76 * t57;
t137 = t47 * t129 + t76 * t58;
t72 = t77 ^ 2;
t135 = -t79 ^ 2 + t72;
t61 = qJD(1) * t67;
t133 = qJD(3) * t77;
t131 = qJD(4) * t52;
t128 = t29 * qJD(4);
t10 = -t76 * t30 + t78 * t33;
t123 = qJD(5) - t10;
t120 = t76 * pkin(7) * t64;
t119 = pkin(7) * t152;
t118 = t33 * t129 + t78 * t31 + t76 * t46;
t117 = pkin(7) * t133;
t116 = pkin(7) * t126;
t115 = t54 * t132;
t114 = t64 * t130;
t113 = t77 * t130;
t112 = t64 * t129;
t108 = pkin(4) + t151;
t62 = t72 * t107;
t104 = t62 - t144;
t102 = -t25 + t131;
t101 = t64 * t111;
t100 = t54 * t111;
t98 = qJ(5) * t105;
t6 = t64 * pkin(4) + t123;
t96 = t6 * t78 - t7 * t76;
t95 = t6 * t76 + t7 * t78;
t94 = t78 * pkin(4) + t76 * qJ(5);
t92 = -t165 * t76 + t78 * t57;
t90 = 0.2e1 * qJD(3) * t61;
t89 = t66 + t93;
t88 = -t11 * t64 - t103;
t87 = t30 * t130 - t118;
t86 = (-qJD(1) * t72 + t64 * t79) * t76;
t85 = t109 - t113;
t1 = -t64 * qJD(5) - t87 + t98;
t83 = t96 * qJD(4) + t1 * t78 + t2 * t76;
t82 = qJD(3) * t86 + t52 * t133 + t101 - t143;
t81 = qJD(1) ^ 2;
t59 = -pkin(3) - t94;
t49 = t64 * t109;
t41 = t54 * t133;
t27 = t89 * t77;
t19 = t54 * pkin(4) + t52 * qJ(5);
t17 = t108 * t79 - t148;
t16 = -t79 * qJ(5) + t136;
t14 = -pkin(4) * t134 - t92;
t13 = qJ(5) * t134 + t138;
t12 = -t25 - t155;
t8 = (t94 * qJD(4) - qJD(5) * t78) * t77 + t89 * t132;
t5 = t136 * qJD(4) - t108 * t133 - t147;
t4 = (-t66 * t130 - qJD(5)) * t79 + (qJ(5) - t150) * t133 + t137;
t15 = [0, 0, 0, 0, 0.2e1 * t79 * t105, -0.2e1 * t135 * t122, t141, -t142, 0, -t66 * t141 + t77 * t90, t66 * t142 + t79 * t90, -t25 * t146 + t85 * t54, -t100 + (-t115 + (t25 + t131) * t77) * t76 + t139, t64 * t113 + t144 + t41 - t49 + t62, t101 + t143 + (-t52 * t77 + t86) * qJD(3), (-t64 - t125) * t133, -(-t47 * t130 + t147) * t64 + (t66 * t112 + (t52 * t66 + t159) * qJD(3) + t103) * t79 + (t78 * t128 + t66 * t26 + t157 + (-t64 * t151 + (-t66 * t149 + t148) * qJD(1) + t10) * qJD(3)) * t77, t137 * t64 + ((-t64 * t66 - t30) * t130 + (t54 * t66 + t158) * qJD(3) + t118) * t79 + (-t76 * t128 - t66 * t25 + t156 + (-t136 * qJD(1) - t64 * t150 - t11) * qJD(3)) * t77, t27 * t26 + t5 * t64 + t8 * t52 + (t9 * t127 + t2) * t79 + (t9 * t129 + t163 + (-qJD(1) * t17 - t6) * qJD(3)) * t77, -t16 * t26 - t17 * t25 - t4 * t52 + t5 * t54 + t96 * t132 + (-qJD(4) * t95 - t1 * t76 + t2 * t78) * t77, t27 * t25 - t4 * t64 - t8 * t54 + (-t9 * t126 - t1) * t79 + (t9 * t130 - t162 + (qJD(1) * t16 + t7) * qJD(3)) * t77, t1 * t16 + t2 * t17 + t3 * t27 + t7 * t4 + t6 * t5 + t9 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t141, 0, 0, 0, 0, 0, t82, t64 * t85 - t104 + t41, t82, t100 + (t102 * t77 + t115) * t76 + t139, -t49 + (-qJD(3) * t54 + t114) * t77 + t104, (qJD(3) * t95 - t3) * t79 + (qJD(3) * t9 + t83) * t77; 0, 0, 0, 0, -t77 * t81 * t79, t135 * t81, 0, 0, 0, t37 * qJD(3) - t61 * t134 - t32, -t61 * t125, -t54 * t152 - t25 * t76, (-t25 + t155) * t78 + (-t26 + t153) * t76, -t112 + (t64 * t145 + (-t54 + t127) * t77) * qJD(1), t114 + (-t64 * t149 + (t52 + t126) * t77) * qJD(1), t64 * t134, -pkin(3) * t26 - t156 + t92 * t64 - t37 * t52 + (t119 + t159) * qJD(4) + (-t10 * t77 + (-t29 * t79 - t117) * t76) * qJD(1), pkin(3) * t25 + t157 - t138 * t64 - t37 * t54 + (-t120 + t158) * qJD(4) + (-t29 * t145 + (t11 - t116) * t77) * qJD(1), -t14 * t64 + t59 * t26 - t162 - t140 * t52 + (t76 * t9 + t119) * qJD(4) + (t6 * t77 + (-t79 * t9 - t117) * t76) * qJD(1), t13 * t52 - t14 * t54 + (t1 - t64 * t6 + (qJD(4) * t54 - t26) * pkin(7)) * t78 + (pkin(7) * t102 + t166) * t76, t13 * t64 + t59 * t25 - t163 + t140 * t54 + (-t78 * t9 + t120) * qJD(4) + (t9 * t145 + (-t7 + t116) * t77) * qJD(1), t83 * pkin(7) - t7 * t13 - t6 * t14 - t140 * t9 + t3 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t52 ^ 2 + t164, t12, -t153 - t26, t105, -t29 * t54 + t88, -t10 * t64 + t29 * t52 + t87, -t19 * t52 - t160 + t88 + 0.2e1 * t99, pkin(4) * t25 - t26 * qJ(5) + (-t11 + t7) * t54 + (t6 - t123) * t52, 0.2e1 * t98 + t19 * t54 - t9 * t52 + (-0.2e1 * qJD(5) + t10) * t64 - t87, -t2 * pkin(4) + t1 * qJ(5) - t6 * t11 + t123 * t7 - t9 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 + t154, t12, -t64 ^ 2 - t164, t160 + t166;];
tauc_reg = t15;
