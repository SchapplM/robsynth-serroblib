% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR13_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:53
% EndTime: 2019-12-31 18:32:57
% DurationCPUTime: 1.08s
% Computational Cost: add. (1241->203), mult. (3330->264), div. (0->0), fcn. (2442->6), ass. (0->114)
t142 = cos(qJ(3));
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t87 = sin(qJ(3));
t66 = t142 * t84 + t87 * t85;
t151 = t66 * qJD(1);
t152 = qJD(5) + t151;
t86 = sin(qJ(5));
t121 = t86 * qJD(3);
t115 = t142 * t85;
t107 = qJD(1) * t115;
t133 = t87 * t84;
t116 = qJD(1) * t133;
t57 = -t107 + t116;
t88 = cos(qJ(5));
t37 = -t88 * t57 + t121;
t154 = t152 * t37;
t153 = t151 * qJD(3);
t132 = pkin(6) + qJ(2);
t70 = t132 * t84;
t67 = qJD(1) * t70;
t71 = t132 * t85;
t68 = qJD(1) * t71;
t131 = -t142 * t67 - t87 * t68;
t149 = qJD(4) - t131;
t110 = t152 ^ 2;
t76 = qJD(3) * t107;
t45 = qJD(3) * t116 - t76;
t41 = t88 * t45;
t99 = -t86 * t110 - t41;
t150 = -qJD(5) + t152;
t111 = qJD(3) * t142;
t112 = qJD(2) * t142;
t24 = (qJD(2) * t84 + qJD(3) * t71) * t87 + t70 * t111 - t85 * t112;
t143 = t57 * pkin(4);
t145 = pkin(3) + pkin(7);
t33 = t142 * t68 - t87 * t67;
t29 = -qJD(3) * qJ(4) - t33;
t15 = -t29 - t143;
t148 = t145 * t45 + (t15 - t33 + t143) * t152;
t147 = t57 ^ 2;
t146 = t151 ^ 2;
t62 = t66 * qJD(3);
t46 = qJD(1) * t62;
t144 = t46 * pkin(3);
t81 = -t85 * pkin(2) - pkin(1);
t100 = -t66 * qJ(4) + t81;
t65 = -t115 + t133;
t18 = t145 * t65 + t100;
t141 = t18 * t45;
t124 = qJD(5) * t88;
t21 = -qJD(5) * t121 + t57 * t124 + t86 * t46;
t140 = t21 * t88;
t69 = t81 * qJD(1) + qJD(2);
t91 = -qJ(4) * t151 + t69;
t23 = t57 * pkin(3) + t91;
t139 = t23 * t151;
t39 = t88 * qJD(3) + t86 * t57;
t138 = t39 * t57;
t137 = t57 * t37;
t136 = t151 * t57;
t135 = t65 * t86;
t134 = t86 * t45;
t129 = t84 ^ 2 + t85 ^ 2;
t128 = t45 * qJ(4);
t127 = t57 * qJ(4);
t126 = qJD(3) * t87;
t125 = qJD(5) * t65;
t123 = t24 * qJD(3);
t97 = t142 * t71 - t87 * t70;
t25 = t66 * qJD(2) + t97 * qJD(3);
t122 = t25 * qJD(3);
t119 = pkin(4) * t151 + t149;
t117 = qJD(1) * qJD(2);
t114 = t129 * qJD(1) ^ 2;
t113 = t87 * t117;
t109 = t152 * t39;
t106 = qJD(1) * t112;
t108 = t85 * t106 - t67 * t111 - t84 * t113 - t68 * t126;
t14 = t84 * t106 + t68 * t111 + t85 * t113 - t67 * t126;
t35 = t142 * t70 + t87 * t71;
t12 = -t145 * qJD(3) + t119;
t8 = t145 * t57 + t91;
t1 = t88 * t12 - t86 * t8;
t2 = t86 * t12 + t88 * t8;
t104 = -qJD(4) * t151 + t128;
t61 = -t85 * t111 + t84 * t126;
t103 = t61 * qJ(4) - t66 * qJD(4);
t102 = 0.2e1 * t129 * t117;
t98 = t65 * t124 + t86 * t62;
t96 = t33 * qJD(3) - t14;
t13 = -qJD(3) * qJD(4) - t108;
t4 = -t46 * pkin(4) - t13;
t95 = t4 + (t152 * t145 + t127) * t152;
t26 = t66 * pkin(4) + t35;
t94 = t15 * t62 + t26 * t45 + t4 * t65;
t93 = -t88 * t110 + t134;
t47 = qJD(3) * t57;
t42 = t88 * t46;
t34 = t45 * t66;
t31 = t65 * pkin(3) + t100;
t30 = pkin(3) * t151 + t127;
t28 = -qJD(3) * pkin(3) + t149;
t27 = -t65 * pkin(4) + t97;
t22 = t39 * qJD(5) - t42;
t17 = t62 * pkin(3) + t103;
t11 = t104 + t144;
t10 = -t61 * pkin(4) + t25;
t9 = -t62 * pkin(4) - t24;
t7 = t145 * t62 + t103;
t6 = -t45 * pkin(4) + t14;
t5 = t88 * t6;
t3 = t145 * t46 + t104;
t16 = [0, 0, 0, 0, 0, t102, qJ(2) * t102, -t151 * t61 - t34, -t151 * t62 + t45 * t65 - t66 * t46 + t61 * t57, -t61 * qJD(3), -t62 * qJD(3), 0, t81 * t46 + t69 * t62 - t122, -t81 * t45 - t69 * t61 + t123, t13 * t65 + t14 * t66 + t151 * t25 + t24 * t57 - t28 * t61 + t29 * t62 - t35 * t45 - t46 * t97, -t11 * t65 - t17 * t57 - t23 * t62 - t31 * t46 + t122, -t11 * t66 - t151 * t17 + t23 * t61 + t31 * t45 - t123, t11 * t31 - t13 * t97 + t14 * t35 + t23 * t17 + t29 * t24 + t28 * t25, t21 * t135 + t98 * t39, (-t37 * t86 + t39 * t88) * t62 + (t140 - t22 * t86 + (-t37 * t88 - t39 * t86) * qJD(5)) * t65, -t65 * t134 + t152 * t98 + t21 * t66 - t39 * t61, -t65 * t41 - t22 * t66 + t37 * t61 + (-t86 * t125 + t88 * t62) * t152, -t152 * t61 - t34, -t1 * t61 + t27 * t22 + t9 * t37 + t5 * t66 + (-t152 * t7 - t3 * t66 + t141) * t86 + (t10 * t152 - t94) * t88 + ((-t88 * t18 - t86 * t26) * t152 - t2 * t66 + t15 * t135) * qJD(5), t2 * t61 + t27 * t21 + t9 * t39 + (-(qJD(5) * t26 + t7) * t152 + t141 - (qJD(5) * t12 + t3) * t66 + t15 * t125) * t88 + (-(-qJD(5) * t18 + t10) * t152 - (-qJD(5) * t8 + t6) * t66 + t94) * t86; 0, 0, 0, 0, 0, -t114, -qJ(2) * t114, 0, 0, 0, 0, 0, 0.2e1 * t153, t76 + (-t57 - t116) * qJD(3), -t146 - t147, -0.2e1 * t153, t45 + t47, t144 + t128 - t29 * t57 + (-qJD(4) - t28) * t151, 0, 0, 0, 0, 0, t93 + t137, t138 - t99; 0, 0, 0, 0, 0, 0, 0, t136, t146 - t147, t76 + (t57 - t116) * qJD(3), 0, 0, -t151 * t69 + t96, qJD(3) * t131 + t69 * t57 - t108, pkin(3) * t45 - qJ(4) * t46 + (-t29 - t33) * t151 + (t28 - t149) * t57, t30 * t57 + t139 - t96, -t23 * t57 + t30 * t151 + (0.2e1 * qJD(4) - t131) * qJD(3) + t108, -t14 * pkin(3) - t13 * qJ(4) - t149 * t29 - t23 * t30 - t28 * t33, -t86 * t109 + t140, (-t22 - t109) * t88 + (-t21 + t154) * t86, t99 + t138, t93 - t137, t152 * t57, qJ(4) * t22 + t1 * t57 + t119 * t37 + t148 * t88 + t95 * t86, qJ(4) * t21 + t119 * t39 - t148 * t86 - t2 * t57 + t95 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 + t47, -t136, -qJD(3) ^ 2 - t146, t29 * qJD(3) + t139 + t14, 0, 0, 0, 0, 0, -qJD(3) * t37 + t99, -qJD(3) * t39 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t21 + t154, t150 * t39 + t42, -t45, -t15 * t39 + t150 * t2 - t86 * t3 + t5, t150 * t1 + t15 * t37 - t88 * t3 - t86 * t6;];
tauc_reg = t16;
