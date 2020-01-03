% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR15_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:29
% DurationCPUTime: 1.44s
% Computational Cost: add. (1234->220), mult. (2771->344), div. (0->0), fcn. (1799->6), ass. (0->127)
t90 = sin(pkin(8));
t91 = cos(pkin(8));
t92 = sin(qJ(5));
t94 = cos(qJ(5));
t164 = -t92 * t90 + t94 * t91;
t166 = t164 * qJD(5);
t95 = cos(qJ(3));
t146 = qJD(1) * t95;
t129 = t90 * t146;
t139 = t91 * qJD(3);
t63 = t129 - t139;
t128 = t91 * t146;
t140 = t90 * qJD(3);
t65 = t128 + t140;
t109 = t92 * t63 - t94 * t65;
t93 = sin(qJ(3));
t138 = t93 * qJD(1);
t84 = qJD(5) + t138;
t168 = t109 * t84;
t21 = t94 * t63 + t92 * t65;
t133 = 0.2e1 * qJD(1);
t96 = -pkin(1) - pkin(6);
t167 = qJD(1) * t96;
t68 = t94 * t90 + t92 * t91;
t103 = t68 * qJD(5);
t165 = -qJD(5) + t84;
t135 = qJD(1) * qJD(3);
t123 = t93 * t135;
t117 = t94 * t123;
t118 = t92 * t123;
t6 = -t109 * qJD(5) - t90 * t117 - t91 * t118;
t163 = t91 * pkin(7);
t162 = t21 * t84;
t83 = qJD(2) + t167;
t74 = t93 * t83;
t159 = t93 * t96;
t157 = t95 * t83;
t97 = qJD(3) ^ 2;
t156 = t97 * t93;
t155 = t97 * t95;
t154 = pkin(7) + qJ(4);
t115 = pkin(3) * t95 + qJ(4) * t93;
t47 = t115 * qJD(3) - t95 * qJD(4) + qJD(2);
t32 = t47 * qJD(1);
t49 = (qJD(4) + t157) * qJD(3);
t11 = t90 * t32 + t91 * t49;
t75 = t93 * pkin(3) - t95 * qJ(4) + qJ(2);
t57 = t75 * qJD(1);
t136 = qJD(3) * qJ(4);
t58 = t74 + t136;
t19 = t90 * t57 + t91 * t58;
t107 = t164 * t93;
t153 = qJD(1) * t107 + t166;
t104 = qJD(1) * t68;
t152 = t93 * t104 + t103;
t70 = t115 * qJD(1);
t29 = t91 * t157 + t90 * t70;
t142 = qJD(3) * t96;
t130 = t95 * t142;
t25 = t91 * t130 + t90 * t47;
t34 = t91 * t159 + t90 * t75;
t150 = t93 ^ 2 - t95 ^ 2;
t98 = qJD(1) ^ 2;
t149 = -t97 - t98;
t148 = qJD(3) * pkin(3);
t147 = t98 * qJ(2);
t145 = qJD(3) * t84;
t144 = qJD(3) * t93;
t143 = qJD(3) * t95;
t137 = qJ(2) * qJD(3);
t134 = t93 * t163;
t132 = t90 * t138;
t131 = t93 * t140;
t127 = qJD(2) * t133;
t10 = t91 * t32 - t90 * t49;
t100 = (pkin(4) * t95 + t134) * qJD(1);
t4 = qJD(3) * t100 + t10;
t119 = t90 * t123;
t8 = pkin(7) * t119 + t11;
t126 = t94 * t4 - t92 * t8;
t125 = -t90 * t96 + pkin(4);
t124 = pkin(4) * t90 - t96;
t122 = t95 * t135;
t18 = t91 * t57 - t90 * t58;
t121 = -t65 + t140;
t120 = -qJD(4) + t148;
t28 = -t90 * t157 + t91 * t70;
t116 = t92 * t4 + t94 * t8;
t7 = pkin(4) * t138 - t65 * pkin(7) + t18;
t9 = -t63 * pkin(7) + t19;
t1 = t94 * t7 - t92 * t9;
t2 = t92 * t7 + t94 * t9;
t114 = -t10 * t90 + t11 * t91;
t113 = -t18 * t91 - t19 * t90;
t112 = -t18 * t90 + t19 * t91;
t61 = t91 * t75;
t20 = t125 * t93 - t95 * t163 + t61;
t27 = -t90 * t95 * pkin(7) + t34;
t111 = t94 * t20 - t92 * t27;
t110 = t92 * t20 + t94 * t27;
t108 = (-t63 - t139) * t93;
t51 = -t120 - t157;
t80 = t154 * t91;
t106 = qJD(4) * t90 + qJD(5) * t80 + t100 + t28;
t79 = t154 * t90;
t105 = pkin(7) * t132 - qJD(4) * t91 + qJD(5) * t79 + t29;
t102 = t164 * qJD(1);
t101 = -t51 + (t83 + t167) * t95;
t99 = -t95 * t136 + (t120 + t51) * t93;
t5 = -t21 * qJD(5) - t91 * t117 + t90 * t118;
t86 = -t91 * pkin(4) - pkin(3);
t69 = t83 * t144;
t62 = t124 * t95;
t52 = t124 * t144;
t46 = t164 * t95;
t45 = t68 * t95;
t40 = -pkin(4) * t132 + t74;
t37 = t91 * t47;
t35 = -pkin(4) * t119 + t69;
t33 = -t90 * t159 + t61;
t26 = t63 * pkin(4) + t51;
t24 = -t90 * t130 + t37;
t16 = pkin(7) * t131 + t25;
t15 = -t92 * t93 * t139 - t94 * t131 + t95 * t166;
t14 = -qJD(3) * t107 - t95 * t103;
t12 = t37 + (t125 * t95 + t134) * qJD(3);
t3 = [0, 0, 0, 0, t127, qJ(2) * t127, -0.2e1 * t93 * t122, 0.2e1 * t150 * t135, -t156, -t155, 0, -t96 * t156 + (qJD(2) * t93 + t95 * t137) * t133, -t96 * t155 + (qJD(2) * t95 - t93 * t137) * t133, (qJD(1) * t24 + t10) * t93 + ((qJD(1) * t33 + t18) * t95 + (t101 * t90 + t63 * t96) * t93) * qJD(3), (-qJD(1) * t25 - t11) * t93 + ((-qJD(1) * t34 - t19) * t95 + (t101 * t91 + t65 * t96) * t93) * qJD(3), -t24 * t65 - t25 * t63 + (-t10 * t91 - t11 * t90) * t95 + ((t33 * t91 + t34 * t90) * qJD(1) - t113) * t144, t10 * t33 + t11 * t34 + t18 * t24 + t19 * t25 + (t51 - t157) * t93 * t142, -t109 * t14 + t5 * t46, t109 * t15 - t14 * t21 - t5 * t45 - t46 * t6, t14 * t84 + t5 * t93 + (qJD(1) * t46 - t109) * t143, -t15 * t84 - t6 * t93 + (-qJD(1) * t45 - t21) * t143, (t84 + t138) * t143, (t94 * t12 - t92 * t16) * t84 + t126 * t93 - t52 * t21 + t62 * t6 + t35 * t45 + t26 * t15 + (-t110 * t84 - t2 * t93) * qJD(5) + (t111 * qJD(1) + t1) * t143, -(t92 * t12 + t94 * t16) * t84 - t116 * t93 + t52 * t109 + t62 * t5 + t35 * t46 + t26 * t14 + (-t1 * t93 - t111 * t84) * qJD(5) + (-qJD(1) * t110 - t2) * t143; 0, 0, 0, 0, -t98, -t147, 0, 0, 0, 0, 0, t149 * t93, t149 * t95, (-t91 * t98 + (t63 - t129) * qJD(3)) * t93, (t90 * t98 + (t65 - t128) * qJD(3)) * t93, (-t63 * t91 + t65 * t90) * t143 + (t63 * t90 + t65 * t91) * qJD(1), t114 * t93 + t113 * qJD(1) + (t51 * t93 + (t112 - t74) * t95) * qJD(3), 0, 0, 0, 0, 0, -t84 * t102 + (-t68 * t145 - t6) * t95 + (-t84 * t166 + (-t68 * t146 + t21) * qJD(3)) * t93, t84 * t104 + (-t145 * t164 - t5) * t95 + (t84 * t103 + (-t102 * t95 - t109) * qJD(3)) * t93; 0, 0, 0, 0, 0, 0, t95 * t98 * t93, -t150 * t98, 0, 0, 0, -t95 * t147, t93 * t147, t83 * t108 + (-t18 * t95 - t28 * t93 + t99 * t90) * qJD(1), t121 * t74 + (t19 * t95 + t29 * t93 + t99 * t91) * qJD(1), t28 * t65 + t29 * t63 + (-qJD(4) * t63 - t18 * t138 + t11) * t91 + (qJD(4) * t65 - t19 * t138 - t10) * t90, -t18 * t28 - t19 * t29 + (-t51 - t148) * t74 + t112 * qJD(4) + t114 * qJ(4), -t109 * t153 + t5 * t68, t109 * t152 - t153 * t21 + t164 * t5 - t68 * t6, t153 * t84 + (qJD(3) * t68 + t109) * t146, -t152 * t84 + (qJD(3) * t164 + t21) * t146, -t84 * t146, -t40 * t21 - t35 * t164 + t86 * t6 + (t105 * t92 - t106 * t94) * t84 + t152 * t26 + ((-t94 * t79 - t92 * t80) * qJD(3) - t1) * t146, t40 * t109 + t35 * t68 + t86 * t5 + (t105 * t94 + t106 * t92) * t84 + t153 * t26 + (-(-t92 * t79 + t94 * t80) * qJD(3) + t2) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t138, qJD(1) * t108, -t63 ^ 2 - t65 ^ 2, t18 * t65 + t19 * t63 + t69, 0, 0, 0, 0, 0, t6 - t168, t5 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t21, t109 ^ 2 - t21 ^ 2, t5 + t162, -t6 - t168, t122, t109 * t26 + t165 * t2 + t126, t165 * t1 + t26 * t21 - t116;];
tauc_reg = t3;
