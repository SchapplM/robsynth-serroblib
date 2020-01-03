% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:50:01
% DurationCPUTime: 0.88s
% Computational Cost: add. (1803->182), mult. (3144->237), div. (0->0), fcn. (1904->6), ass. (0->126)
t167 = 2 * qJD(3);
t166 = (qJ(4) + pkin(7));
t101 = sin(pkin(8));
t102 = cos(pkin(8));
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t123 = qJD(3) * t166;
t111 = -t103 * qJD(4) - t105 * t123;
t106 = cos(qJ(2));
t154 = pkin(1) * qJD(1);
t135 = t106 * t154;
t94 = t105 * qJD(4);
t67 = -t103 * t123 + t94;
t147 = t102 * t105;
t148 = t101 * t103;
t74 = -t147 + t148;
t156 = t101 * t111 + t102 * t67 + t74 * t135;
t75 = t101 * t105 + t102 * t103;
t98 = qJD(1) + qJD(2);
t63 = t75 * t98;
t58 = t63 ^ 2;
t137 = t98 * t147;
t61 = t98 * t148 - t137;
t165 = -t61 ^ 2 - t58;
t104 = sin(qJ(2));
t90 = t104 * pkin(1) + pkin(7);
t150 = -qJ(4) - t90;
t122 = t150 * t103;
t96 = t105 * qJ(4);
t73 = t105 * t90 + t96;
t40 = t101 * t73 - t102 * t122;
t153 = pkin(1) * qJD(2);
t134 = qJD(1) * t153;
t120 = t106 * t134;
t113 = qJD(4) * t98 + t120;
t136 = t104 * t154;
t124 = t166 * t98 + t136;
t118 = qJD(3) * t124;
t108 = -t113 * t103 - t105 * t118;
t36 = -t103 * t118 + t113 * t105;
t8 = t101 * t36 - t102 * t108;
t164 = t8 * t40;
t125 = t166 * t103;
t83 = t105 * pkin(7) + t96;
t48 = t101 * t83 + t102 * t125;
t163 = t8 * t48;
t162 = t8 * t75;
t69 = t75 * qJD(3);
t55 = t98 * t69;
t140 = qJD(3) * t105;
t128 = t102 * t140;
t141 = qJD(3) * t103;
t129 = t101 * t141;
t80 = t98 * t129;
t56 = t98 * t128 - t80;
t133 = t98 * t141;
t86 = t104 * t134;
t68 = pkin(3) * t133 + t86;
t116 = t55 * pkin(4) - t56 * qJ(5) + t68;
t11 = -t63 * qJD(5) + t116;
t131 = -t105 * pkin(3) - pkin(2);
t60 = t131 * t98 + qJD(4) - t135;
t20 = t61 * pkin(4) - t63 * qJ(5) + t60;
t161 = t11 * t74 + t20 * t69;
t70 = t128 - t129;
t160 = -t11 * t75 - t20 * t70;
t159 = t106 * pkin(1);
t9 = t101 * t108 + t102 * t36;
t157 = t101 * t67 - t102 * t111 - t75 * t135;
t52 = t124 * t105;
t46 = t102 * t52;
t51 = t124 * t103;
t50 = qJD(3) * pkin(3) - t51;
t27 = t101 * t50 + t46;
t79 = -t98 * pkin(2) - t135;
t155 = t103 * t86 + t79 * t140;
t152 = t101 * t52;
t151 = t103 * t98;
t149 = -t103 ^ 2 + t105 ^ 2;
t146 = t104 * t105;
t107 = qJD(3) ^ 2;
t145 = t107 * t103;
t95 = t107 * t105;
t144 = -qJD(1) - t98;
t143 = -qJD(2) + t98;
t30 = -t102 * t51 - t152;
t142 = qJD(5) - t30;
t139 = qJD(3) * t106;
t138 = t106 * t153;
t93 = t104 * t153;
t92 = pkin(3) * t141;
t132 = t98 * t140;
t130 = t103 * t139;
t26 = t102 * t50 - t152;
t23 = -qJD(3) * pkin(4) + qJD(5) - t26;
t24 = qJD(3) * qJ(5) + t27;
t7 = qJD(3) * qJD(5) + t9;
t127 = t23 * t70 - t24 * t69 - t7 * t74 + t162;
t126 = -t26 * t70 - t27 * t69 - t9 * t74 + t162;
t121 = qJD(3) * t150;
t28 = t69 * pkin(4) - t70 * qJ(5) - t75 * qJD(5) + t92;
t119 = -t28 + t136;
t117 = t20 * t63 + t8;
t109 = (-qJD(4) - t138) * t103 + t105 * t121;
t44 = t103 * t121 + t105 * t138 + t94;
t18 = t101 * t44 - t102 * t109;
t19 = t101 * t109 + t102 * t44;
t41 = t101 * t122 + t102 * t73;
t115 = t18 * t63 - t19 * t61 + t40 * t56 - t41 * t55;
t114 = -t79 * t98 - t120;
t112 = -t104 * t151 + t105 * t139;
t42 = t74 * pkin(4) - t75 * qJ(5) + t131;
t49 = -t101 * t125 + t102 * t83;
t110 = -t156 * t61 + t157 * t63 + t48 * t56 - t49 * t55;
t97 = t98 ^ 2;
t91 = -pkin(2) - t159;
t89 = -t102 * pkin(3) - pkin(4);
t87 = t101 * pkin(3) + qJ(5);
t77 = 0.2e1 * t103 * t132;
t71 = t79 * t141;
t65 = t149 * t98 * t167;
t39 = t42 - t159;
t31 = pkin(3) * t151 + t63 * pkin(4) + t61 * qJ(5);
t29 = -t101 * t51 + t46;
t21 = t28 + t93;
t1 = [0, 0, 0, 0, -t98 * t93 - t86, t144 * t138, t77, t65, t95, -t145, 0, t91 * t133 - t90 * t95 + t71 + (t144 * t146 - t130) * t153, -t112 * t153 + t91 * t132 + t90 * t145 + t155, t115 + t126, t9 * t41 + t27 * t19 + t164 - t26 * t18 + t68 * (t131 - t159) + t60 * (t93 + t92), -t18 * qJD(3) + t21 * t61 + t39 * t55 + t161, t115 + t127, t19 * qJD(3) - t21 * t63 - t39 * t56 + t160, t11 * t39 + t23 * t18 + t24 * t19 + t20 * t21 + t7 * t41 + t164; 0, 0, 0, 0, t98 * t136 - t86, t143 * t135, t77, t65, t95, -t145, 0, -pkin(2) * t133 - pkin(7) * t95 + t71 + (t143 * t146 + t130) * t154, -pkin(2) * t132 + pkin(7) * t145 + t112 * t154 + t155, t110 + t126, t9 * t49 + t163 + t68 * t131 + (-t136 + t92) * t60 + t156 * t27 - t157 * t26, -t157 * qJD(3) - t119 * t61 + t42 * t55 + t161, t110 + t127, t156 * qJD(3) + t119 * t63 - t42 * t56 + t160, t11 * t42 - t119 * t20 + t156 * t24 + t157 * t23 + t7 * t49 + t163; 0, 0, 0, 0, 0, 0, -t103 * t97 * t105, -t149 * t97, 0, 0, 0, t114 * t103, t114 * t105, (t27 - t29) * t63 + (-t26 + t30) * t61 + (-t101 * t55 - t102 * t56) * pkin(3), t26 * t29 - t27 * t30 + (t101 * t9 - t102 * t8 - t60 * t151) * pkin(3), t29 * qJD(3) - t31 * t61 - t117, -t87 * t55 + t89 * t56 + (t24 - t29) * t63 + (t23 - t142) * t61, -t20 * t61 + t31 * t63 + (0.2e1 * qJD(5) - t30) * qJD(3) + t9, t142 * t24 - t20 * t31 - t23 * t29 + t7 * t87 + t8 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t26 * t63 + t27 * t61 + t68, t63 * t167, t165, t80 + (t61 - t137) * qJD(3), t24 * t61 + (-qJD(5) - t23) * t63 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t80 + (t61 + t137) * qJD(3), -t58 - t107, -t24 * qJD(3) + t117;];
tauc_reg = t1;
