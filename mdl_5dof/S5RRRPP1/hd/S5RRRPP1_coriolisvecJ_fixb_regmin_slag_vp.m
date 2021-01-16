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
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 22:14:48
% EndTime: 2021-01-15 22:14:53
% DurationCPUTime: 0.92s
% Computational Cost: add. (2019->204), mult. (3532->258), div. (0->0), fcn. (2144->6), ass. (0->136)
t171 = -qJ(4) - pkin(7);
t109 = sin(pkin(8));
t110 = cos(pkin(8));
t111 = sin(qJ(3));
t113 = cos(qJ(3));
t133 = qJD(3) * t171;
t120 = -t111 * qJD(4) + t113 * t133;
t114 = cos(qJ(2));
t165 = pkin(1) * qJD(1);
t144 = t114 * t165;
t102 = t113 * qJD(4);
t72 = t111 * t133 + t102;
t159 = t110 * t113;
t160 = t109 * t111;
t79 = -t159 + t160;
t169 = t109 * t120 + t110 * t72 + t79 * t144;
t106 = qJD(1) + qJD(2);
t80 = t109 * t113 + t110 * t111;
t68 = t80 * t106;
t63 = t68 ^ 2;
t142 = t106 * t159;
t66 = t106 * t160 - t142;
t180 = -t66 ^ 2 - t63;
t112 = sin(qJ(2));
t97 = t112 * pkin(1) + pkin(7);
t162 = -qJ(4) - t97;
t132 = t162 * t111;
t104 = t113 * qJ(4);
t78 = t113 * t97 + t104;
t43 = t109 * t78 - t110 * t132;
t164 = pkin(1) * qJD(2);
t143 = qJD(1) * t164;
t129 = t114 * t143;
t122 = qJD(4) * t106 + t129;
t145 = t112 * t165;
t131 = -t171 * t106 + t145;
t127 = qJD(3) * t131;
t116 = -t122 * t111 - t113 * t127;
t37 = -t111 * t127 + t122 * t113;
t8 = t109 * t37 - t110 * t116;
t179 = t8 * t43;
t134 = t171 * t111;
t90 = t113 * pkin(7) + t104;
t53 = t109 * t90 - t110 * t134;
t178 = t8 * t53;
t177 = t8 * t80;
t74 = t80 * qJD(3);
t60 = t106 * t74;
t151 = qJD(3) * t113;
t137 = t110 * t151;
t152 = qJD(3) * t111;
t138 = t109 * t152;
t86 = t106 * t138;
t61 = t106 * t137 - t86;
t141 = t106 * t152;
t93 = t112 * t143;
t73 = pkin(3) * t141 + t93;
t125 = t60 * pkin(4) - t61 * qJ(5) + t73;
t11 = -t68 * qJD(5) + t125;
t99 = -t113 * pkin(3) - pkin(2);
t65 = t99 * t106 + qJD(4) - t144;
t20 = t66 * pkin(4) - t68 * qJ(5) + t65;
t176 = t11 * t79 + t20 * t74;
t75 = t137 - t138;
t175 = -t11 * t80 - t20 * t75;
t174 = pkin(3) * t111;
t173 = t114 * pkin(1);
t172 = t20 * t68;
t9 = t109 * t116 + t110 * t37;
t170 = t109 * t72 - t110 * t120 - t80 * t144;
t168 = t65 * t74 + t73 * t79;
t167 = t65 * t75 + t73 * t80;
t57 = t131 * t113;
t51 = t110 * t57;
t56 = t131 * t111;
t55 = qJD(3) * pkin(3) - t56;
t27 = t109 * t55 + t51;
t85 = -t106 * pkin(2) - t144;
t166 = t111 * t93 + t85 * t151;
t163 = t109 * t57;
t161 = t106 * t111;
t158 = t112 * t113;
t115 = qJD(3) ^ 2;
t157 = t115 * t111;
t103 = t115 * t113;
t130 = qJD(3) * t162;
t146 = t114 * t164;
t117 = (-qJD(4) - t146) * t111 + t113 * t130;
t49 = t111 * t130 + t113 * t146 + t102;
t18 = t109 * t49 - t110 * t117;
t156 = t18 * qJD(3);
t19 = t109 * t117 + t110 * t49;
t155 = t19 * qJD(3);
t30 = -t110 * t56 - t163;
t154 = qJD(5) - t30;
t153 = t111 ^ 2 - t113 ^ 2;
t150 = qJD(3) * t114;
t149 = -qJD(1) - t106;
t148 = -qJD(2) + t106;
t147 = pkin(3) * t161;
t101 = t112 * t164;
t100 = pkin(3) * t152;
t140 = t106 * t151;
t139 = t111 * t150;
t26 = t110 * t55 - t163;
t23 = -qJD(3) * pkin(4) + qJD(5) - t26;
t24 = qJD(3) * qJ(5) + t27;
t7 = qJD(3) * qJD(5) + t9;
t136 = t23 * t75 - t24 * t74 - t7 * t79 + t177;
t135 = -t26 * t75 - t27 * t74 - t9 * t79 + t177;
t28 = t74 * pkin(4) - t75 * qJ(5) - t80 * qJD(5) + t100;
t128 = -t28 + t145;
t29 = -t109 * t56 + t51;
t126 = t29 * qJD(3) - t8;
t44 = t109 * t132 + t110 * t78;
t124 = t18 * t68 - t19 * t66 + t43 * t61 - t44 * t60;
t123 = -t106 * t85 - t129;
t121 = -t112 * t161 + t113 * t150;
t45 = t79 * pkin(4) - t80 * qJ(5) + t99;
t119 = 0.2e1 * t68 * qJD(3);
t54 = t109 * t134 + t110 * t90;
t118 = -t169 * t66 + t170 * t68 + t53 * t61 - t54 * t60;
t105 = t106 ^ 2;
t98 = -pkin(2) - t173;
t96 = -t110 * pkin(3) - pkin(4);
t94 = t109 * pkin(3) + qJ(5);
t87 = t99 - t173;
t83 = 0.2e1 * t111 * t140;
t82 = t101 + t100;
t76 = t85 * t152;
t70 = -0.2e1 * t153 * t106 * qJD(3);
t40 = t45 - t173;
t36 = -t86 + (-t66 + t142) * qJD(3);
t31 = t68 * pkin(4) + t66 * qJ(5) + t147;
t21 = t101 + t28;
t1 = [0, 0, 0, 0, -t106 * t101 - t93, t149 * t146, t83, t70, t103, -t157, 0, t98 * t141 - t97 * t103 + t76 + (t149 * t158 - t139) * t164, -t121 * t164 + t98 * t140 + t97 * t157 + t166, t87 * t60 + t82 * t66 - t156 + t168, t87 * t61 + t82 * t68 - t155 + t167, t124 + t135, -t26 * t18 + t27 * t19 + t9 * t44 + t65 * t82 + t73 * t87 + t179, t21 * t66 + t40 * t60 - t156 + t176, t124 + t136, -t21 * t68 - t40 * t61 + t155 + t175, t11 * t40 + t23 * t18 + t24 * t19 + t20 * t21 + t7 * t44 + t179; 0, 0, 0, 0, t106 * t145 - t93, t148 * t144, t83, t70, t103, -t157, 0, -pkin(2) * t141 - pkin(7) * t103 + t76 + (t148 * t158 + t139) * t165, -pkin(2) * t140 + pkin(7) * t157 + t121 * t165 + t166, -t66 * t145 + t99 * t60 + (t66 * t174 - t170) * qJD(3) + t168, -t68 * t145 + t99 * t61 + (t68 * t174 - t169) * qJD(3) + t167, t118 + t135, t178 + t9 * t54 + t73 * t99 + (-t145 + t100) * t65 + t169 * t27 - t170 * t26, -t170 * qJD(3) - t128 * t66 + t45 * t60 + t176, t118 + t136, t169 * qJD(3) + t128 * t68 - t45 * t61 + t175, t11 * t45 - t128 * t20 + t169 * t24 + t170 * t23 + t7 * t54 + t178; 0, 0, 0, 0, 0, 0, -t111 * t105 * t113, t153 * t105, 0, 0, 0, t123 * t111, t123 * t113, -t66 * t147 - t65 * t68 + t126, t30 * qJD(3) - t68 * t147 + t65 * t66 - t9, (t27 - t29) * t68 + (-t26 + t30) * t66 + (-t109 * t60 - t110 * t61) * pkin(3), t26 * t29 - t27 * t30 + (t109 * t9 - t110 * t8 - t65 * t161) * pkin(3), -t31 * t66 + t126 - t172, -t94 * t60 + t96 * t61 + (t24 - t29) * t68 + (t23 - t154) * t66, -t20 * t66 + t31 * t68 + (0.2e1 * qJD(5) - t30) * qJD(3) + t9, t154 * t24 - t20 * t31 - t23 * t29 + t7 * t94 + t8 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t36, t180, t26 * t68 + t27 * t66 + t73, t119, t180, -t36, t24 * t66 + (-qJD(5) - t23) * t68 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t66, -t86 + (t66 + t142) * qJD(3), -t63 - t115, -t24 * qJD(3) + t172 + t8;];
tauc_reg = t1;
