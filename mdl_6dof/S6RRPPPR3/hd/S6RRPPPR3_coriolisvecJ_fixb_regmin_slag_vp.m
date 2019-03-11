% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:38
% EndTime: 2019-03-09 08:15:44
% DurationCPUTime: 2.27s
% Computational Cost: add. (1991->327), mult. (4527->445), div. (0->0), fcn. (2645->6), ass. (0->189)
t133 = sin(qJ(2));
t196 = qJD(1) * t133;
t104 = qJD(6) + t196;
t127 = sin(pkin(9));
t132 = sin(qJ(6));
t128 = cos(pkin(9));
t134 = cos(qJ(6));
t207 = t134 * t128;
t84 = t127 * t132 - t207;
t146 = t84 * t133;
t76 = t84 * qJD(6);
t219 = qJD(1) * t146 + t76;
t227 = t219 * t104;
t226 = qJD(6) - t104;
t135 = cos(qJ(2));
t195 = qJD(1) * t135;
t188 = qJD(1) * qJD(2);
t225 = -0.2e1 * t188;
t124 = qJD(2) * qJ(3);
t112 = pkin(7) * t195;
t91 = -qJ(4) * t195 + t112;
t79 = -t124 - t91;
t179 = t127 * t195;
t191 = t128 * qJD(2);
t82 = t179 - t191;
t178 = t128 * t195;
t194 = qJD(2) * t127;
t83 = t178 + t194;
t40 = t132 * t82 - t134 * t83;
t224 = qJD(6) * t40;
t85 = t127 * t134 + t128 * t132;
t58 = t85 * t196;
t77 = t85 * qJD(6);
t218 = -t77 - t58;
t176 = t218 * t104;
t136 = -pkin(2) - pkin(3);
t180 = qJD(2) * t136;
t189 = t135 * qJD(4);
t193 = qJD(2) * t133;
t149 = pkin(7) * t193 + t189;
t175 = t133 * t188;
t101 = qJ(4) * t175;
t123 = qJD(2) * qJD(3);
t200 = t101 + t123;
t51 = qJD(1) * t149 - t200;
t118 = t135 * qJ(4);
t98 = pkin(7) * t135 - t118;
t223 = t51 * t98;
t222 = pkin(7) - qJ(4);
t122 = -qJ(5) + t136;
t221 = pkin(8) - t122;
t131 = qJ(3) + pkin(4);
t145 = pkin(4) * t135 + t122 * t133;
t140 = qJD(2) * t145 + qJD(5) * t135;
t116 = t133 * qJD(3);
t174 = t135 * t188;
t199 = qJ(3) * t174 + qJD(1) * t116;
t25 = qJD(1) * t140 + t199;
t103 = pkin(7) * t174;
t192 = qJD(2) * t135;
t177 = qJ(4) * t192;
t190 = t133 * qJD(4);
t60 = t103 + (-t177 - t190) * qJD(1);
t52 = -qJD(2) * qJD(5) + t60;
t8 = t127 * t25 + t128 * t52;
t157 = pkin(4) * t133 + qJ(5) * t135;
t80 = -qJD(1) * pkin(1) - pkin(2) * t195 - qJ(3) * t196;
t57 = pkin(3) * t195 + qJD(4) - t80;
t41 = qJD(1) * t157 + t57;
t107 = qJ(4) * t196;
t186 = pkin(7) * t196;
t162 = qJD(3) + t186;
t151 = -t107 + t162;
t55 = t122 * qJD(2) + t151;
t13 = t127 * t41 + t128 * t55;
t198 = qJ(3) * t192 + t116;
t33 = t140 + t198;
t72 = t222 * t192 - t190;
t16 = t127 * t33 + t128 * t72;
t74 = t134 * t82;
t220 = qJD(6) * t74 + t175 * t207;
t109 = qJ(3) * t195;
t44 = qJD(1) * t145 + t109;
t24 = t127 * t44 + t128 * t91;
t95 = -t135 * pkin(2) - t133 * qJ(3) - pkin(1);
t81 = t135 * pkin(3) - t95;
t54 = t157 + t81;
t97 = t222 * t133;
t32 = t127 * t54 + t128 * t97;
t217 = qJD(2) * pkin(2);
t216 = t104 * t40;
t12 = -t127 * t55 + t128 * t41;
t215 = t12 * t135;
t214 = t13 * t135;
t213 = qJD(2) * t84;
t125 = t133 ^ 2;
t138 = qJD(1) ^ 2;
t212 = t125 * t138;
t211 = t127 * t133;
t210 = t127 * t135;
t209 = t128 * t133;
t208 = t128 * t135;
t206 = t135 * t138;
t137 = qJD(2) ^ 2;
t205 = t137 * t133;
t204 = t137 * t135;
t182 = pkin(5) * t127 - pkin(7);
t163 = t182 * t133;
t203 = -qJD(1) * t163 + qJD(3) - t107;
t89 = -t107 + t186;
t202 = qJD(3) + t89;
t201 = -qJD(4) - t57;
t126 = t135 ^ 2;
t197 = t125 - t126;
t187 = pkin(8) * t211;
t185 = t127 * t212;
t184 = t128 * t212;
t183 = t133 * t206;
t152 = pkin(5) * t135 - pkin(8) * t209;
t144 = t152 * qJD(2);
t7 = -t127 * t52 + t128 * t25;
t4 = qJD(1) * t144 + t7;
t164 = t127 * t175;
t5 = -pkin(8) * t164 + t8;
t181 = -t132 * t5 + t134 * t4;
t94 = t162 - t217;
t173 = -t94 - t217;
t15 = -t127 * t72 + t128 * t33;
t23 = -t127 * t91 + t128 * t44;
t31 = -t127 * t97 + t128 * t54;
t166 = t133 * t180;
t46 = qJD(1) * t166 + t199;
t56 = t166 + t198;
t172 = qJD(1) * t56 + t46;
t171 = qJD(1) * t81 + t57;
t170 = qJD(1) * t95 + t80;
t169 = t201 * t133;
t168 = pkin(1) * t225;
t167 = 0.2e1 * t174;
t165 = -t85 * t174 + t227;
t161 = t8 * t127 + t7 * t128;
t160 = -t7 * t127 + t8 * t128;
t159 = t132 * t4 + t134 * t5;
t6 = pkin(5) * t196 + pkin(8) * t83 + t12;
t9 = pkin(8) * t82 + t13;
t158 = t132 * t9 - t134 * t6;
t2 = t132 * t6 + t134 * t9;
t156 = -t12 * t128 - t127 * t13;
t155 = t12 * t127 - t128 * t13;
t22 = pkin(5) * t133 + pkin(8) * t208 + t31;
t26 = pkin(8) * t210 + t32;
t154 = -t132 * t26 + t134 * t22;
t153 = t132 * t22 + t134 * t26;
t61 = pkin(2) * t175 - t199;
t73 = pkin(2) * t193 - t198;
t150 = -pkin(7) * t137 - qJD(1) * t73 - t61;
t88 = t221 * t128;
t148 = -qJD(1) * t152 + qJD(5) * t127 + qJD(6) * t88 - t23;
t87 = t221 * t127;
t147 = -qJD(1) * t187 + qJD(5) * t128 - qJD(6) * t87 + t24;
t143 = t85 * t193;
t64 = qJD(2) * pkin(4) + qJD(5) - t79;
t142 = qJD(2) * t163 - t189;
t141 = -t122 * t192 + (qJD(2) * t131 + qJD(5) - t64) * t133;
t17 = (qJD(6) * t83 - t164) * t132 + t220;
t93 = -pkin(7) * t175 + t123;
t96 = t112 + t124;
t139 = t93 * t135 + (t135 * t94 + (-t96 + t112) * t133) * qJD(2);
t18 = qJD(2) * t58 + t224;
t119 = 0.2e1 * t123;
t110 = qJ(4) * t193;
t100 = pkin(5) * t128 + t131;
t99 = -t137 - t212;
t90 = pkin(2) * t196 - t109;
t75 = -t182 * t135 - t118;
t71 = -t110 + t149;
t70 = t136 * t196 + t109;
t68 = t84 * t135;
t67 = t85 * t135;
t65 = t180 + t151;
t50 = t110 + t142;
t38 = -t132 * t83 - t74;
t37 = qJD(1) * t142 + t200;
t36 = -pkin(5) * t82 + t64;
t30 = -t135 * t76 - t143;
t29 = -qJD(2) * t146 + t135 * t77;
t11 = -qJD(2) * t187 + t16;
t10 = t144 + t15;
t1 = [0, 0, 0, t133 * t167, t197 * t225, t204, -t205, 0, -pkin(7) * t204 + t133 * t168, pkin(7) * t205 + t135 * t168, t135 * t150 + t170 * t193, t139, t133 * t150 - t170 * t192, pkin(7) * t139 + t61 * t95 + t73 * t80, t172 * t133 + (t135 * t171 - t71) * qJD(2), -t172 * t135 + (t133 * t171 + t72) * qJD(2), -t133 * t60 + t135 * t51 + (-t133 * t79 - t135 * t65) * qJD(2) + (-t133 * t72 + t135 * t71 + (t133 * t98 - t135 * t97) * qJD(2)) * qJD(1), t46 * t81 + t56 * t57 + t60 * t97 + t65 * t72 + t71 * t79 - t223, t51 * t210 + t71 * t82 + (qJD(1) * t15 + t7) * t133 + (t64 * t211 + t215 + (t135 * t31 + t98 * t211) * qJD(1)) * qJD(2), t51 * t208 + t71 * t83 + (-qJD(1) * t16 - t8) * t133 + (t64 * t209 - t214 + (-t135 * t32 + t98 * t209) * qJD(1)) * qJD(2), t15 * t83 + t16 * t82 + t161 * t135 + ((-t127 * t32 - t128 * t31) * qJD(1) + t156) * t193, t12 * t15 + t13 * t16 + t31 * t7 + t32 * t8 - t64 * t71 - t223, t17 * t68 + t29 * t40, t17 * t67 - t18 * t68 - t29 * t38 + t30 * t40, t29 * t104 + t17 * t133 + (qJD(1) * t68 + t40) * t192, t30 * t104 - t18 * t133 + (qJD(1) * t67 - t38) * t192 (t104 + t196) * t192 (t134 * t10 - t132 * t11) * t104 + t181 * t133 + t50 * t38 + t75 * t18 - t37 * t67 - t36 * t30 + (-t104 * t153 - t133 * t2) * qJD(6) + (qJD(1) * t154 - t158) * t192 -(t132 * t10 + t134 * t11) * t104 - t159 * t133 + t50 * t40 + t75 * t17 + t37 * t68 + t36 * t29 + (-t104 * t154 + t133 * t158) * qJD(6) + (-qJD(1) * t153 - t2) * t192; 0, 0, 0, -t183, t197 * t138, 0, 0, 0, t138 * pkin(1) * t133, pkin(1) * t206 (-t133 * t80 + t135 * t90) * qJD(1) ((t96 - t124) * t133 + (qJD(3) + t173) * t135) * qJD(1), t119 + (t133 * t90 + t135 * t80) * qJD(1), qJ(3) * t93 + qJD(3) * t96 - t80 * t90 + (t133 * t96 + t135 * t173) * qJD(1) * pkin(7), qJD(2) * t89 + t101 + t119 + (t201 * t135 + (-pkin(7) * qJD(2) - t70) * t133) * qJD(1), -qJD(2) * t91 + t103 + ((-qJ(4) * qJD(2) + t70) * t135 + t169) * qJD(1) (-t202 + t65 - t180) * t195, -qJ(3) * t51 + t136 * t60 - t202 * t79 - t57 * t70 - t65 * t91, -t128 * t51 - t202 * t82 + (t127 * t141 - t133 * t23 - t215) * qJD(1), t127 * t51 - t202 * t83 + (t128 * t141 + t133 * t24 + t214) * qJD(1), -t23 * t83 - t24 * t82 + (-qJD(5) * t82 + t12 * t196 - t8) * t128 + (qJD(5) * t83 + t13 * t196 + t7) * t127, t155 * qJD(5) - t12 * t23 + t160 * t122 - t13 * t24 - t131 * t51 + t202 * t64, -t17 * t85 + t219 * t40, t17 * t84 + t18 * t85 - t218 * t40 - t219 * t38, -t40 * t195 + t165, -t176 + (t38 + t213) * t195, -t104 * t195, t100 * t18 - t37 * t84 + t203 * t38 + t218 * t36 + (t132 * t147 + t134 * t148) * t104 + ((t132 * t88 + t134 * t87) * qJD(2) + t158) * t195, t100 * t17 - t37 * t85 + t203 * t40 + t219 * t36 + (-t132 * t148 + t134 * t147) * t104 + (-(t132 * t87 - t134 * t88) * qJD(2) + t2) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, 0, t99, -qJD(2) * t96 + t80 * t196 + t103, t99, t183, 0, t79 * qJD(2) + t103 + (t169 - t177) * qJD(1), -t184 + (t82 - t179) * qJD(2), t185 + (t83 - t178) * qJD(2) (-t127 * t82 - t128 * t83) * t196, -t64 * qJD(2) + t156 * t196 + t160, 0, 0, 0, 0, 0, -qJD(2) * t38 + t165, -t176 + (t195 * t84 - t40) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, 0.2e1 * t175 (-t125 - t126) * t138 (-t135 * t79 + (t65 + t180) * t133) * qJD(1) + t199, -t185 + (-t82 + t191) * t195, -t184 + (-t83 - t194) * t195 (-t127 * t83 + t128 * t82 + (-t127 ^ 2 - t128 ^ 2) * qJD(2)) * t196 (-t133 * t155 + t135 * t64) * qJD(1) + t161, 0, 0, 0, 0, 0, t176 + (t38 - t213) * t195, t227 + (-qJD(2) * t85 + t40) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t83 + t194) * t196 (t82 + t191) * t196, -t82 ^ 2 - t83 ^ 2, -t12 * t83 - t13 * t82 - t51, 0, 0, 0, 0, 0, t18 + t216, t74 * t104 + (-t164 + (qJD(6) + t104) * t83) * t132 + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t38, -t38 ^ 2 + t40 ^ 2, t38 * t104 + t17, -qJD(1) * t143 + t216 - t224, t174, -t226 * t2 - t36 * t40 + t181, t226 * t158 + t36 * t38 - t159;];
tauc_reg  = t1;
