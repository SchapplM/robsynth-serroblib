% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR14V3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:06
% EndTime: 2019-04-12 15:10:19
% DurationCPUTime: 3.66s
% Computational Cost: add. (2787->374), mult. (7352->558), div. (0->0), fcn. (5669->8), ass. (0->186)
t88 = cos(qJ(4));
t171 = qJD(4) * t88;
t89 = cos(qJ(2));
t180 = qJD(2) * t89;
t84 = sin(qJ(4));
t85 = sin(qJ(2));
t100 = t85 * t171 + t84 * t180;
t163 = qJD(2) * qJD(4);
t41 = t100 * qJD(1) + t84 * t163;
t242 = t41 * qJ(3);
t189 = t88 * t89;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t46 = (t87 * t189 + t83 * t85) * qJD(1);
t124 = t87 * t171 - t46;
t169 = qJD(5) * t83;
t241 = -t84 * t169 + t124;
t183 = qJD(1) * t85;
t78 = t88 * qJD(2);
t232 = -t84 * t183 + t78;
t240 = qJD(3) * t232;
t174 = qJD(4) * t84;
t152 = t85 * t174;
t128 = qJD(1) * t152;
t165 = t89 * qJD(1);
t178 = qJD(3) * t84;
t184 = qJ(3) * t88;
t239 = qJ(3) * t128 - ((qJD(4) + t165) * t184 + t178) * qJD(2);
t75 = -qJD(4) + t165;
t115 = t75 * t83;
t156 = t88 * t183;
t182 = qJD(2) * t84;
t57 = t156 + t182;
t39 = t57 * t87 - t115;
t55 = qJD(5) - t232;
t82 = sin(qJ(6));
t86 = cos(qJ(6));
t13 = t39 * t82 - t86 * t55;
t238 = t13 * t55;
t59 = t87 * t75;
t37 = t57 * t83 + t59;
t237 = t37 * t75;
t30 = qJD(6) + t37;
t168 = qJD(5) * t87;
t194 = t86 * t87;
t26 = t194 * t232 + t57 * t82;
t122 = t86 * t168 - t26;
t167 = qJD(6) * t82;
t98 = -t83 * t167 + t122;
t236 = t98 * t30;
t206 = t75 * t89;
t80 = t85 ^ 2;
t119 = qJD(1) * t80 - t206;
t235 = t119 * t84;
t104 = t55 * t169 - t87 * t41;
t209 = t232 * t75;
t153 = t89 * t78;
t40 = qJD(1) * t153 + t88 * t163 - t128;
t234 = t40 - t209;
t208 = t57 * t75;
t233 = t41 - t208;
t181 = qJD(2) * t85;
t140 = qJD(1) * t181;
t11 = t39 * qJD(5) - t87 * t140 + t83 * t40;
t10 = t83 * t140 - t75 * t168 - t57 * t169 + t87 * t40;
t15 = t39 * t86 + t55 * t82;
t4 = t15 * qJD(6) + t82 * t10 - t86 * t41;
t54 = t232 * qJ(3);
t118 = t87 * qJD(3) - t54 * t83;
t23 = t240 - t242;
t7 = t118 * qJD(5) + t87 * t23;
t43 = qJD(3) * t83 + t54 * t87;
t50 = t57 * qJ(3);
t120 = t43 * t82 - t50 * t86;
t24 = qJD(3) * t156 - t239;
t1 = -t120 * qJD(6) + t82 * t24 + t86 * t7;
t230 = 0.2e1 * qJD(1);
t166 = qJD(6) * t86;
t3 = t86 * t10 + t55 * t166 - t39 * t167 + t82 * t41;
t229 = t3 * t82;
t228 = t3 * t83;
t227 = t4 * t83;
t203 = t83 * t23;
t8 = t43 * qJD(5) + t203;
t226 = t8 * t82;
t225 = t8 * t86;
t224 = t87 * t4;
t223 = t10 * t83;
t222 = t11 * t87;
t221 = t13 * t30;
t220 = t15 * t30;
t219 = t15 * t232;
t218 = t24 * t83;
t217 = t24 * t87;
t135 = t30 * t86;
t216 = t37 * t55;
t215 = t37 * t57;
t214 = t39 * t55;
t213 = t39 * t57;
t212 = t40 * t84;
t211 = t41 * t88;
t210 = t55 * t83;
t207 = t57 * t88;
t205 = t82 * t11;
t204 = t82 * t87;
t202 = t83 * t41;
t201 = t83 * t84;
t200 = t83 * t88;
t199 = t84 * t85;
t198 = t84 * t86;
t197 = t84 * t87;
t196 = t85 * t88;
t195 = t86 * t11;
t193 = t86 * t88;
t191 = t87 * t88;
t190 = t87 * t89;
t90 = qJD(2) ^ 2;
t188 = t90 * t85;
t172 = qJD(4) * t87;
t130 = -qJD(6) + t172;
t131 = -qJD(6) * t87 + qJD(4);
t146 = t86 * t169;
t157 = t84 * t165;
t187 = t130 * t193 + (t131 * t82 - t146) * t84 - t82 * t157 - t46 * t86;
t173 = qJD(4) * t86;
t186 = t86 * t157 - t88 * t167 + (t166 * t87 - t173) * t84 + t241 * t82;
t185 = -t89 ^ 2 + t80;
t179 = qJD(3) * t80;
t177 = qJD(3) * t85;
t176 = qJD(3) * t89;
t175 = qJD(4) * t83;
t170 = qJD(5) * t13;
t164 = qJD(4) + t75;
t162 = qJD(3) * qJD(2);
t161 = t84 * t206;
t160 = t75 * t189;
t91 = qJD(1) ^ 2;
t159 = t85 * t91 * t89;
t154 = t85 * t180;
t150 = t75 * t171;
t148 = t82 * t169;
t145 = 0.2e1 * t162;
t144 = qJ(3) * qJD(4) * t75;
t143 = t15 * t169 - t3 * t87;
t142 = t185 * t91;
t137 = -t80 * t91 - t90;
t136 = qJD(2) * t185;
t134 = t55 * t87;
t133 = qJD(3) * t164;
t129 = t88 * t144;
t45 = t165 * t200 - t87 * t183;
t125 = -t83 * t171 + t45;
t25 = t204 * t232 - t86 * t57;
t123 = t82 * t168 - t25;
t17 = t43 * t86 + t50 * t82;
t117 = qJD(5) * t15 - t195;
t116 = t55 * t75;
t114 = t75 * t39;
t112 = -qJD(3) * t55 + qJD(4) * t50;
t111 = t210 * t232 - t104;
t110 = t55 * t201 + t37 * t88;
t109 = t55 * t197 + t39 * t88;
t52 = t85 * t191 - t83 * t89;
t32 = t82 * t199 + t52 * t86;
t53 = t84 * t194 - t82 * t88;
t49 = t82 * t197 + t193;
t108 = t123 * t30;
t107 = -t84 * t176 + t50 * t85;
t106 = -t88 * t176 + t54 * t85;
t105 = -t55 * t168 - t202;
t103 = -t30 * t166 - t205;
t101 = -t131 + t165;
t99 = t84 * t168 - t125;
t97 = t30 * (-qJD(6) - t59);
t96 = t88 * t133 - t84 * t144;
t2 = -t17 * qJD(6) + t86 * t24 - t82 * t7;
t95 = -t15 * t201 + t53 * t30;
t94 = -t13 * t201 + t49 * t30;
t92 = qJ(3) ^ 2;
t51 = t83 * t196 + t190;
t31 = -t85 * t198 + t52 * t82;
t22 = (-qJD(5) + t78) * t190 + (-t84 * t172 + (-qJD(5) * t88 + qJD(2)) * t83) * t85;
t21 = -t83 * t152 - t89 * t169 - t87 * t181 + (t85 * t168 + t83 * t180) * t88;
t6 = (qJD(6) * t199 + t22) * t86 + (-qJD(6) * t52 + t100) * t82;
t5 = t32 * qJD(6) - t100 * t86 + t22 * t82;
t9 = [0, 0, 0, 0.2e1 * t89 * t140, -t136 * t230, t90 * t89, -t188, 0, 0, 0 (-qJ(3) * t136 + t85 * t176) * t230, -qJ(3) * t188 + t89 * t145 (0.4e1 * qJ(3) * t154 + 0.2e1 * t179) * qJD(1) (qJ(3) * t179 + t92 * t154) * t230, t40 * t196 + (-t152 + t153) * t57 (t232 * t88 - t57 * t84) * t180 + (-t212 - t211 + (-t232 * t84 - t207) * qJD(4)) * t85, t75 * t152 - t40 * t89 + (t119 * t88 + t57 * t85) * qJD(2), t85 * t150 + t41 * t89 + (t232 * t85 - t235) * qJD(2) (-t75 - t165) * t181, t24 * t89 + t96 * t85 + (-t119 * t184 - t107) * qJD(2), t23 * t89 + (-t164 * t178 - t129) * t85 + (qJ(3) * t235 - t106) * qJD(2), t10 * t52 + t22 * t39, -t10 * t51 - t11 * t52 - t21 * t39 - t22 * t37, t10 * t199 + t100 * t39 + t22 * t55 + t52 * t41, -t100 * t37 - t11 * t199 - t21 * t55 - t51 * t41, t100 * t55 + t41 * t199, t50 * t21 + t24 * t51 + (qJ(3) * t110 + t118 * t84) * t180 + (t118 * t171 - t8 * t84 + t110 * qJD(3) + ((t175 * t55 + t11) * t88 + (-qJD(4) * t37 - t105) * t84) * qJ(3)) * t85, t50 * t22 + t24 * t52 + (qJ(3) * t109 - t43 * t84) * t180 + (-t43 * t171 - t7 * t84 + t109 * qJD(3) + ((t172 * t55 + t10) * t88 + (-qJD(4) * t39 - t104) * t84) * qJ(3)) * t85, t15 * t6 + t3 * t32, -t13 * t6 - t15 * t5 - t3 * t31 - t32 * t4, t11 * t32 + t15 * t21 + t3 * t51 + t30 * t6, -t11 * t31 - t13 * t21 - t30 * t5 - t4 * t51, t11 * t51 + t21 * t30, -t120 * t21 + t2 * t51 + t8 * t31 - t118 * t5 + t94 * t177 + (t94 * t180 + ((t130 * t30 * t82 - t13 * t175 + t195) * t88 + ((-t148 - t173) * t30 - t227 + (-t103 - t170) * t87) * t84) * t85) * qJ(3), -t1 * t51 - t17 * t21 + t8 * t32 - t118 * t6 + t95 * t177 + (t95 * t180 + ((t130 * t135 - t15 * t175 - t205) * t88 + (-(-qJD(4) * t82 + t146) * t30 - t228 + (-t167 * t30 - t117) * t87) * t84) * t85) * qJ(3); 0, 0, 0, -t159, t142, 0, 0, 0, 0, 0, qJ(3) * t142, 0, -0.2e1 * qJ(3) * t159 + 0.2e1 * t162, qJ(3) * t145 - t92 * t159, -t75 * t207 + t212, -t233 * t84 + t234 * t88, -t150 + (t160 + (-t57 + t182) * t85) * qJD(1), t75 * t174 + (-t161 + (-t232 + t78) * t85) * qJD(1), t75 * t183, t129 + t84 * t133 + ((-t84 * t181 - t160) * qJ(3) + t107) * qJD(1) ((-t78 * t85 + t161) * qJ(3) + t106) * qJD(1) + t96, t10 * t197 + t241 * t39, t46 * t37 + t39 * t45 + (-t37 * t87 - t39 * t83) * t171 + (-t223 - t222 + (t37 * t83 - t39 * t87) * qJD(5)) * t84, -t10 * t88 + t124 * t55 + (-t114 - t104) * t84, t11 * t88 + t125 * t55 + (t105 + t237) * t84, -t84 * t116 - t211, -t45 * t50 + (t8 + t112 * t83 + (t105 - t237) * qJ(3)) * t88 + (t50 * t168 + qJD(3) * t37 + t218 - t75 * t118 + (-t116 * t83 + t11) * qJ(3)) * t84, -t46 * t50 + (t7 + t112 * t87 + (-t114 + t104) * qJ(3)) * t88 + (-t50 * t169 + qJD(3) * t39 + t217 + t75 * t43 + (-t55 * t59 + t10) * qJ(3)) * t84, t15 * t187 + t3 * t53, -t13 * t187 - t15 * t186 - t3 * t49 - t4 * t53, t53 * t11 + t15 * t99 + t187 * t30 + t3 * t201, -t49 * t11 - t13 * t99 - t186 * t30 - t4 * t201, t11 * t201 + t30 * t99, t2 * t201 + t8 * t49 - t186 * t118 - t99 * t120 + ((-t191 * t82 + t198) * t30 + t13 * t200) * qJD(3) + ((t115 * t13 + t82 * t97 + t195) * t84 + (t227 + (t170 - t205) * t87 + (-t101 * t86 + t148) * t30) * t88) * qJ(3), -t1 * t201 + t8 * t53 - t187 * t118 - t99 * t17 + (-(t191 * t86 + t82 * t84) * t30 + t15 * t200) * qJD(3) + ((t115 * t15 + t86 * t97 - t205) * t84 + (t228 + t117 * t87 + (t101 * t82 + t146) * t30) * t88) * qJ(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, 0, t137, t137 * qJ(3), 0, 0, 0, 0, 0, t233, t234, 0, 0, 0, 0, 0, t111 - t215, -t55 ^ 2 * t87 - t202 - t213, 0, 0, 0, 0, 0, -t224 - t108 + (t103 + t238) * t83 (-t195 - t219) * t83 - t236 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t232, -t232 ^ 2 + t57 ^ 2, t40 + t209, -t41 - t208, t140, -t54 * t75 + (-t57 - t156) * qJD(3) + t239, t50 * t75 - 0.2e1 * t240 + t242, t134 * t39 + t223 (t10 - t216) * t87 + (-t11 - t214) * t83, t134 * t55 + t202 - t213, t111 + t215, -t55 * t57, -t118 * t57 - t54 * t37 - t217, -t54 * t39 + t43 * t57 + t218, t15 * t98 + t86 * t228, t26 * t13 + t15 * t25 + (-t13 * t86 - t15 * t82) * t168 + (-t229 - t4 * t86 + (t13 * t82 - t15 * t86) * qJD(6)) * t83 (t195 - t219) * t83 + t236 + t143, t224 - t108 + (t103 - t238) * t83, t30 * t210 - t222, -t2 * t87 - (t50 * t204 + t54 * t86) * t30 - t123 * t118 + (-t118 * t166 - t120 * t55 + t50 * t13 + t226) * t83, t1 * t87 + (-t194 * t50 + t54 * t82) * t30 - t122 * t118 + (t118 * t167 + t50 * t15 - t17 * t55 + t225) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t10 + t216, -t11 + t214, t41, -t50 * t39 - t203 + (-qJD(5) + t55) * t43, t118 * t55 + t50 * t37 - t7, t135 * t15 + t229 (t3 - t221) * t86 + (-t4 - t220) * t82, t135 * t30 - t15 * t39 + t205, -t30 ^ 2 * t82 + t13 * t39 + t195, -t30 * t39, t120 * t39 - t43 * t13 - t225, -t43 * t15 + t17 * t39 + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t13, -t13 ^ 2 + t15 ^ 2, t3 + t221, t220 - t4, t11, t118 * t15 + t17 * t30 + t2, -t118 * t13 - t120 * t30 - t1;];
tauc_reg  = t9;
