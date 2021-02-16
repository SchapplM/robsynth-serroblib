% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:38:48
% EndTime: 2021-01-16 02:39:01
% DurationCPUTime: 3.17s
% Computational Cost: add. (4232->353), mult. (10952->491), div. (0->0), fcn. (8353->10), ass. (0->188)
t136 = sin(qJ(2));
t132 = sin(pkin(6));
t196 = qJD(1) * t132;
t179 = t136 * t196;
t135 = sin(qJ(3));
t190 = t135 * qJD(3);
t244 = pkin(3) * t190 - t179;
t131 = sin(pkin(11));
t138 = cos(qJ(3));
t226 = qJ(4) + pkin(8);
t173 = qJD(3) * t226;
t148 = -t135 * qJD(4) - t138 * t173;
t210 = cos(pkin(11));
t99 = t138 * qJD(4) - t135 * t173;
t58 = t131 * t148 + t210 * t99;
t172 = t210 * t138;
t150 = -t131 * t135 + t172;
t139 = cos(qJ(2));
t178 = t139 * t196;
t82 = t150 * t178;
t219 = t58 - t82;
t111 = t131 * t138 + t135 * t210;
t103 = t111 * qJD(3);
t106 = t150 * qJD(3);
t243 = t103 * pkin(4) - t106 * pkin(9) + t244;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t113 = qJD(2) * pkin(8) + t179;
t168 = qJ(4) * qJD(2) + t113;
t133 = cos(pkin(6));
t195 = qJD(1) * t133;
t177 = t135 * t195;
t80 = t138 * t168 + t177;
t176 = t210 * t80;
t121 = t138 * t195;
t79 = -t135 * t168 + t121;
t75 = qJD(3) * pkin(3) + t79;
t32 = t131 * t75 + t176;
t29 = qJD(3) * pkin(9) + t32;
t120 = qJD(2) * t172;
t194 = qJD(2) * t135;
t101 = t131 * t194 - t120;
t104 = t111 * qJD(2);
t127 = -t138 * pkin(3) - pkin(2);
t97 = qJD(2) * t127 + qJD(4) - t178;
t43 = t101 * pkin(4) - t104 * pkin(9) + t97;
t16 = t134 * t43 + t137 * t29;
t161 = qJD(4) + t178;
t47 = (-t135 * t113 + t121) * qJD(3) + (-qJ(4) * t190 + t138 * t161) * qJD(2);
t48 = (-t138 * t113 - t177) * qJD(3) + (-qJD(3) * t138 * qJ(4) - t135 * t161) * qJD(2);
t21 = t131 * t48 + t210 * t47;
t188 = qJD(2) * qJD(3);
t175 = t135 * t188;
t100 = pkin(3) * t175 + qJD(2) * t179;
t118 = t131 * t175;
t149 = qJD(3) * t120 - t118;
t96 = qJD(2) * t103;
t40 = t96 * pkin(4) - pkin(9) * t149 + t100;
t38 = t137 * t40;
t144 = -qJD(5) * t16 - t134 * t21 + t38;
t189 = t137 * qJD(3);
t192 = qJD(5) * t134;
t157 = qJD(5) * t189 - t104 * t192 + t137 * t149;
t143 = -qJ(6) * t157 + t144;
t235 = t96 * pkin(5);
t88 = t134 * qJD(3) + t137 * t104;
t1 = -t88 * qJD(6) + t143 + t235;
t86 = t134 * t104 - t189;
t10 = -t86 * qJ(6) + t16;
t98 = qJD(5) + t101;
t231 = t10 * t98;
t242 = t1 + t231;
t241 = t134 * t82 + t243 * t137;
t220 = -t111 * t178 + t131 * t99 - t210 * t148;
t240 = t137 * t96 - t98 * t192;
t191 = qJD(5) * t137;
t70 = -pkin(4) * t150 - t111 * pkin(9) + t127;
t239 = t243 * t134 + t219 * t137 + t70 * t191;
t50 = qJD(5) * t88 + t134 * t149;
t238 = t88 ^ 2;
t15 = -t134 * t29 + t137 * t43;
t9 = -t88 * qJ(6) + t15;
t7 = t98 * pkin(5) + t9;
t237 = t7 - t9;
t236 = t86 * pkin(5);
t159 = -qJ(6) * t106 - qJD(6) * t111;
t209 = qJ(6) * t111;
t116 = t226 * t135;
t117 = t226 * t138;
t84 = -t131 * t116 + t117 * t210;
t76 = t137 * t84;
t234 = t103 * pkin(5) - t134 * t58 + t159 * t137 + (-t76 + (-t70 + t209) * t134) * qJD(5) + t241;
t180 = t111 * t191;
t233 = -qJ(6) * t180 + (-qJD(5) * t84 + t159) * t134 + t239;
t232 = pkin(3) * t135;
t20 = t131 * t47 - t210 * t48;
t8 = t50 * pkin(5) + t20;
t230 = t8 * t134;
t229 = t8 * t137;
t228 = t86 * t98;
t227 = t88 * t98;
t124 = t131 * pkin(3) + pkin(9);
t198 = qJ(6) + t124;
t170 = qJD(5) * t198;
t208 = t101 * t137;
t71 = t131 * t80;
t36 = t210 * t79 - t71;
t185 = pkin(3) * t194;
t59 = t104 * pkin(4) + t101 * pkin(9) + t185;
t54 = t137 * t59;
t225 = t104 * pkin(5) + qJ(6) * t208 + t137 * t170 + t54 + (qJD(6) - t36) * t134;
t202 = t134 * t101;
t222 = t134 * t59 + t137 * t36;
t224 = qJ(6) * t202 - t137 * qJD(6) + t134 * t170 + t222;
t201 = t134 * t106;
t223 = (t180 + t201) * pkin(5) + t220;
t221 = -t134 * t50 - t86 * t191;
t218 = t134 * t70 + t76;
t217 = qJD(2) * pkin(2);
t216 = t104 * t86;
t215 = t134 * t88;
t214 = t134 * t96;
t213 = t137 * t86;
t212 = t157 * t134;
t211 = t88 * t104;
t207 = t106 * t137;
t206 = t111 * t134;
t205 = t132 * t136;
t204 = t132 * t139;
t141 = qJD(2) ^ 2;
t203 = t132 * t141;
t140 = qJD(3) ^ 2;
t200 = t140 * t135;
t199 = t140 * t138;
t197 = t135 ^ 2 - t138 ^ 2;
t193 = qJD(2) * t136;
t183 = t136 * t203;
t182 = t132 * t193;
t181 = qJD(2) * t204;
t174 = qJD(6) + t236;
t34 = t131 * t79 + t176;
t171 = -t134 * t40 - t137 * t21 - t43 * t191 + t29 * t192;
t83 = t210 * t116 + t131 * t117;
t169 = t137 * t98;
t167 = t135 * t181;
t166 = t138 * t181;
t125 = -pkin(3) * t210 - pkin(4);
t165 = qJD(5) * t124 * t98 + t20;
t156 = t50 * qJ(6) + t171;
t4 = -t86 * qJD(6) - t156;
t164 = -t7 * t98 + t4;
t31 = t210 * t75 - t71;
t162 = t20 * t111 - t84 * t96;
t160 = -t213 - t215;
t158 = -t98 * t202 + t240;
t155 = t137 * t157 - t192 * t88;
t107 = t133 * t135 + t138 * t205;
t153 = t133 * t138 - t135 * t205;
t65 = t107 * t210 + t131 * t153;
t45 = -t134 * t65 - t137 * t204;
t154 = t134 * t204 - t137 * t65;
t152 = t217 * qJD(2);
t28 = -qJD(3) * pkin(4) - t31;
t151 = -t124 * t96 + t28 * t98;
t147 = t107 * qJD(3);
t146 = -0.2e1 * qJD(3) * t217;
t142 = -t147 - t167;
t115 = -t137 * pkin(5) + t125;
t109 = t198 * t137;
t108 = t198 * t134;
t85 = t86 ^ 2;
t78 = qJD(3) * t153 + t166;
t67 = t137 * t70;
t64 = t131 * t107 - t153 * t210;
t56 = pkin(5) * t206 + t83;
t35 = t131 * t142 + t210 * t78;
t33 = t131 * t78 - t142 * t210;
t25 = -pkin(5) * t202 + t34;
t24 = -qJ(6) * t206 + t218;
t23 = t174 + t28;
t22 = -pkin(5) * t150 - t134 * t84 - t137 * t209 + t67;
t18 = -t137 * t98 ^ 2 - t211 - t214;
t17 = t158 - t216;
t14 = qJD(5) * t154 - t134 * t35 + t137 * t182;
t13 = qJD(5) * t45 + t134 * t182 + t137 * t35;
t3 = t14 * t98 + t33 * t86 + t45 * t96 + t64 * t50;
t2 = -t13 * t98 + t154 * t96 + t157 * t64 + t33 * t88;
t5 = [0, 0, -t183, -t139 * t203, 0, 0, 0, 0, 0, -t138 * t183 + (-t147 - 0.2e1 * t167) * qJD(3), t135 * t183 + (-t78 - t166) * qJD(3), -t33 * qJD(3) + (t101 * t193 - t139 * t96) * t132, -t35 * qJD(3) + (t104 * t193 - t139 * t149) * t132, -t35 * t101 + t33 * t104 + t149 * t64 - t65 * t96, t20 * t64 + t21 * t65 - t31 * t33 + t32 * t35 + (-t100 * t139 + t193 * t97) * t132, 0, 0, 0, 0, 0, t3, t2, t3, t2, -t13 * t86 - t14 * t88 + t154 * t50 - t157 * t45, t1 * t45 + t10 * t13 + t7 * t14 - t154 * t4 + t23 * t33 + t8 * t64; 0, 0, 0, 0, 0.2e1 * t138 * t175, -0.2e1 * t197 * t188, t199, -t200, 0, -pkin(8) * t199 + t135 * t146, pkin(8) * t200 + t138 * t146, -t101 * t179 - t100 * t150 + t97 * t103 + t127 * t96 + (t101 * t232 - t220) * qJD(3), -t104 * t179 + t100 * t111 + t97 * t106 - t127 * t118 + (t104 * t232 + t120 * t127 - t219) * qJD(3), -t101 * t219 - t32 * t103 + t104 * t220 - t31 * t106 + t149 * t83 + t150 * t21 + t162, t100 * t127 + t20 * t83 + t21 * t84 + t219 * t32 - t220 * t31 + t244 * t97, t111 * t155 + t207 * t88, t160 * t106 + (-t212 - t137 * t50 + (t134 * t86 - t137 * t88) * qJD(5)) * t111, t88 * t103 + t111 * t240 - t150 * t157 + t98 * t207, -t98 * t201 - t86 * t103 + t50 * t150 + (-t191 * t98 - t214) * t111, t98 * t103 - t150 * t96, t67 * t96 - (-t191 * t29 + t38) * t150 + t15 * t103 + t83 * t50 + t28 * t180 + (-t191 * t84 + t241) * t98 + t220 * t86 + ((-qJD(5) * t70 - t58) * t98 - (-qJD(5) * t43 - t21) * t150 + t28 * t106 + t162) * t134, -t218 * t96 - t171 * t150 - t16 * t103 + t83 * t157 + t28 * t207 + (t192 * t84 - t239) * t98 + t220 * t88 + (t20 * t137 - t192 * t28) * t111, t23 * t201 - t1 * t150 + t7 * t103 + t22 * t96 + t56 * t50 + t234 * t98 + t223 * t86 + (t191 * t23 + t230) * t111, t23 * t207 - t10 * t103 + t4 * t150 - t24 * t96 + t56 * t157 - t233 * t98 + t223 * t88 + (-t192 * t23 + t229) * t111, -t22 * t157 - t24 * t50 - t234 * t88 - t233 * t86 + (-t10 * t134 - t137 * t7) * t106 + (-t1 * t137 - t4 * t134 + (-t10 * t137 + t134 * t7) * qJD(5)) * t111, t1 * t22 + t233 * t10 + t223 * t23 + t234 * t7 + t4 * t24 + t8 * t56; 0, 0, 0, 0, -t135 * t141 * t138, t197 * t141, 0, 0, 0, t135 * t152, t138 * t152, t34 * qJD(3) - t101 * t185 - t97 * t104 - t20, t36 * qJD(3) + t97 * t101 - t104 * t185 - t21, (t32 - t34) * t104 + (t36 - t31) * t101 + (-t131 * t96 - t149 * t210) * pkin(3), t31 * t34 - t32 * t36 + (t131 * t21 - t194 * t97 - t20 * t210) * pkin(3), t169 * t88 + t212, t101 * t160 + t155 + t221, t169 * t98 - t211 + t214, t158 + t216, -t98 * t104, -t15 * t104 + t125 * t50 - t34 * t86 - t54 * t98 - t165 * t137 + (t36 * t98 + t151) * t134, t16 * t104 + t125 * t157 + t134 * t165 + t137 * t151 + t222 * t98 - t34 * t88, -t7 * t104 - t108 * t96 + t115 * t50 - t229 - t25 * t86 - t225 * t98 + (t101 * t23 + (t23 + t236) * qJD(5)) * t134, t23 * t208 + t10 * t104 - t109 * t96 + t115 * t157 + t230 - t25 * t88 + t224 * t98 + (pkin(5) * t215 + t137 * t23) * qJD(5), t108 * t157 - t109 * t50 - t242 * t134 + t164 * t137 + t224 * t86 + t225 * t88, -t1 * t108 + t4 * t109 + t8 * t115 - t225 * t7 + (pkin(5) * t192 - t25) * t23 - t224 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t104 * qJD(3), -t118 + (t120 - t101) * qJD(3), -t101 ^ 2 - t104 ^ 2, t32 * t101 + t31 * t104 + t100, 0, 0, 0, 0, 0, t17, t18, t17, t18, (-t213 + t215) * t101 - t155 + t221, -t23 * t104 + t164 * t134 + t242 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t86, -t85 + t238, t157 + t228, t227 - t50, t96, t16 * t98 - t28 * t88 + t144, t15 * t98 + t28 * t86 + t171, 0.2e1 * t235 + t231 + (-t174 - t23) * t88 + t143, -t238 * pkin(5) + t9 * t98 + (qJD(6) + t23) * t86 + t156, -pkin(5) * t157 - t237 * t86, t237 * t10 + (-t23 * t88 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 + t227, t157 - t228, -t85 - t238, t10 * t86 + t7 * t88 + t8;];
tauc_reg = t5;
