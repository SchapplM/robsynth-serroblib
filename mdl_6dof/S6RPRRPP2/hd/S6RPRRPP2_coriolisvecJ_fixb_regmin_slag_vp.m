% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:33:21
% EndTime: 2019-03-09 04:33:29
% DurationCPUTime: 2.77s
% Computational Cost: add. (2966->378), mult. (6688->468), div. (0->0), fcn. (3950->6), ass. (0->193)
t129 = cos(qJ(3));
t199 = qJD(1) * t129;
t110 = -qJD(4) + t199;
t128 = cos(qJ(4));
t126 = sin(qJ(4));
t198 = qJD(3) * t126;
t127 = sin(qJ(3));
t200 = qJD(1) * t127;
t91 = t128 * t200 + t198;
t225 = t110 * t91;
t194 = qJD(4) * t128;
t176 = t127 * t194;
t196 = qJD(3) * t129;
t180 = t126 * t196;
t188 = qJD(3) * qJD(4);
t49 = (t176 + t180) * qJD(1) + t126 * t188;
t259 = -t49 + t225;
t189 = qJD(1) * qJD(3);
t173 = t127 * t189;
t245 = pkin(4) + pkin(5);
t258 = t245 * t173;
t257 = qJ(5) * t173 - t110 * qJD(5);
t106 = t110 * qJ(5);
t118 = t127 * qJD(2);
t112 = sin(pkin(9)) * pkin(1) + pkin(7);
t98 = t112 * qJD(1);
t65 = t129 * t98 + t118;
t57 = qJD(3) * pkin(8) + t65;
t113 = -cos(pkin(9)) * pkin(1) - pkin(2);
t79 = -pkin(3) * t129 - pkin(8) * t127 + t113;
t60 = t79 * qJD(1);
t24 = t126 * t60 + t128 * t57;
t18 = -t106 + t24;
t165 = pkin(4) * t173;
t195 = qJD(4) * t126;
t249 = t129 * qJD(2) - t127 * t98;
t58 = t249 * qJD(3);
t160 = pkin(3) * t127 - pkin(8) * t129;
t96 = t160 * qJD(3);
t78 = qJD(1) * t96;
t168 = t126 * t58 - t128 * t78 + t57 * t194 + t60 * t195;
t5 = -t165 + t168;
t256 = t110 * t18 + t5;
t190 = t128 * qJD(3);
t89 = t126 * t200 - t190;
t255 = t49 * qJ(6) + t89 * qJD(6);
t161 = qJD(3) * pkin(3) + t249;
t141 = qJ(5) * t91 + t161;
t14 = -t245 * t89 + qJD(6) + t141;
t254 = (qJD(6) + t14) * t91;
t86 = t91 ^ 2;
t253 = -t110 ^ 2 - t86;
t252 = -qJD(5) * t126 - t65;
t226 = t110 * t89;
t175 = t128 * t189;
t177 = t127 * t195;
t48 = qJD(1) * t177 - t128 * t188 - t129 * t175;
t251 = -t48 + t226;
t219 = t129 * t49;
t84 = t110 * t176;
t250 = t84 - t219;
t23 = -t126 * t57 + t128 * t60;
t204 = qJD(5) - t23;
t248 = 0.2e1 * t257;
t216 = qJ(5) * t126;
t247 = -t245 * t128 - t216;
t246 = t89 ^ 2;
t59 = qJD(3) * t118 + t98 * t196;
t8 = t49 * pkin(4) + t48 * qJ(5) - t91 * qJD(5) + t59;
t244 = t126 * t8;
t243 = t128 * t8;
t22 = pkin(4) * t89 - t141;
t242 = t22 * t91;
t241 = t91 * t89;
t240 = pkin(8) - qJ(6);
t102 = t240 * t128;
t95 = t160 * qJD(1);
t171 = -t126 * t249 + t128 * t95;
t208 = t128 * t129;
t239 = (-qJ(6) * t208 - t245 * t127) * qJD(1) - t171 - qJD(4) * t102 + qJD(6) * t126;
t191 = qJD(6) * t128;
t234 = t126 * t95 + t128 * t249;
t28 = qJ(5) * t200 + t234;
t238 = qJ(6) * t126 * t199 + t240 * t195 + t191 + t28;
t215 = qJ(5) * t128;
t145 = -t245 * t126 + t215;
t237 = t110 * t145 + t252;
t156 = pkin(4) * t126 - t215;
t236 = t110 * t156 - t252;
t179 = t129 * t190;
t209 = t127 * t128;
t235 = -t89 * t179 - t49 * t209;
t94 = t112 * t208;
t233 = qJD(4) * t94 + t79 * t195;
t232 = t126 * t96 + t79 * t194;
t231 = t126 * t79 + t94;
t230 = qJ(5) * t49;
t229 = qJ(5) * t89;
t16 = qJ(6) * t89 + t24;
t10 = -t106 + t16;
t228 = t10 * t110;
t224 = t126 * t48;
t223 = t126 * t161;
t222 = t127 * t89;
t221 = t128 * t161;
t220 = t129 * t48;
t218 = t59 * t126;
t217 = t59 * t128;
t214 = qJ(6) * t127;
t213 = qJD(4) * t89;
t212 = t110 * t128;
t211 = t112 * t126;
t210 = t126 * t129;
t131 = qJD(3) ^ 2;
t207 = t131 * t127;
t206 = t131 * t129;
t15 = qJ(6) * t91 + t23;
t205 = qJD(5) - t15;
t122 = t127 ^ 2;
t202 = -t129 ^ 2 + t122;
t99 = qJD(1) * t113;
t201 = qJD(1) * t122;
t197 = qJD(3) * t127;
t192 = qJD(5) * t128;
t187 = t126 * t78 + t128 * t58 + t60 * t194;
t186 = pkin(8) * t110 * t126;
t185 = pkin(8) * t212;
t184 = pkin(8) * t197;
t183 = pkin(8) * t190;
t182 = t91 * t196;
t181 = t126 * t201;
t178 = t110 * t195;
t172 = -pkin(4) - t211;
t93 = t112 * t210;
t170 = t128 * t79 - t93;
t69 = t91 * t197;
t169 = t69 + t220;
t167 = -t48 + t213;
t104 = t122 * t175;
t166 = t104 - t220;
t164 = t91 * t176;
t163 = t128 * t96 - t233;
t162 = t110 * t180 + t250;
t33 = -qJ(5) * t129 + t231;
t9 = t245 * t110 + t205;
t159 = t10 * t128 + t126 * t9;
t158 = t10 * t126 - t128 * t9;
t157 = pkin(4) * t128 + t216;
t17 = pkin(4) * t110 + t204;
t155 = -t126 * t18 + t128 * t17;
t154 = t126 * t17 + t128 * t18;
t152 = 0.2e1 * qJD(3) * t99;
t151 = qJ(5) * t197 - qJD(5) * t129 + t232;
t3 = -pkin(5) * t49 - t8;
t150 = -t126 * t3 - t14 * t194;
t149 = t128 * t3 - t14 * t195;
t148 = t112 + t156;
t147 = qJ(6) * t48 + t168;
t146 = -t110 * t24 - t168;
t144 = -t173 + t241;
t143 = t57 * t195 - t187;
t83 = t110 * t179;
t142 = t110 * t177 + t104 - t83;
t140 = (t110 * t129 - t201) * t126;
t139 = -t112 + t145;
t4 = -t143 + t257;
t137 = -t110 * t23 + t143;
t136 = t48 + t226;
t135 = t155 * qJD(4) + t126 * t5 + t128 * t4;
t134 = t147 - t258;
t132 = qJD(1) ^ 2;
t121 = t129 * pkin(4);
t101 = t240 * t126;
t97 = -pkin(3) - t157;
t80 = pkin(3) - t247;
t70 = t89 * t197;
t50 = t148 * t127;
t37 = pkin(4) * t91 + t229;
t35 = t139 * t127;
t34 = t121 - t170;
t31 = t126 * t214 + t33;
t30 = -pkin(4) * t200 - t171;
t27 = -t245 * t91 - t229;
t25 = pkin(5) * t129 + t121 + t93 + (-t79 - t214) * t128;
t20 = (t157 * qJD(4) - t192) * t127 + t148 * t196;
t13 = t172 * t197 - t163;
t12 = (t247 * qJD(4) + t192) * t127 + t139 * t196;
t11 = (-t127 * t190 - t129 * t195) * t112 + t151;
t7 = (qJ(6) * qJD(4) - qJD(3) * t112) * t209 + (qJD(6) * t127 + (qJ(6) * qJD(3) - qJD(4) * t112) * t129) * t126 + t151;
t6 = (-qJ(6) * t196 - t96) * t128 + (qJ(6) * t195 - t191 + (-pkin(5) + t172) * qJD(3)) * t127 + t233;
t2 = -qJD(6) * t91 + t134;
t1 = t4 + t255;
t19 = [0, 0, 0, 0, 0.2e1 * t129 * t173, -0.2e1 * t202 * t189, t206, -t207, 0, -t112 * t206 + t127 * t152, t112 * t207 + t129 * t152, t91 * t179 + (-t128 * t48 - t91 * t195) * t127, -t164 + (-t182 + (t48 + t213) * t127) * t126 + t235, t142 + t169, t84 + t219 + (t140 - t222) * qJD(3) (-t110 - t199) * t197, -t163 * t110 + ((t112 * t89 - t223) * qJD(3) + t168) * t129 + (-t161 * t194 + t112 * t49 + t218 + (qJD(1) * t170 - t110 * t211 + t23) * qJD(3)) * t127, t232 * t110 + ((-t110 * t112 - t57) * t195 + (t112 * t91 - t221) * qJD(3) + t187) * t129 + (t161 * t195 - t112 * t48 + t217 + (-t231 * qJD(1) - t112 * t212 - t24) * qJD(3)) * t127, t110 * t13 + t20 * t89 + t49 * t50 + (t198 * t22 + t5) * t129 + (t22 * t194 + t244 + (-qJD(1) * t34 - t17) * qJD(3)) * t127, -t11 * t89 + t13 * t91 - t33 * t49 - t34 * t48 + t155 * t196 + (-qJD(4) * t154 - t126 * t4 + t128 * t5) * t127, -t11 * t110 - t20 * t91 + t48 * t50 + (-t190 * t22 - t4) * t129 + (t22 * t195 - t243 + (qJD(1) * t33 + t18) * qJD(3)) * t127, t11 * t18 + t13 * t17 + t20 * t22 + t33 * t4 + t34 * t5 + t50 * t8, t110 * t6 - t12 * t89 - t35 * t49 + (-t14 * t198 + t2) * t129 + ((-qJD(1) * t25 - t9) * qJD(3) + t150) * t127, -t110 * t7 + t12 * t91 - t35 * t48 + (t14 * t190 - t1) * t129 + ((qJD(1) * t31 + t10) * qJD(3) + t149) * t127, t25 * t48 + t31 * t49 - t6 * t91 + t7 * t89 + t158 * t196 + (qJD(4) * t159 + t1 * t126 - t128 * t2) * t127, t1 * t31 + t10 * t7 + t12 * t14 + t2 * t25 + t3 * t35 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -t206, 0, 0, 0, 0, 0, qJD(3) * t140 + t250 + t70, t69 + (-t177 + t179) * t110 - t166, -qJD(3) * t181 + t162 + t70, t164 + (t127 * t167 + t182) * t126 + t235, -t83 + (-qJD(3) * t91 + t178) * t127 + t166 (qJD(3) * t154 - t8) * t129 + (qJD(3) * t22 + t135) * t127 (-t181 + t222) * qJD(3) + t162, t142 - t169 (-t126 * t91 + t128 * t89) * t196 + (t224 + t128 * t49 + (-t126 * t89 - t128 * t91) * qJD(4)) * t127 (qJD(3) * t159 + t3) * t129 + (-qJD(3) * t14 - qJD(4) * t158 + t1 * t128 + t126 * t2) * t127; 0, 0, 0, 0, -t127 * t132 * t129, t202 * t132, 0, 0, 0, qJD(3) * t65 - t99 * t200 - t59, -t99 * t199, -t91 * t212 - t224, t126 * t259 + t251 * t128, -t110 * t194 + (t110 * t208 + (-t91 + t198) * t127) * qJD(1), t178 + (-t110 * t210 + (t89 + t190) * t127) * qJD(1), t110 * t200, -pkin(3) * t49 - t217 + t171 * t110 - t65 * t89 + (t185 - t223) * qJD(4) + (-t23 * t127 + (t129 * t161 - t184) * t126) * qJD(1), pkin(3) * t48 + t218 - t234 * t110 - t65 * t91 + (-t186 - t221) * qJD(4) + (t161 * t208 + (t24 - t183) * t127) * qJD(1), -t110 * t30 - t243 + t49 * t97 - t236 * t89 + (t126 * t22 + t185) * qJD(4) + (t127 * t17 + (-t129 * t22 - t184) * t126) * qJD(1), t28 * t89 - t30 * t91 + (t4 - t110 * t17 + (qJD(4) * t91 - t49) * pkin(8)) * t128 + (pkin(8) * t167 + t256) * t126, t110 * t28 - t244 + t48 * t97 + t236 * t91 + (-t128 * t22 + t186) * qJD(4) + (t22 * t208 + (-t18 + t183) * t127) * qJD(1), t135 * pkin(8) - t17 * t30 - t18 * t28 - t236 * t22 + t8 * t97, -t49 * t80 + t237 * t89 - t239 * t110 + (t14 * t210 + (-qJD(3) * t101 + t9) * t127) * qJD(1) + t149, -t48 * t80 - t237 * t91 + t238 * t110 + (-t14 * t208 + (qJD(3) * t102 - t10) * t127) * qJD(1) - t150, t101 * t48 + t102 * t49 + t239 * t91 - t238 * t89 + (t110 * t9 - t1) * t128 + (-t2 - t228) * t126, t1 * t102 - t238 * t10 + t101 * t2 - t237 * t14 - t239 * t9 + t3 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t86 - t246, -t136, -t49 - t225, t173, t161 * t91 + t146, -t161 * t89 + t137, -t37 * t89 + t146 + 0.2e1 * t165 - t242, pkin(4) * t48 - t230 + (t18 - t24) * t91 + (t17 - t204) * t89, -t22 * t89 + t37 * t91 - t137 + t248, -pkin(4) * t5 + qJ(5) * t4 - t17 * t24 + t18 * t204 - t22 * t37, -t110 * t16 + t27 * t89 - t147 + t254 + 0.2e1 * t258, t110 * t15 + t14 * t89 - t27 * t91 - t143 + t248 + t255, t230 - t245 * t48 + (-t10 + t16) * t91 + (-t9 + t205) * t89, qJ(5) * t1 + t10 * t205 - t14 * t27 - t16 * t9 - t2 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t136, t253, t242 + t256, t144, t253, t136, t134 + t228 - t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t251, -t86 - t246, -t10 * t89 + t9 * t91 + t3;];
tauc_reg  = t19;
