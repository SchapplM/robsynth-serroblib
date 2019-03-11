% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:06
% EndTime: 2019-03-09 05:57:15
% DurationCPUTime: 2.68s
% Computational Cost: add. (5621->336), mult. (12784->436), div. (0->0), fcn. (8731->8), ass. (0->187)
t143 = sin(qJ(4));
t144 = sin(qJ(3));
t146 = cos(qJ(3));
t247 = cos(qJ(4));
t121 = t143 * t146 + t247 * t144;
t202 = qJD(3) + qJD(4);
t91 = t202 * t121;
t262 = t91 * qJD(1);
t145 = cos(qJ(5));
t205 = qJD(5) * t145;
t193 = qJD(1) * t247;
t209 = qJD(1) * t144;
t113 = t143 * t209 - t146 * t193;
t216 = t113 * t145;
t261 = t205 + t216;
t142 = sin(qJ(5));
t206 = qJD(5) * t142;
t217 = t113 * t142;
t260 = t206 + t217;
t131 = sin(pkin(10)) * pkin(1) + pkin(7);
t243 = pkin(8) + t131;
t112 = qJD(5) + t113;
t188 = t243 * qJD(1);
t104 = t144 * qJD(2) + t188 * t146;
t192 = qJD(4) * t247;
t207 = qJD(4) * t143;
t103 = t146 * qJD(2) - t188 * t144;
t94 = t103 * qJD(3);
t95 = t104 * qJD(3);
t236 = qJD(3) * pkin(3);
t98 = t103 + t236;
t21 = -t104 * t207 - t143 * t95 + t98 * t192 + t247 * t94;
t164 = -t143 * t144 + t247 * t146;
t160 = t164 * qJD(4);
t40 = t262 * pkin(4) + (-pkin(9) * t160 + (t144 * pkin(3) - t164 * pkin(9)) * qJD(3)) * qJD(1);
t97 = t247 * t104;
t59 = t143 * t98 + t97;
t54 = t202 * pkin(9) + t59;
t208 = qJD(1) * t146;
t115 = -t143 * t208 - t144 * t193;
t132 = -cos(pkin(10)) * pkin(1) - pkin(2);
t122 = -pkin(3) * t146 + t132;
t116 = t122 * qJD(1);
t71 = t113 * pkin(4) + t115 * pkin(9) + t116;
t166 = t142 * t40 + t145 * t21 + t71 * t205 - t54 * t206;
t223 = t262 * qJ(6);
t2 = qJD(6) * t112 + t166 + t223;
t189 = t142 * t21 - t145 * t40 + t54 * t205 + t71 * t206;
t248 = pkin(5) * t262;
t4 = t189 - t248;
t259 = t4 * t142 + t2 * t145;
t228 = t145 * t262;
t258 = (t112 * t206 - t228) * pkin(9);
t90 = t164 * qJD(3) + t160;
t150 = t90 * qJD(1);
t162 = t145 * t115 - t142 * t202;
t49 = -t162 * qJD(5) + t142 * t150;
t186 = t145 * t202;
t99 = -t115 * t142 - t186;
t257 = t164 * t49 - t91 * t99;
t256 = t90 * t202;
t27 = t142 * t71 + t145 * t54;
t16 = qJ(6) * t112 + t27;
t96 = t143 * t104;
t58 = t247 * t98 - t96;
t53 = -t202 * pkin(4) - t58;
t28 = t99 * pkin(5) + qJ(6) * t162 + t53;
t255 = t16 * t115 - t28 * t216;
t22 = t104 * t192 + t143 * t94 + t98 * t207 + t247 * t95;
t48 = -qJD(5) * t186 - t115 * t206 - t145 * t150;
t6 = pkin(5) * t49 + qJ(6) * t48 + qJD(6) * t162 + t22;
t254 = -t6 * t145 + t28 * t206;
t253 = t22 * t145 - t53 * t206;
t65 = t143 * t103 + t97;
t183 = pkin(3) * t207 - t65;
t252 = pkin(5) * t260 - qJ(6) * t261 - qJD(6) * t142;
t118 = t243 * t144;
t119 = t243 * t146;
t251 = -t247 * t118 - t143 * t119;
t250 = t162 ^ 2;
t249 = t112 ^ 2;
t246 = pkin(5) * t115;
t245 = t6 * t142;
t214 = t121 * t145;
t227 = t145 * t90;
t242 = -t49 * t214 - t99 * t227;
t84 = -pkin(4) * t115 + pkin(9) * t113;
t241 = t142 * t84 + t145 * t58;
t66 = t247 * t103 - t96;
t76 = pkin(3) * t209 + t84;
t240 = t142 * t76 + t145 * t66;
t239 = -t162 * t91 + t164 * t48;
t79 = -pkin(4) * t164 - pkin(9) * t121 + t122;
t82 = -t143 * t118 + t247 * t119;
t238 = t142 * t79 + t145 * t82;
t237 = pkin(3) * qJD(4);
t235 = t162 * t28;
t234 = t162 * t99;
t233 = t112 * t99;
t134 = pkin(3) * t143 + pkin(9);
t232 = t134 * t262;
t231 = t142 * t262;
t230 = t142 * t90;
t229 = t142 * t99;
t224 = t48 * t142;
t222 = t252 + t183;
t221 = -t59 + t252;
t220 = t162 * t112;
t219 = t162 * t145;
t218 = t112 * t115;
t187 = t112 * t145;
t215 = t115 * t113;
t147 = qJD(3) ^ 2;
t213 = t147 * t144;
t212 = t147 * t146;
t26 = -t142 * t54 + t145 * t71;
t211 = qJD(6) - t26;
t210 = t144 ^ 2 - t146 ^ 2;
t125 = qJD(1) * t132;
t203 = qJD(1) * qJD(3);
t201 = t247 * pkin(3);
t199 = t144 * t236;
t198 = t162 * t230;
t69 = t262 * t214;
t197 = t142 * t247;
t196 = t145 * t247;
t194 = t121 * t206;
t191 = t144 * t203;
t190 = qJD(3) * t243;
t185 = pkin(3) * t192;
t184 = -t27 * t115 + t22 * t142 + t53 * t205;
t181 = t145 * pkin(5) + t142 * qJ(6);
t180 = pkin(5) * t142 - qJ(6) * t145;
t179 = t113 * t53 - t232;
t13 = -pkin(5) * t112 + t211;
t178 = t13 * t145 - t142 * t16;
t177 = t13 * t142 + t145 * t16;
t176 = -t142 * t58 + t145 * t84;
t174 = 0.2e1 * qJD(3) * t125;
t123 = -pkin(4) - t181;
t173 = -t13 * t115 + t254;
t172 = t26 * t115 - t253;
t171 = t28 * t205 + t245;
t170 = t112 * t27 - t189;
t169 = t116 * t115 - t22;
t167 = t194 - t227;
t109 = t144 * t190;
t110 = t146 * t190;
t43 = t251 * qJD(4) - t247 * t109 - t143 * t110;
t47 = pkin(4) * t91 - pkin(9) * t90 + t199;
t165 = t142 * t47 + t145 * t43 + t79 * t205 - t82 * t206;
t163 = -t112 * t194 + t90 * t187 + t69;
t161 = (-t112 * t205 - t231) * pkin(9);
t159 = -t134 * t206 + t145 * t185;
t158 = -t224 + (-t219 + t229) * qJD(5);
t157 = t13 * t261 - t16 * t260 + t259;
t156 = (-t121 * t205 - t230) * t112 - t121 * t231;
t155 = t178 * qJD(5) + t259;
t154 = -t145 * t49 + t158;
t153 = t156 - t257;
t152 = t116 * t113 - t21;
t44 = t82 * qJD(4) - t143 * t109 + t247 * t110;
t148 = qJD(1) ^ 2;
t135 = -t201 - pkin(4);
t117 = -t201 + t123;
t108 = t115 * qJ(6);
t85 = t91 * t202;
t72 = -t113 ^ 2 + t115 ^ 2;
t64 = -t115 * t202 - t262;
t63 = t113 * t202 + t150;
t62 = -pkin(5) * t162 + qJ(6) * t99;
t46 = t180 * t121 - t251;
t34 = pkin(5) * t164 + t142 * t82 - t145 * t79;
t33 = -qJ(6) * t164 + t238;
t31 = -t48 + t233;
t30 = -t176 + t246;
t29 = -t108 + t241;
t25 = t142 * t66 - t145 * t76 + t246;
t24 = -t108 + t240;
t15 = t112 * t187 - t115 * t162 + t231;
t14 = -t99 * t115 - t249 * t142 + t228;
t11 = -t162 * t187 - t224;
t9 = t180 * t90 + (t181 * qJD(5) - qJD(6) * t145) * t121 + t44;
t8 = -t91 * pkin(5) + t238 * qJD(5) + t142 * t43 - t145 * t47;
t7 = qJ(6) * t91 - qJD(6) * t164 + t165;
t5 = (-t48 - t233) * t145 + (-t49 + t220) * t142;
t1 = [0, 0, 0, 0, 0.2e1 * t146 * t191, -0.2e1 * t210 * t203, t212, -t213, 0, -t131 * t212 + t144 * t174, t131 * t213 + t146 * t174, -t115 * t90 + t121 * t150, -t90 * t113 + t115 * t91 - t121 * t262 + t150 * t164, t256, -t85, 0, t122 * t262 + t116 * t91 - t44 * t202 + (-qJD(1) * t164 + t113) * t199, pkin(3) * t121 * t191 - t115 * t199 + t116 * t90 + t122 * t150 - t43 * t202, t162 * t167 - t48 * t214, t198 + (t224 + (t219 + t229) * qJD(5)) * t121 + t242, t163 + t239, t156 + t257, t112 * t91 - t164 * t262, t189 * t164 + t26 * t91 + t44 * t99 - t251 * t49 + ((-qJD(5) * t82 + t47) * t112 + t79 * t262 + t53 * qJD(5) * t121) * t145 + ((-qJD(5) * t79 - t43) * t112 - t82 * t262 + t22 * t121 + t53 * t90) * t142, -t165 * t112 + t253 * t121 - t162 * t44 + t164 * t166 + t53 * t227 - t238 * t262 + t251 * t48 - t27 * t91, -t8 * t112 + t121 * t171 - t13 * t91 + t164 * t4 + t28 * t230 - t262 * t34 + t46 * t49 + t9 * t99, -t8 * t162 - t33 * t49 - t34 * t48 - t7 * t99 + t178 * t90 + (-qJD(5) * t177 - t142 * t2 + t145 * t4) * t121, t7 * t112 + t254 * t121 + t16 * t91 + t162 * t9 - t164 * t2 - t28 * t227 + t262 * t33 + t46 * t48, t13 * t8 + t16 * t7 + t2 * t33 + t28 * t9 + t34 * t4 + t46 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, -t212, 0, 0, 0, 0, 0, -t85, -t256, 0, 0, 0, 0, 0, t153, t112 * t167 + t239 - t69, t153, t121 * t158 - t198 + t242, t163 - t239, t121 * t155 - t164 * t6 + t177 * t90 + t28 * t91; 0, 0, 0, 0, -t144 * t148 * t146, t210 * t148, 0, 0, 0, -t125 * t209, -t125 * t208, -t215, t72, t63, t64, 0, t65 * t202 + (-t113 * t209 - t202 * t207) * pkin(3) + t169, t66 * t202 + (t115 * t209 - t202 * t192) * pkin(3) + t152, t11, t5, t15, t14, t218, t135 * t49 + t183 * t99 + t179 * t142 + ((-qJD(5) * t134 - t76) * t145 + (-t185 + t66) * t142) * t112 + t172, -t135 * t48 + t179 * t145 - t183 * t162 + (-t159 + t240) * t112 + t184, t117 * t49 + t222 * t99 + (t113 * t28 - t232) * t142 + (-t134 * t205 - t142 * t185 + t25) * t112 + t173, t25 * t162 + t24 * t99 + (-t162 * t197 - t196 * t99) * t237 + t154 * t134 + t157, t117 * t48 - t245 + (-qJD(5) * t28 + t232) * t145 + t222 * t162 + (-t24 + t159) * t112 + t255, t6 * t117 - t13 * t25 - t16 * t24 + t222 * t28 + (t13 * t197 + t16 * t196) * t237 + t155 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t72, t63, t64, 0, t59 * t202 + t169, t58 * t202 + t152, t11, t5, t15, t14, t218, -pkin(4) * t49 - t176 * t112 + t53 * t217 - t59 * t99 + t161 + t172, pkin(4) * t48 + t241 * t112 + t162 * t59 + t53 * t216 + t184 + t258, t30 * t112 + t123 * t49 + t217 * t28 + t221 * t99 + t161 + t173, pkin(9) * t154 + t162 * t30 + t29 * t99 + t157, -t29 * t112 + t123 * t48 + t162 * t221 - t171 + t255 - t258, pkin(9) * t155 + t6 * t123 - t13 * t30 - t16 * t29 + t221 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t99 ^ 2 + t250, t31, -t49 - t220, t262, t162 * t53 + t170, t112 * t26 + t53 * t99 - t166, -t62 * t99 + t170 + t235 + 0.2e1 * t248, pkin(5) * t48 - t49 * qJ(6) - (t16 - t27) * t162 + (t13 - t211) * t99, 0.2e1 * t223 - t62 * t162 - t28 * t99 + (0.2e1 * qJD(6) - t26) * t112 + t166, -t4 * pkin(5) + t2 * qJ(6) - t13 * t27 + t16 * t211 - t28 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262 - t234, t31, -t249 - t250, -t112 * t16 - t235 + t4;];
tauc_reg  = t1;
