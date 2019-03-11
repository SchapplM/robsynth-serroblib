% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:27:47
% EndTime: 2019-03-09 08:27:58
% DurationCPUTime: 4.10s
% Computational Cost: add. (6885->373), mult. (17655->496), div. (0->0), fcn. (13212->8), ass. (0->184)
t181 = sin(pkin(9));
t184 = sin(qJ(2));
t186 = cos(qJ(2));
t242 = cos(pkin(9));
t161 = t181 * t186 + t242 * t184;
t148 = t161 * qJD(1);
t180 = sin(pkin(10));
t182 = cos(pkin(10));
t121 = t180 * qJD(2) + t182 * t148;
t183 = sin(qJ(5));
t185 = cos(qJ(5));
t234 = t180 * t148;
t273 = t182 * qJD(2) - t234;
t198 = t185 * t273;
t73 = -t183 * t121 + t198;
t279 = t73 ^ 2;
t213 = t242 * t186;
t171 = qJD(1) * t213;
t222 = t184 * qJD(1);
t145 = t181 * t222 - t171;
t141 = qJD(5) + t145;
t278 = t73 * t141;
t276 = t185 * t121 + t183 * t273;
t262 = t276 ^ 2;
t147 = t161 * qJD(2);
t137 = qJD(1) * t147;
t162 = t185 * t180 + t183 * t182;
t231 = t185 * t182;
t160 = t183 * t180 - t231;
t224 = qJD(5) * t185;
t225 = qJD(5) * t183;
t265 = -t180 * t225 + t182 * t224;
t272 = -t160 * t145 + t265;
t208 = t162 * t137 + t272 * t141;
t248 = t276 * t148;
t277 = t208 + t248;
t221 = qJD(1) * qJD(2);
t216 = t184 * t221;
t138 = qJD(2) * t171 - t181 * t216;
t235 = t180 * t138;
t40 = (qJD(5) * t121 + t235) * t183 - qJD(5) * t198 - t138 * t231;
t275 = t40 - t278;
t152 = t162 * qJD(5);
t243 = t162 * t145 + t152;
t207 = -t160 * t137 - t141 * t243;
t250 = t148 * t73;
t274 = t207 + t250;
t269 = -0.2e1 * t221;
t174 = t181 * pkin(2) + qJ(4);
t257 = pkin(8) + t174;
t156 = t257 * t180;
t157 = t257 * t182;
t104 = -t183 * t156 + t185 * t157;
t259 = t182 * pkin(8);
t256 = -qJ(3) - pkin(7);
t169 = t256 * t186;
t165 = qJD(1) * t169;
t153 = t181 * t165;
t168 = t256 * t184;
t164 = qJD(1) * t168;
t114 = t242 * t164 + t153;
t93 = pkin(2) * t222 + t148 * pkin(3) + t145 * qJ(4);
t55 = -t180 * t114 + t182 * t93;
t38 = t148 * pkin(4) + t145 * t259 + t55;
t238 = t145 * t180;
t56 = t182 * t114 + t180 * t93;
t48 = pkin(8) * t238 + t56;
t268 = -t162 * qJD(4) - t104 * qJD(5) + t183 * t48 - t185 * t38;
t199 = -t185 * t156 - t183 * t157;
t267 = t160 * qJD(4) - t199 * qJD(5) + t183 * t38 + t185 * t48;
t266 = t162 * t138;
t264 = t160 * t40 - t243 * t276;
t193 = -t181 * t184 + t213;
t219 = -t186 * pkin(2) - pkin(1);
t106 = -pkin(3) * t193 - t161 * qJ(4) + t219;
t119 = t181 * t168 - t242 * t169;
t59 = t182 * t106 - t180 * t119;
t47 = -pkin(4) * t193 - t161 * t259 + t59;
t236 = t161 * t180;
t60 = t180 * t106 + t182 * t119;
t52 = -pkin(8) * t236 + t60;
t201 = t183 * t47 + t185 * t52;
t150 = t193 * qJD(2);
t251 = qJD(2) * pkin(2);
t220 = t184 * t251;
t75 = t147 * pkin(3) - t150 * qJ(4) - t161 * qJD(4) + t220;
t215 = qJD(2) * t256;
t142 = t186 * qJD(3) + t184 * t215;
t143 = -t184 * qJD(3) + t186 * t215;
t95 = t242 * t142 + t181 * t143;
t43 = -t180 * t95 + t182 * t75;
t25 = t147 * pkin(4) - t150 * t259 + t43;
t237 = t150 * t180;
t44 = t180 * t75 + t182 * t95;
t33 = -pkin(8) * t237 + t44;
t263 = -t201 * qJD(5) - t183 * t33 + t185 * t25;
t144 = t145 ^ 2;
t261 = t137 * pkin(5);
t158 = t164 + t251;
t107 = t242 * t158 + t153;
t97 = -qJD(2) * pkin(3) + qJD(4) - t107;
t64 = -pkin(4) * t273 + t97;
t16 = -pkin(5) * t73 - qJ(6) * t276 + t64;
t260 = t16 * t276;
t258 = t276 * t73;
t255 = t148 * qJ(6) + t267;
t254 = -t148 * pkin(5) + t268;
t214 = t242 * t165;
t113 = t181 * t164 - t214;
t79 = -pkin(4) * t238 + t113;
t253 = -t243 * pkin(5) + t272 * qJ(6) + t162 * qJD(6) + t79;
t173 = pkin(2) * t216;
t63 = t137 * pkin(3) - t138 * qJ(4) - t148 * qJD(4) + t173;
t130 = t142 * qJD(1);
t131 = t143 * qJD(1);
t82 = t242 * t130 + t181 * t131;
t78 = qJD(2) * qJD(4) + t82;
t31 = t180 * t63 + t182 * t78;
t108 = t181 * t158 - t214;
t102 = qJD(2) * qJ(4) + t108;
t206 = t219 * qJD(1);
t167 = qJD(3) + t206;
t85 = t145 * pkin(3) - t148 * qJ(4) + t167;
t50 = t182 * t102 + t180 * t85;
t118 = -t242 * t168 - t181 * t169;
t81 = t181 * t130 - t242 * t131;
t246 = t81 * t118;
t49 = -t180 * t102 + t182 * t85;
t28 = t145 * pkin(4) - t121 * pkin(8) + t49;
t39 = pkin(8) * t273 + t50;
t9 = -t183 * t39 + t185 * t28;
t245 = qJD(6) - t9;
t241 = t199 * t137;
t240 = t104 * t137;
t239 = t137 * qJ(6);
t233 = t182 * t138;
t188 = qJD(1) ^ 2;
t230 = t186 * t188;
t187 = qJD(2) ^ 2;
t229 = t187 * t184;
t228 = t187 * t186;
t226 = t184 ^ 2 - t186 ^ 2;
t58 = pkin(4) * t235 + t81;
t30 = -t180 * t78 + t182 * t63;
t19 = t137 * pkin(4) - pkin(8) * t233 + t30;
t22 = -pkin(8) * t235 + t31;
t212 = t183 * t22 - t185 * t19 + t39 * t224 + t28 * t225;
t211 = pkin(1) * t269;
t94 = t181 * t142 - t242 * t143;
t41 = qJD(5) * t276 + t266;
t209 = -t162 * t41 + t272 * t73;
t177 = -t242 * pkin(2) - pkin(3);
t65 = pkin(4) * t237 + t94;
t91 = pkin(4) * t236 + t118;
t205 = -t180 * t49 + t182 * t50;
t10 = t183 * t28 + t185 * t39;
t202 = -t183 * t52 + t185 * t47;
t200 = t118 * t138 + t81 * t161;
t2 = t212 - t261;
t196 = t10 * t141 - t212;
t195 = t183 * t19 + t185 * t22 + t28 * t224 - t39 * t225;
t194 = t183 * t25 + t185 * t33 + t47 * t224 - t52 * t225;
t166 = -t182 * pkin(4) + t177;
t192 = t97 * t150 + t200;
t191 = -t174 * t137 + t177 * t138 + (-qJD(4) + t97) * t145;
t5 = t41 * pkin(5) + t40 * qJ(6) - qJD(6) * t276 + t58;
t189 = t141 * t276 + t41;
t99 = t160 * t161;
t98 = t162 * t161;
t96 = t160 * pkin(5) - t162 * qJ(6) + t166;
t54 = t162 * t150 + t265 * t161;
t53 = t160 * t150 + t161 * t152;
t37 = pkin(5) * t276 - qJ(6) * t73;
t32 = t98 * pkin(5) + t99 * qJ(6) + t91;
t15 = -t40 - t278;
t14 = pkin(5) * t193 - t202;
t13 = -qJ(6) * t193 + t201;
t8 = t54 * pkin(5) + t53 * qJ(6) + t99 * qJD(6) + t65;
t7 = t141 * qJ(6) + t10;
t6 = -t141 * pkin(5) + t245;
t4 = -t147 * pkin(5) - t263;
t3 = t147 * qJ(6) - qJD(6) * t193 + t194;
t1 = t141 * qJD(6) + t195 + t239;
t11 = [0, 0, 0, 0.2e1 * t186 * t216, t226 * t269, t228, -t229, 0, -pkin(7) * t228 + t184 * t211, pkin(7) * t229 + t186 * t211, -t107 * t150 - t108 * t147 - t119 * t137 - t95 * t145 + t94 * t148 + t193 * t82 + t200, -t107 * t94 + t108 * t95 + t246 + t82 * t119 + (t167 + t206) * t220, t59 * t137 + t43 * t145 + t49 * t147 + t192 * t180 - t193 * t30 - t273 * t94, t94 * t121 - t60 * t137 - t44 * t145 - t50 * t147 + t182 * t192 + t193 * t31, -t43 * t121 - t44 * t234 + (t44 * qJD(2) - t59 * t138 - t49 * t150 - t30 * t161) * t182 + (-t60 * t138 - t50 * t150 - t31 * t161) * t180, t30 * t59 + t31 * t60 + t49 * t43 + t50 * t44 + t97 * t94 + t246, -t276 * t53 + t40 * t99, -t276 * t54 + t40 * t98 + t41 * t99 - t53 * t73, -t99 * t137 - t53 * t141 + t147 * t276 + t193 * t40, -t98 * t137 - t54 * t141 + t147 * t73 + t193 * t41, -t137 * t193 + t141 * t147, t202 * t137 + t263 * t141 + t9 * t147 + t193 * t212 + t91 * t41 + t64 * t54 + t58 * t98 - t65 * t73, -t10 * t147 - t137 * t201 - t141 * t194 + t193 * t195 + t276 * t65 - t91 * t40 - t64 * t53 - t58 * t99, -t14 * t137 - t4 * t141 - t6 * t147 + t16 * t54 + t193 * t2 + t32 * t41 + t5 * t98 - t73 * t8, -t1 * t98 - t13 * t41 - t14 * t40 - t2 * t99 + t276 * t4 + t3 * t73 - t53 * t6 - t54 * t7, -t1 * t193 + t13 * t137 + t3 * t141 + t7 * t147 + t16 * t53 - t276 * t8 + t32 * t40 + t5 * t99, t1 * t13 + t14 * t2 + t16 * t8 + t3 * t7 + t32 * t5 + t4 * t6; 0, 0, 0, -t184 * t230, t226 * t188, 0, 0, 0, t188 * pkin(1) * t184, pkin(1) * t230 (t108 - t113) * t148 + (-t107 + t114) * t145 + (-t137 * t181 - t242 * t138) * pkin(2), t107 * t113 - t108 * t114 + (-t167 * t222 + t181 * t82 - t242 * t81) * pkin(2), t113 * t273 - t55 * t145 - t49 * t148 + t191 * t180 - t81 * t182, -t113 * t121 + t56 * t145 + t50 * t148 + t81 * t180 + t182 * t191, t55 * t121 + t56 * t234 + (-qJD(4) * t234 - t49 * t145 + t31 + (t182 * qJD(4) - t56) * qJD(2)) * t182 + (qJD(4) * t121 - t50 * t145 - t30) * t180, -t97 * t113 + t81 * t177 - t49 * t55 - t50 * t56 + (-t30 * t180 + t31 * t182) * t174 + t205 * qJD(4), -t40 * t162 + t272 * t276, t209 + t264, t208 - t248, t207 - t250, -t141 * t148, t268 * t141 - t9 * t148 + t58 * t160 + t166 * t41 + t243 * t64 + t73 * t79 + t241, t10 * t148 + t267 * t141 + t58 * t162 - t166 * t40 + t272 * t64 - t276 * t79 - t240, t141 * t254 + t6 * t148 + t16 * t243 + t5 * t160 + t253 * t73 + t96 * t41 + t241, -t1 * t160 - t104 * t41 + t2 * t162 + t199 * t40 - t243 * t7 - t254 * t276 - t255 * t73 + t272 * t6, -t141 * t255 - t7 * t148 - t16 * t272 - t5 * t162 + t253 * t276 + t96 * t40 + t240, t1 * t104 - t16 * t253 - t199 * t2 - t254 * t6 - t255 * t7 + t5 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 ^ 2 - t144, t107 * t148 + t108 * t145 + t173, t182 * t137 - t180 * t144 + t148 * t273, -t148 * t121 - t180 * t137 - t182 * t144 (t180 * t121 + t182 * t273) * t145 + (-t180 ^ 2 - t182 ^ 2) * t138, t145 * t205 - t97 * t148 + t31 * t180 + t30 * t182, 0, 0, 0, 0, 0, t274, -t277, t274, t209 - t264, t277, t1 * t162 - t16 * t148 + t2 * t160 + t243 * t6 + t272 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t145 + t235, t145 * t273 + t233, -t121 ^ 2 - t273 ^ 2, t121 * t49 - t273 * t50 + t81, 0, 0, 0, 0, 0, t189, -t275, t189, -t262 - t279, t275, -t276 * t6 - t7 * t73 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t262 - t279, t15, -t266 + (-qJD(5) + t141) * t276, t137, -t276 * t64 + t196, t9 * t141 - t64 * t73 - t195, t37 * t73 + t196 - t260 + 0.2e1 * t261, pkin(5) * t40 - qJ(6) * t41 + (-t10 + t7) * t276 - (t6 - t245) * t73, 0.2e1 * t239 + t16 * t73 + t37 * t276 + (0.2e1 * qJD(6) - t9) * t141 + t195, -pkin(5) * t2 + qJ(6) * t1 - t10 * t6 - t16 * t37 + t245 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t148 - t258, t15, -t141 ^ 2 - t262, -t7 * t141 + t2 + t260;];
tauc_reg  = t11;
