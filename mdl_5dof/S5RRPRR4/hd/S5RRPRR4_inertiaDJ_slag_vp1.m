% Calculate time derivative of joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:04
% EndTime: 2022-01-20 10:48:13
% DurationCPUTime: 3.20s
% Computational Cost: add. (9289->348), mult. (5874->496), div. (0->0), fcn. (4392->10), ass. (0->216)
t171 = qJ(1) + qJ(2);
t164 = sin(t171);
t166 = cos(t171);
t124 = t166 * rSges(3,1) - rSges(3,2) * t164;
t162 = pkin(9) + t171;
t158 = sin(t162);
t159 = cos(t162);
t174 = cos(qJ(4));
t236 = qJD(4) * t174;
t169 = qJD(1) + qJD(2);
t172 = sin(qJ(4));
t245 = t169 * t172;
t284 = -t158 * t236 - t159 * t245;
t261 = rSges(5,2) * t172;
t264 = rSges(5,1) * t174;
t283 = -t261 + t264;
t255 = Icges(5,4) * t174;
t202 = -Icges(5,2) * t172 + t255;
t193 = t202 * t159;
t88 = Icges(5,6) * t158 + t193;
t256 = Icges(5,4) * t172;
t204 = Icges(5,1) * t174 - t256;
t195 = t204 * t159;
t90 = Icges(5,5) * t158 + t195;
t205 = t172 * t88 - t174 * t90;
t282 = t205 * t158;
t87 = -Icges(5,6) * t159 + t202 * t158;
t89 = -Icges(5,5) * t159 + t204 * t158;
t207 = t172 * t87 - t174 * t89;
t281 = t207 * t159;
t170 = qJ(4) + qJ(5);
t163 = sin(t170);
t165 = cos(t170);
t253 = Icges(6,4) * t165;
t201 = -Icges(6,2) * t163 + t253;
t192 = t201 * t159;
t76 = Icges(6,6) * t158 + t192;
t254 = Icges(6,4) * t163;
t203 = Icges(6,1) * t165 - t254;
t194 = t203 * t159;
t78 = Icges(6,5) * t158 + t194;
t211 = t163 * t76 - t165 * t78;
t280 = t211 * t158;
t75 = -Icges(6,6) * t159 + t201 * t158;
t77 = -Icges(6,5) * t159 + t203 * t158;
t212 = t163 * t75 - t165 * t77;
t279 = t212 * t159;
t147 = t158 * rSges(6,3);
t263 = rSges(6,1) * t165;
t278 = t159 * t263 + t147;
t160 = pkin(2) * t166;
t277 = -t159 * rSges(4,1) - t160;
t168 = qJD(4) + qJD(5);
t120 = Icges(6,2) * t165 + t254;
t121 = Icges(6,1) * t163 + t253;
t198 = t120 * t163 - t121 * t165;
t199 = Icges(6,5) * t165 - Icges(6,6) * t163;
t276 = t199 * t168 + t198 * t169;
t143 = Icges(5,2) * t174 + t256;
t144 = Icges(5,1) * t172 + t255;
t197 = t143 * t172 - t144 * t174;
t200 = Icges(5,5) * t174 - Icges(5,6) * t172;
t275 = t200 * qJD(4) + t197 * t169;
t274 = 2 * m(3);
t273 = 2 * m(4);
t272 = 2 * m(5);
t271 = 2 * m(6);
t270 = t158 / 0.2e1;
t269 = -t159 / 0.2e1;
t132 = t283 * qJD(4);
t268 = m(5) * t132;
t146 = rSges(5,1) * t172 + rSges(5,2) * t174;
t267 = m(5) * t146;
t266 = pkin(2) * t164;
t152 = t158 * pkin(7);
t173 = sin(qJ(1));
t265 = t173 * pkin(1);
t260 = rSges(6,2) * t163;
t130 = t158 * t260;
t241 = t159 * rSges(6,3) + t130;
t79 = t158 * t263 - t241;
t230 = t159 * t260;
t80 = -t230 + t278;
t39 = t158 * t79 + t159 * t80;
t259 = pkin(1) * qJD(1);
t148 = t158 * rSges(5,3);
t252 = t120 * t168;
t251 = t121 * t168;
t250 = t158 * t169;
t249 = t159 * t169;
t176 = -pkin(8) - pkin(7);
t248 = t159 * t176;
t247 = t163 * t168;
t246 = t165 * t168;
t244 = t169 * t176;
t243 = rSges(6,3) * t249 + t169 * t130;
t237 = qJD(4) * t172;
t225 = t158 * t237;
t242 = -pkin(4) * t225 - t158 * t244;
t240 = t159 * rSges(5,3) + t158 * t261;
t239 = -t159 * pkin(3) - t152;
t238 = t158 ^ 2 + t159 ^ 2;
t234 = rSges(6,2) * t246;
t227 = -t169 * t230 + (-t247 * rSges(6,1) - t234) * t158;
t235 = t158 * (t278 * t169 + t227) + t159 * (-t159 * t234 + (-t159 * t247 - t165 * t250) * rSges(6,1) + t243) + t79 * t249;
t233 = pkin(4) * t237;
t232 = t173 * t259;
t175 = cos(qJ(1));
t231 = t175 * t259;
t229 = t158 * t245;
t226 = -rSges(5,1) * t225 + t284 * rSges(5,2);
t223 = t250 / 0.2e1;
t222 = t249 / 0.2e1;
t221 = -pkin(3) - t264;
t122 = rSges(6,1) * t163 + rSges(6,2) * t165;
t220 = -pkin(4) * t172 - t122;
t186 = Icges(6,6) * t169 - t252;
t46 = t186 * t159 - t201 * t250;
t219 = t168 * t78 + t46;
t47 = t186 * t158 + t169 * t192;
t218 = t168 * t77 + t47;
t187 = Icges(6,5) * t169 - t251;
t48 = t187 * t159 - t203 * t250;
t217 = -t168 * t76 + t48;
t49 = t187 * t158 + t169 * t194;
t216 = -t168 * t75 + t49;
t161 = pkin(4) * t174 + pkin(3);
t215 = -t161 - t263;
t214 = -t158 * t176 + t159 * t161;
t107 = t124 * t169;
t100 = -rSges(4,2) * t158 - t277;
t73 = -Icges(6,3) * t159 + t199 * t158;
t17 = -t212 * t158 - t159 * t73;
t190 = t199 * t159;
t74 = Icges(6,3) * t158 + t190;
t18 = -t159 * t74 - t280;
t19 = t158 * t73 - t279;
t20 = t158 * t74 - t211 * t159;
t119 = Icges(6,5) * t163 + Icges(6,6) * t165;
t185 = Icges(6,3) * t169 - t119 * t168;
t44 = t185 * t159 - t199 * t250;
t45 = t185 * t158 + t169 * t190;
t213 = -t159 * ((t159 * t45 + (t18 + t279) * t169) * t159 + (t17 * t169 + (-t163 * t46 + t165 * t48 - t76 * t246 - t78 * t247) * t158 + (-t44 + (t169 * t78 - t216) * t165 + (-t169 * t76 + t218) * t163) * t159) * t158) + t158 * ((t158 * t44 + (t19 + t280) * t169) * t158 + (t20 * t169 + (t163 * t47 - t165 * t49 + t75 * t246 + t77 * t247) * t159 + (-t45 + (t169 * t77 + t217) * t165 + (-t169 * t75 - t219) * t163) * t158) * t159) + (t158 * t18 - t159 * t17) * t250 + (t158 * t20 - t159 * t19) * t249;
t123 = -rSges(3,1) * t164 - rSges(3,2) * t166;
t208 = t172 * t89 + t174 * t87;
t206 = t172 * t90 + t174 * t88;
t142 = Icges(5,5) * t172 + Icges(5,6) * t174;
t92 = t283 * t159 + t148;
t106 = t123 * t169;
t103 = t201 * t168;
t104 = t203 * t168;
t178 = t119 * t169 + (t104 - t252) * t165 + (-t103 - t251) * t163;
t196 = (t276 * t158 + t178 * t159 + t217 * t163 + t219 * t165) * t270 + (t178 * t158 - t276 * t159 + t216 * t163 + t218 * t165) * t269 + (-t159 * t119 - t198 * t158 + t163 * t77 + t165 * t75) * t223 + (t158 * t119 - t198 * t159 + t163 * t78 + t165 * t76) * t222;
t191 = t200 * t159;
t99 = -rSges(4,1) * t158 - rSges(4,2) * t159 - t266;
t189 = t221 * t158 - t266;
t82 = rSges(4,2) * t250 + t277 * t169;
t188 = t215 * t158 - t266;
t65 = t160 + t92 - t239;
t81 = t99 * t169;
t184 = Icges(5,5) * t169 - t144 * qJD(4);
t183 = Icges(5,6) * t169 - t143 * qJD(4);
t182 = Icges(5,3) * t169 - t142 * qJD(4);
t61 = t160 + t214 + t80;
t181 = -t146 * t159 * qJD(4) + rSges(5,2) * t229 + rSges(5,3) * t249;
t153 = t159 * pkin(7);
t64 = t153 + t189 + t240;
t60 = t188 + t241 - t248;
t127 = t202 * qJD(4);
t128 = t204 * qJD(4);
t180 = t165 * t103 + t163 * t104 - t120 * t247 + t121 * t246 + t174 * t127 + t172 * t128 - t143 * t237 + t144 * t236;
t177 = -t127 * t172 + t128 * t174 + t142 * t169 + (-t143 * t174 - t144 * t172) * qJD(4);
t179 = t196 + (-t205 * qJD(4) + t275 * t158 + t177 * t159 + t172 * (t184 * t159 - t204 * t250) + t174 * (t183 * t159 - t202 * t250)) * t270 + (-t207 * qJD(4) + t177 * t158 - t275 * t159 + t172 * (t184 * t158 + t169 * t195) + t174 * (t183 * t158 + t169 * t193)) * t269 + (-t159 * t142 - t197 * t158 + t208) * t223 + (t158 * t142 - t197 * t159 + t206) * t222;
t32 = (-t160 + t221 * t159 + (-rSges(5,3) - pkin(7)) * t158) * t169 - t226;
t24 = (t215 * t159 - t147 - t160) * t169 - t227 - t242;
t139 = pkin(7) * t249;
t31 = t189 * t169 + t139 + t181;
t23 = t188 * t169 + (-t122 * t168 - t233 - t244) * t159 + t243;
t167 = t175 * pkin(1);
t109 = t124 + t167;
t108 = t123 - t265;
t105 = (-t260 + t263) * t168;
t96 = -t107 - t231;
t95 = t106 - t232;
t94 = t100 + t167;
t93 = t99 - t265;
t91 = t158 * t264 - t240;
t86 = Icges(5,3) * t158 + t191;
t85 = -Icges(5,3) * t159 + t200 * t158;
t84 = t220 * t159;
t83 = t220 * t158;
t72 = t214 + t239;
t71 = t248 + t153 + t158 * (-pkin(3) + t161);
t70 = t82 - t231;
t69 = t81 - t232;
t63 = t167 + t65;
t62 = t64 - t265;
t55 = t182 * t158 + t169 * t191;
t54 = t182 * t159 - t200 * t250;
t53 = t167 + t61;
t52 = t60 - t265;
t43 = t284 * pkin(4) - t158 * t105 - t122 * t249;
t42 = t122 * t250 - t159 * t105 + (-t159 * t236 + t229) * pkin(4);
t30 = t32 - t231;
t29 = t31 - t232;
t28 = t158 * t86 - t205 * t159;
t27 = t158 * t85 - t281;
t26 = -t159 * t86 - t282;
t25 = -t207 * t158 - t159 * t85;
t22 = t24 - t231;
t21 = t23 - t232;
t16 = t158 * t71 + t159 * t72 + t39;
t11 = ((-t92 + t148) * t169 + t226) * t158 + (t169 * t91 + t181) * t159;
t8 = -t80 * t250 + t235;
t3 = t158 * t242 + t159 * (-t159 * t233 - t139) + ((t71 - t248) * t159 + (-t72 - t80 - t152) * t158) * t169 + t235;
t1 = [(t21 * t53 + t22 * t52) * t271 + (t29 * t63 + t30 * t62) * t272 + (t69 * t94 + t70 * t93) * t273 + (t108 * t96 + t109 * t95) * t274 + t180; m(6) * (t21 * t61 + t22 * t60 + t23 * t53 + t24 * t52) + m(5) * (t29 * t65 + t30 * t64 + t31 * t63 + t32 * t62) + m(4) * (t100 * t69 + t70 * t99 + t81 * t94 + t82 * t93) + m(3) * (t106 * t109 - t107 * t108 + t123 * t96 + t124 * t95) + t180; (t23 * t61 + t24 * t60) * t271 + (t31 * t65 + t32 * t64) * t272 + (t100 * t81 + t82 * t99) * t273 + (t106 * t124 - t107 * t123) * t274 + t180; 0; 0; 0; ((-t169 * t63 - t30) * t159 + (t169 * t62 - t29) * t158) * t267 + m(6) * (t21 * t83 + t22 * t84 + t42 * t52 + t43 * t53) + (-t158 * t63 - t159 * t62) * t268 + t179; m(6) * (t23 * t83 + t24 * t84 + t42 * t60 + t43 * t61) + ((-t169 * t65 - t32) * t159 + (t169 * t64 - t31) * t158) * t267 + (-t158 * t65 - t159 * t64) * t268 + t179; m(5) * t11 + m(6) * t3; ((t158 * t91 + t159 * t92) * t11 + t238 * t146 * t132) * t272 + (t28 * t158 - t159 * t27) * t249 + t158 * ((t158 * t54 + (t27 + t282) * t169) * t158 + (t28 * t169 + (t87 * t236 + t89 * t237) * t159 + (-t206 * qJD(4) - t207 * t169 - t55) * t158) * t159) + (t158 * t26 - t159 * t25) * t250 - t159 * ((t159 * t55 + (t26 + t281) * t169) * t159 + (t25 * t169 + (-t88 * t236 - t90 * t237) * t158 + (t208 * qJD(4) - t205 * t169 - t54) * t159) * t158) + (t16 * t3 + t42 * t84 + t43 * t83) * t271 + t213; m(6) * ((-t158 * t53 - t159 * t52) * t105 + ((-t169 * t53 - t22) * t159 + (t169 * t52 - t21) * t158) * t122) + t196; m(6) * ((-t158 * t61 - t159 * t60) * t105 + ((-t169 * t61 - t24) * t159 + (t169 * t60 - t23) * t158) * t122) + t196; m(6) * t8; m(6) * (t16 * t8 + t3 * t39 + (-t158 * t83 - t159 * t84) * t105 + ((-t169 * t83 - t42) * t159 + (t169 * t84 - t43) * t158) * t122) + t213; (t238 * t122 * t105 + t39 * t8) * t271 + t213;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
