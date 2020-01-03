% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m_mdh [6x1]
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:05
% EndTime: 2020-01-03 11:20:21
% DurationCPUTime: 6.56s
% Computational Cost: add. (8366->408), mult. (8328->550), div. (0->0), fcn. (7746->10), ass. (0->206)
t177 = qJ(1) + pkin(7);
t172 = cos(t177);
t170 = sin(t177);
t277 = qJ(3) * t170;
t130 = pkin(2) * t172 + t277;
t175 = cos(qJ(1)) * pkin(1);
t168 = qJD(1) * t175;
t253 = qJD(1) * t172;
t254 = qJD(1) * t170;
t258 = pkin(2) * t253 + qJ(3) * t254;
t322 = -qJD(1) * t130 + t168 + t258;
t176 = pkin(9) + qJ(5);
t169 = sin(t176);
t171 = cos(t176);
t179 = sin(pkin(8));
t181 = cos(pkin(8));
t271 = t172 * t181;
t103 = t169 * t271 - t170 * t171;
t104 = t169 * t170 + t171 * t271;
t272 = t172 * t179;
t58 = Icges(6,5) * t104 - Icges(6,6) * t103 + Icges(6,3) * t272;
t91 = Icges(6,4) * t104;
t61 = -Icges(6,2) * t103 + Icges(6,6) * t272 + t91;
t90 = Icges(6,4) * t103;
t65 = -Icges(6,1) * t104 - Icges(6,5) * t272 + t90;
t24 = (t169 * t61 + t171 * t65) * t179 + t181 * t58;
t321 = rSges(5,3) + qJ(4);
t274 = t170 * t181;
t101 = -t169 * t274 - t171 * t172;
t102 = -t169 * t172 + t171 * t274;
t275 = t170 * t179;
t57 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t275;
t280 = Icges(6,4) * t102;
t60 = Icges(6,2) * t101 + Icges(6,6) * t275 + t280;
t89 = Icges(6,4) * t101;
t63 = Icges(6,1) * t102 + Icges(6,5) * t275 + t89;
t14 = t101 * t60 + t102 * t63 + t275 * t57;
t211 = -t103 * t61 - t104 * t65;
t320 = -t211 + t14;
t299 = pkin(3) * t181;
t213 = qJ(4) * t179 + t299;
t124 = t213 * t172;
t249 = qJD(4) * t179;
t143 = t170 * t249;
t319 = -qJD(1) * t124 + t143 + t322;
t252 = qJD(1) * t179;
t238 = t172 * t252;
t250 = qJD(3) * t172;
t180 = cos(pkin(9));
t268 = t180 * t181;
t178 = sin(pkin(9));
t270 = t178 * t181;
t276 = t170 * t178;
t310 = -rSges(5,2) * (-t170 * t180 + t172 * t270) + rSges(5,1) * (t172 * t268 + t276);
t318 = rSges(5,3) * t238 + qJD(1) * t310 - t250;
t316 = t101 * t61 - t102 * t65;
t100 = -rSges(6,3) * t181 + (rSges(6,1) * t171 - rSges(6,2) * t169) * t179;
t247 = qJD(5) * t179;
t313 = t100 * t247;
t161 = -qJD(5) * t181 + qJD(1);
t97 = -Icges(6,3) * t181 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t179;
t278 = Icges(6,4) * t171;
t98 = -Icges(6,6) * t181 + (-Icges(6,2) * t169 + t278) * t179;
t279 = Icges(6,4) * t169;
t99 = -Icges(6,5) * t181 + (Icges(6,1) * t171 - t279) * t179;
t189 = -t103 * t98 + t104 * t99 + t272 * t97;
t312 = t189 * t161;
t311 = rSges(4,1) * t274 - rSges(4,2) * t275;
t273 = t172 * t178;
t309 = -(t170 * t268 - t273) * rSges(5,1) - (-t170 * t270 - t172 * t180) * rSges(5,2);
t67 = rSges(6,1) * t104 - rSges(6,2) * t103 + rSges(6,3) * t272;
t308 = t161 * t67 + t143 + t168;
t256 = rSges(3,1) * t170 + rSges(3,2) * t172;
t289 = rSges(4,3) * t172;
t303 = t289 - t311;
t302 = pkin(3) * t274 + t275 * t321 - t309;
t190 = t170 * (-Icges(6,2) * t102 + t63 + t89) - t172 * (Icges(6,2) * t104 + t65 + t90);
t301 = t170 * (-Icges(6,1) * t101 + t280 + t60) - t172 * (-Icges(6,1) * t103 - t61 - t91);
t185 = qJD(1) ^ 2;
t74 = qJD(1) * t101 + qJD(5) * t104;
t75 = qJD(1) * t102 + qJD(5) * t103;
t220 = rSges(6,1) * t75 + rSges(6,2) * t74;
t239 = t170 * t252;
t42 = rSges(6,3) * t239 + t220;
t76 = -qJD(1) * t103 - qJD(5) * t102;
t77 = qJD(1) * t104 + qJD(5) * t101;
t43 = rSges(6,1) * t77 + rSges(6,2) * t76 + rSges(6,3) * t238;
t208 = t170 * t42 + t172 * t43;
t66 = rSges(6,1) * t102 + rSges(6,2) * t101 + rSges(6,3) * t275;
t9 = ((-t170 * t66 - t172 * t67) * qJD(1) + t208) * t247;
t300 = m(6) * t9;
t298 = pkin(4) * t178;
t174 = sin(qJ(1)) * pkin(1);
t297 = -t103 * t60 + t104 * t63;
t290 = rSges(4,3) * t170;
t265 = t124 + t130;
t166 = pkin(4) * t180 + pkin(3);
t263 = -pkin(4) * t276 - t166 * t271;
t182 = -pkin(6) - qJ(4);
t266 = qJ(4) + t182;
t73 = (t179 * t266 + t299) * t172 + t263;
t22 = (-qJD(3) - t313) * t172 + (-t73 + t265) * qJD(1) + t308;
t288 = t170 * t22;
t287 = t170 * t58;
t165 = t170 * pkin(2);
t129 = -qJ(3) * t172 + t165;
t202 = pkin(4) * t273 - t174;
t269 = t179 * t182;
t218 = t166 * t274 - t170 * t269;
t236 = t170 * t247;
t144 = t172 * t249;
t162 = qJD(3) * t170;
t259 = -t144 - t162;
t21 = -t100 * t236 + t161 * t66 + (t129 - t202 + t218) * qJD(1) + t259;
t286 = t172 * t21;
t285 = t172 * t57;
t257 = qJ(3) * t253 + t162;
t105 = pkin(2) * t254 - t257;
t284 = -t213 * t254 - t105 + t144;
t120 = (-Icges(6,2) * t171 - t279) * t179;
t283 = t120 + t99;
t121 = (-Icges(6,1) * t169 - t278) * t179;
t282 = t121 - t98;
t167 = t185 * t175;
t281 = t167 + qJD(1) * (-t250 + t258);
t251 = qJD(1) * t181;
t237 = t172 * t251;
t243 = qJD(1) * t298;
t264 = t166 * t237 + t170 * t243;
t262 = pkin(3) * t237 + qJ(4) * t238;
t261 = rSges(4,1) * t237 + rSges(4,3) * t254;
t255 = t165 + t174;
t248 = qJD(5) * t100;
t15 = -t275 * t58 - t316;
t246 = t185 * t174;
t241 = t144 + t257;
t235 = qJD(1) * t249;
t233 = -t247 / 0.2e1;
t232 = t247 / 0.2e1;
t231 = t129 + t174;
t230 = t175 + t277;
t229 = t170 * t235 + qJD(1) * (t143 + t262) + t281;
t228 = t168 - t250;
t226 = t170 * t233;
t225 = t170 * t232;
t224 = t172 * t233;
t223 = t172 * t232;
t222 = qJD(1) * t162 - t246;
t221 = qJD(1) * t232;
t217 = t143 + t228;
t131 = rSges(3,1) * t172 - rSges(3,2) * t170;
t215 = rSges(4,1) * t181 - rSges(4,2) * t179;
t214 = t309 * qJD(1);
t194 = (-rSges(6,1) * t169 - rSges(6,2) * t171) * t179;
t111 = qJD(5) * t194;
t205 = t172 * t235 + t222;
t10 = -t172 * t111 * t247 - t161 * t42 + (t172 * t243 + ((pkin(3) - t166) * t251 + (qJD(1) * t266 + t248) * t179) * t170 + t284) * qJD(1) + t205;
t11 = -t111 * t236 + t161 * t43 + ((-qJD(3) + (-qJD(1) * t182 - t248) * t179) * t172 - t262 + t264) * qJD(1) + t229;
t212 = -t10 * t172 - t11 * t170;
t210 = t14 * t172 + t15 * t170;
t209 = -t170 * t21 - t172 * t22;
t207 = -t170 * t67 + t172 * t66;
t206 = t170 * (Icges(6,5) * t101 - Icges(6,6) * t102) - t172 * (Icges(6,5) * t103 + Icges(6,6) * t104);
t204 = t170 * t221;
t203 = t172 * t221;
t201 = pkin(2) + t215;
t199 = rSges(5,3) * t272 + t310;
t193 = (t14 * t170 - t15 * t172) * t179;
t16 = -t272 * t57 - t297;
t17 = t272 * t58 + t211;
t192 = (t16 * t170 - t17 * t172) * t179;
t119 = (-Icges(6,5) * t169 - Icges(6,6) * t171) * t179;
t36 = Icges(6,5) * t75 + Icges(6,6) * t74 + Icges(6,3) * t239;
t37 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t238;
t38 = Icges(6,4) * t75 + Icges(6,2) * t74 + Icges(6,6) * t239;
t39 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t238;
t40 = Icges(6,1) * t75 + Icges(6,4) * t74 + Icges(6,5) * t239;
t41 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t238;
t188 = ((t103 * t39 - t104 * t41 + t60 * t74 + t63 * t75 + (-t172 * t37 + t254 * t57) * t179) * t170 - t172 * (t103 * t38 - t104 * t40 - t61 * t74 + t65 * t75 + (-t172 * t36 - t254 * t58) * t179) + (t16 * t172 + t17 * t170) * qJD(1)) * t179;
t187 = (qJD(1) * t210 + t170 * (t101 * t39 + t102 * t41 + t60 * t76 + t63 * t77 + (t170 * t37 + t253 * t57) * t179) - t172 * (t101 * t38 + t102 * t40 - t61 * t76 + t65 * t77 + (t170 * t36 - t253 * t58) * t179)) * t179;
t23 = -t181 * t57 + (-t169 * t60 + t171 * t63) * t179;
t7 = -t181 * t37 + (-t169 * t39 + t171 * t41 + (-t169 * t63 - t171 * t60) * qJD(5)) * t179;
t8 = -t181 * t36 + (-t169 * t38 + t171 * t40 + (-t169 * t65 + t171 * t61) * qJD(5)) * t179;
t186 = (t170 * t7 - t172 * t8 + (t24 * t170 + t172 * t23) * qJD(1)) * t179;
t110 = qJD(5) * t121;
t109 = qJD(5) * t120;
t108 = qJD(5) * t119;
t88 = t172 * t215 + t290;
t85 = rSges(6,1) * t103 + rSges(6,2) * t104;
t84 = rSges(6,1) * t101 - rSges(6,2) * t102;
t70 = (t130 + t88) * qJD(1) + t228;
t53 = ((-rSges(4,2) * t252 - qJD(3)) * t172 + t261) * qJD(1) + t281;
t52 = (qJD(1) * t303 - t105) * qJD(1) + t222;
t45 = (t199 + t265) * qJD(1) + t217;
t31 = t101 * t98 + t102 * t99 + t275 * t97;
t30 = t31 * t161;
t29 = qJD(1) * t318 + t229;
t28 = (-rSges(5,3) * t239 + t214 + t284) * qJD(1) + t205;
t27 = -qJD(4) * t181 + t207 * t247 + qJD(2);
t26 = -t108 * t181 + (-t109 * t169 + t110 * t171 + (-t169 * t99 - t171 * t98) * qJD(5)) * t179;
t25 = t26 * t161;
t13 = t101 * t109 + t102 * t110 + t76 * t98 + t77 * t99 + (t108 * t170 + t253 * t97) * t179;
t12 = t103 * t109 - t104 * t110 + t74 * t98 + t75 * t99 + (-t108 * t172 + t254 * t97) * t179;
t6 = qJD(5) * t192 - t312;
t5 = qJD(5) * t193 + t30;
t1 = [(t30 + ((t17 + t320) * t170 + (t16 + (t285 - t287) * t179 - t15 + t297) * t172) * t247) * t223 + t25 + m(3) * ((-t185 * t256 - t246) * (t131 + t175) + (-t167 + (-0.2e1 * rSges(3,1) * t253 + 0.2e1 * rSges(3,2) * t254 + qJD(1) * t131) * qJD(1)) * (-t256 - t174)) + (t7 + t13) * t225 + (t24 - t189) * t204 + (t23 + t31) * t203 + (t10 * (t230 + t67 - t263) + t22 * (-t220 + t241) + t11 * (t218 + t66 + t255) + (t22 * t202 - t269 * t286 + (-t166 * t181 - pkin(2) + (-rSges(6,3) + t182) * t179) * t288) * qJD(1) + (qJD(1) * t73 + t22 + t264 - t308 + t319 + t43) * t21 + (t10 * (pkin(2) - t269) + t11 * (-qJ(3) - t298) + t21 * t313) * t172) * m(6) + (t28 * (t199 + t230) + t29 * (t255 + t302) + (t28 * (pkin(2) + t213) - t29 * qJ(3)) * t172 + (t214 + t241 + (-t174 + (-t179 * t321 - pkin(2) - t299) * t170) * qJD(1)) * t45 + (-qJD(1) * t199 - t217 + t262 + t318 + t319 + t45) * ((t231 + t302) * qJD(1) + t259)) * m(5) + (t52 * (t230 + t290) + t53 * (t255 + t311) + (t52 * t201 + t53 * (-rSges(4,3) - qJ(3))) * t172 + (t257 + (-t170 * t201 - t174 + t289) * qJD(1)) * t70 + (-t250 - t228 + t261 + t70 + (-rSges(4,2) * t272 - t88) * qJD(1) + t322) * (-t162 + (t231 - t303) * qJD(1))) * m(4) + (t6 + t312 + ((t297 + t316) * t170 - t320 * t172 + ((t285 + t287) * t170 + t172 ^ 2 * t58) * t179 + t210) * t247) * t226 + (t8 + t5 + t12) * t224; t300; m(4) * (-t170 * t53 - t172 * t52) + m(5) * (-t170 * t29 - t172 * t28) + m(6) * t212; -t181 * t300 + 0.2e1 * (m(5) * (t170 * t28 - t172 * t29) / 0.2e1 + m(6) * (t10 * t170 - t11 * t172) / 0.2e1) * t179; -t181 * (qJD(5) * t186 + t25) / 0.2e1 + t161 * (-t26 * t181 + t186) / 0.2e1 + (qJD(5) * t187 + t13 * t161) * t275 / 0.2e1 + (-t181 * t31 + t193) * t203 + (-t13 * t181 + t187) * t225 - (qJD(5) * t188 + t12 * t161) * t272 / 0.2e1 + (t181 * t189 + t192) * t204 + (-t12 * t181 + t188) * t224 - t161 * (-t181 * t119 * t161 + ((-t169 * t283 + t171 * t282) * t161 + ((-t169 * t190 - t171 * t301) * t179 - t206 * t181) * qJD(5)) * t179) / 0.2e1 + ((t101 * t283 + t102 * t282 + t119 * t275) * t161 + (t101 * t190 - t102 * t301 + t206 * t275) * t247) * t226 + ((t103 * t283 - t104 * t282 - t119 * t272) * t161 + (t190 * t103 + t104 * t301 - t206 * t272) * t247) * t223 + (t170 * t6 + t172 * t5) * t252 / 0.2e1 + ((-t10 * t67 - t11 * t66 - t21 * t43 + t22 * t42) * t181 + (t9 * t207 + t27 * (-t253 * t67 - t254 * t66 + t208) + t209 * t111 + ((-t286 + t288) * qJD(1) + t212) * t100) * t179 - (t21 * t84 - t22 * t85) * t161 - (t27 * (t170 * t85 + t172 * t84) + t209 * t194) * t247) * m(6);];
tauc = t1(:);
