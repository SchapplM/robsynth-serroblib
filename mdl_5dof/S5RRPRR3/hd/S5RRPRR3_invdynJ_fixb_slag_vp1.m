% Calculate vector of inverse dynamics joint torques for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:59
% EndTime: 2022-01-20 10:34:07
% DurationCPUTime: 3.99s
% Computational Cost: add. (11057->369), mult. (6222->464), div. (0->0), fcn. (4592->10), ass. (0->222)
t182 = qJ(1) + qJ(2);
t174 = pkin(9) + t182;
t168 = qJ(4) + t174;
t162 = cos(t168);
t157 = t162 * pkin(8);
t161 = sin(t168);
t110 = pkin(4) * t161 - t157;
t181 = qJD(1) + qJD(2);
t173 = qJD(4) + t181;
t256 = t162 * t173;
t120 = pkin(8) * t256;
t183 = sin(qJ(5));
t258 = t161 * t183;
t134 = rSges(6,2) * t258;
t185 = cos(qJ(5));
t276 = rSges(6,2) * t185;
t241 = qJD(5) * t276;
t215 = rSges(6,3) * t256 + t173 * t134 - t162 * t241;
t245 = qJD(5) * t183;
t238 = t162 * t245;
t252 = t162 * rSges(6,3) + t134;
t257 = t161 * t185;
t79 = rSges(6,1) * t257 - t252;
t68 = t173 * t79;
t311 = -rSges(6,1) * t238 + t173 * t110 + t120 + t215 + t68;
t107 = rSges(5,1) * t161 + rSges(5,2) * t162;
t96 = t173 * t107;
t184 = sin(qJ(1));
t275 = pkin(1) * qJD(1);
t242 = t184 * t275;
t175 = sin(t182);
t176 = cos(t182);
t121 = rSges(3,1) * t175 + rSges(3,2) * t176;
t263 = t121 * t181;
t100 = -t242 - t263;
t278 = rSges(6,1) * t185;
t148 = -rSges(6,2) * t183 + t278;
t128 = t148 * qJD(5);
t146 = rSges(6,1) * t183 + t276;
t180 = qJDD(1) + qJDD(2);
t172 = qJDD(4) + t180;
t165 = sin(t174);
t166 = cos(t174);
t179 = t181 ^ 2;
t186 = cos(qJ(1));
t187 = qJD(1) ^ 2;
t211 = (-qJDD(1) * t184 - t186 * t187) * pkin(1);
t199 = t211 + (-t175 * t180 - t176 * t179) * pkin(2);
t189 = (-t165 * t180 - t166 * t179) * pkin(3) + t199;
t247 = qJD(5) * t162;
t268 = -t110 - t79;
t297 = -t162 * pkin(4) - t161 * pkin(8);
t255 = t162 * t183;
t244 = rSges(6,2) * t255;
t240 = t173 * t244 + (rSges(6,1) * t245 + t241) * t161;
t254 = t162 * t185;
t298 = rSges(6,1) * t254 + t161 * rSges(6,3);
t48 = t173 * t298 - t240;
t246 = qJD(5) * t173;
t89 = -qJDD(5) * t162 + t161 * t246;
t11 = -t128 * t247 + t89 * t146 + t268 * t172 + (t173 * t297 - t48) * t173 + t189;
t309 = -g(1) + t11;
t116 = rSges(4,1) * t165 + rSges(4,2) * t166;
t277 = rSges(4,2) * t165;
t137 = t181 * t277;
t159 = t166 * rSges(4,1);
t308 = -g(1) - t180 * t116 - t181 * (t159 * t181 - t137) + t199;
t164 = t176 * rSges(3,1);
t253 = t175 * t181;
t106 = -rSges(3,2) * t253 + t164 * t181;
t307 = -t106 * t181 - t121 * t180 - g(1) + t211;
t259 = t161 * t173;
t118 = rSges(5,2) * t259;
t85 = rSges(5,1) * t256 - t118;
t306 = t107 * t172 + t173 * t85 + g(1) - t189;
t160 = pkin(3) * t166;
t167 = pkin(2) * t176;
t178 = t186 * pkin(1);
t286 = pkin(1) * t184;
t229 = qJDD(1) * t178 - t187 * t286;
t216 = t180 * t167 + t229;
t285 = pkin(2) * t175;
t228 = -pkin(3) * t165 - t285;
t202 = t180 * t160 + t179 * t228 + t216;
t248 = qJD(5) * t161;
t47 = (-t173 * t257 - t238) * rSges(6,1) + t215;
t80 = -t244 + t298;
t62 = t80 - t297;
t88 = qJDD(5) * t161 + t162 * t246;
t12 = -t128 * t248 - t88 * t146 + (-pkin(4) * t259 + t120 + t47) * t173 + t62 * t172 + t202;
t305 = -g(2) + t12;
t155 = t162 * rSges(5,1);
t108 = -rSges(5,2) * t161 + t155;
t304 = t108 * t172 - t173 * t96 - g(2) + t202;
t102 = -t116 - t285;
t117 = t159 - t277;
t303 = t102 * t179 + t180 * t117 - g(2) + t216;
t122 = -rSges(3,2) * t175 + t164;
t302 = t122 * t180 - t181 * t263 - g(2) + t229;
t204 = t181 * t228 - t242;
t59 = t204 - t96;
t299 = t117 + t167;
t296 = t160 + t167;
t177 = Icges(6,4) * t185;
t219 = -Icges(6,2) * t183 + t177;
t143 = Icges(6,1) * t183 + t177;
t237 = -pkin(4) - t278;
t239 = t146 * t247;
t200 = t204 - t239;
t30 = t173 * t268 + t200;
t273 = t162 * t30;
t243 = t186 * t275;
t205 = t181 * t296 + t243;
t294 = t146 * t248 - t173 * t62;
t31 = t205 - t294;
t190 = (t237 * t273 + (t30 * (-rSges(6,3) - pkin(8)) + t31 * t237) * t161) * t173;
t295 = t190 + (t239 + t311) * t31 + (t240 - t294) * t30;
t140 = Icges(6,5) * t185 - Icges(6,6) * t183;
t139 = Icges(6,5) * t183 + Icges(6,6) * t185;
t206 = Icges(6,3) * t173 - qJD(5) * t139;
t213 = t219 * t162;
t76 = Icges(6,6) * t161 + t213;
t270 = t183 * t76;
t266 = Icges(6,4) * t183;
t144 = Icges(6,1) * t185 - t266;
t214 = t144 * t162;
t78 = Icges(6,5) * t161 + t214;
t221 = -t185 * t78 + t270;
t293 = -t140 * t259 + t162 * t206 + t173 * t221;
t212 = t140 * t162;
t75 = Icges(6,4) * t257 - Icges(6,2) * t258 - Icges(6,6) * t162;
t271 = t183 * t75;
t133 = Icges(6,4) * t258;
t77 = Icges(6,1) * t257 - Icges(6,5) * t162 - t133;
t222 = -t185 * t77 + t271;
t292 = t161 * t206 + (t212 + t222) * t173;
t141 = Icges(6,2) * t185 + t266;
t217 = t141 * t183 - t143 * t185;
t291 = t140 * qJD(5) + t173 * t217;
t73 = Icges(6,5) * t257 - Icges(6,6) * t258 - Icges(6,3) * t162;
t24 = -t161 * t222 - t162 * t73;
t280 = -Icges(6,2) * t257 - t133 + t77;
t282 = t143 * t161 + t75;
t290 = -t183 * t280 - t185 * t282;
t289 = t88 / 0.2e1;
t288 = t89 / 0.2e1;
t287 = m(4) + m(5);
t284 = -t161 * t73 - t77 * t254;
t74 = Icges(6,3) * t161 + t212;
t283 = t161 * t74 + t78 * t254;
t281 = -t143 * t162 - t76;
t279 = -t141 * t162 + t78;
t261 = t139 * t162;
t55 = -t161 * t217 - t261;
t269 = t55 * t173;
t264 = t108 * t173;
t262 = t139 * t161;
t260 = t140 * t173;
t251 = -t141 + t144;
t250 = t143 + t219;
t236 = -t248 / 0.2e1;
t235 = t248 / 0.2e1;
t234 = -t247 / 0.2e1;
t233 = t247 / 0.2e1;
t63 = t78 * t257;
t231 = t162 * t74 - t63;
t230 = -t73 + t270;
t149 = rSges(2,1) * t186 - rSges(2,2) * t184;
t147 = rSges(2,1) * t184 + rSges(2,2) * t186;
t25 = -t258 * t76 - t231;
t226 = t161 * t25 - t24 * t162;
t26 = -t255 * t75 - t284;
t27 = -t255 * t76 + t283;
t225 = t27 * t161 - t26 * t162;
t224 = -t161 * t31 - t273;
t223 = t161 * t79 + t162 * t80;
t45 = t183 * t77 + t185 * t75;
t46 = t183 * t78 + t185 * t76;
t218 = t141 * t185 + t143 * t183;
t83 = t108 + t296;
t210 = -t183 * t279 + t185 * t281;
t61 = t161 * t237 + t157 + t252;
t82 = -t107 + t228;
t54 = t62 + t296;
t209 = (-t183 * t250 + t185 * t251) * t173;
t208 = Icges(6,5) * t173 - qJD(5) * t143;
t207 = Icges(6,6) * t173 - qJD(5) * t141;
t53 = t228 + t61;
t56 = -t162 * t217 + t262;
t52 = t56 * t173;
t10 = qJD(5) * t225 + t52;
t126 = t219 * qJD(5);
t127 = t144 * qJD(5);
t42 = t161 * t207 + t173 * t213;
t44 = t161 * t208 + t173 * t214;
t16 = -qJD(5) * t222 + t183 * t44 + t185 * t42;
t41 = t162 * t207 - t219 * t259;
t43 = -t144 * t259 + t162 * t208;
t17 = -qJD(5) * t221 + t183 * t43 + t185 * t41;
t193 = -qJD(5) * t218 - t126 * t183 + t127 * t185 + t139 * t173;
t20 = t161 * t291 + t193 * t162;
t21 = t193 * t161 - t162 * t291;
t9 = qJD(5) * t226 + t269;
t201 = (t52 + ((t25 - t63 + (t74 + t271) * t162 + t284) * t162 + t283 * t161) * qJD(5)) * t233 + (-qJD(5) * t217 + t126 * t185 + t127 * t183) * t173 + (t46 + t56) * t289 + (t45 + t55) * t288 + (-t269 + ((t162 * t230 + t27 - t283) * t162 + (t161 * t230 + t231 + t26) * t161) * qJD(5) + t9) * t236 + (t17 + t20) * t235 + (Icges(5,3) + t218) * t172 + (t16 + t21 + t10) * t234;
t197 = -qJD(5) * t45 + t173 * t73 - t183 * t42 + t185 * t44;
t196 = -qJD(5) * t46 + t173 * t74 - t183 * t41 + t185 * t43;
t192 = t201 + (Icges(3,3) + Icges(4,3)) * t180;
t71 = t102 * t181 - t242;
t72 = t181 * t299 + t243;
t191 = (t71 * (-t159 - t167) + t72 * t102) * t181;
t104 = t181 * t116;
t101 = t122 * t181 + t243;
t98 = t146 * t162;
t97 = t146 * t161;
t60 = t205 + t264;
t34 = qJD(5) * t223 + qJD(3);
t13 = t79 * t88 - t80 * t89 + qJDD(3) + (t161 * t48 + t162 * t47) * qJD(5);
t6 = t196 * t161 - t162 * t293;
t5 = t197 * t161 - t162 * t292;
t4 = t161 * t293 + t196 * t162;
t3 = t161 * t292 + t197 * t162;
t1 = [Icges(2,3) * qJDD(1) + t192 + (t302 * (t122 + t178) + t307 * (-t121 - t286) + (-t106 - t243 + t101) * t100) * m(3) + (g(1) * t147 - g(2) * t149 + (t147 ^ 2 + t149 ^ 2) * qJDD(1)) * m(2) + (t30 * (t240 - t243) + (t228 * t31 - t296 * t30) * t181 + t190 + t305 * (t178 + t54) + t309 * (t53 - t286) + (-t242 + t30 - t200 + t311) * t31) * m(6) + (t59 * (-t85 - t243) + (t228 * t60 - t296 * t59) * t181 + t304 * (t178 + t83) - t306 * (t82 - t286) + t60 * (-rSges(5,1) * t259 - rSges(5,2) * t256 - t242)) * m(5) + (t71 * (t137 - t243) + t191 - (-pkin(2) * t253 - t104 - t71) * t72 + t303 * (t299 + t178) + t308 * (t102 - t286)) * m(4); t192 + (t72 * t104 - (-t285 * t72 - t299 * t71) * t181 + t71 * t137 + t191 + t303 * t299 + t308 * t102) * m(4) + (-t100 * t106 - t101 * t263 + (t100 * t181 + t302) * t122 + (t101 * t181 - t307) * t121) * m(3) + (t305 * t54 + t309 * t53 + t295) * m(6) + (t304 * t83 - t306 * t82 + (-t155 * t173 + t118 + t264) * t59) * m(5); m(6) * t13 + t287 * qJDD(3) + (-m(6) - t287) * g(3); t201 + (t305 * t62 + t309 * t61 + t295) * m(6) + (-t59 * t85 + (t173 * t59 + t304) * t108 + t306 * t107) * m(5); t10 * t256 / 0.2e1 + t161 * (t56 * t172 + t20 * t173 + t26 * t89 + t27 * t88 + (t161 * t4 - t3 * t162) * qJD(5)) / 0.2e1 + t225 * t289 + ((t173 * t27 - t3) * t162 + (t173 * t26 + t4) * t161) * t235 + t9 * t259 / 0.2e1 - t162 * (t55 * t172 + t21 * t173 + t24 * t89 + t25 * t88 + (t161 * t6 - t162 * t5) * qJD(5)) / 0.2e1 + t226 * t288 + ((t173 * t25 - t5) * t162 + (t173 * t24 + t6) * t161) * t234 + t172 * (t46 * t161 - t45 * t162) / 0.2e1 + t173 * ((t173 * t46 - t16) * t162 + (t173 * t45 + t17) * t161) / 0.2e1 + ((-t248 * t261 + t260) * t161 + (t209 + (-t290 * t162 + (t262 + t210) * t161) * qJD(5)) * t162) * t236 + ((-t247 * t262 - t260) * t162 + (t209 + (t210 * t161 + (-t290 + t261) * t162) * qJD(5)) * t161) * t233 - t173 * ((t251 * t183 + t250 * t185) * t173 + ((t161 * t279 - t162 * t280) * t185 + (t161 * t281 + t162 * t282) * t183) * qJD(5)) / 0.2e1 + (t13 * t223 + t34 * ((t47 + t68) * t162 + (-t173 * t80 + t48) * t161) + t224 * t128 + ((-t173 * t31 - t11) * t162 + (t173 * t30 - t12) * t161) * t146 - (t30 * t97 - t31 * t98) * t173 - (t34 * (-t161 * t97 - t162 * t98) + t224 * t148) * qJD(5) + g(1) * t98 + g(2) * t97 - g(3) * t148) * m(6);];
tau = t1;
