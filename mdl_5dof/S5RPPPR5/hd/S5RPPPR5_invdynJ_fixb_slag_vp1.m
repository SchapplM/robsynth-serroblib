% Calculate vector of inverse dynamics joint torques for
% S5RPPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:28
% DurationCPUTime: 5.77s
% Computational Cost: add. (5440->420), mult. (9322->500), div. (0->0), fcn. (9482->8), ass. (0->215)
t299 = sin(pkin(7));
t300 = cos(pkin(7));
t317 = sin(qJ(1));
t318 = cos(qJ(1));
t153 = -t299 * t317 - t300 * t318;
t154 = t318 * t299 - t317 * t300;
t190 = pkin(8) + qJ(5);
t178 = sin(t190);
t179 = cos(t190);
t297 = Icges(6,4) * t179;
t228 = -Icges(6,2) * t178 + t297;
t61 = Icges(6,6) * t153 + t154 * t228;
t298 = Icges(6,4) * t178;
t230 = Icges(6,1) * t179 - t298;
t64 = Icges(6,5) * t153 + t154 * t230;
t335 = t178 * t61 - t179 * t64;
t226 = Icges(6,5) * t179 - Icges(6,6) * t178;
t58 = Icges(6,3) * t153 + t154 * t226;
t20 = t153 * t58 - t335 * t154;
t60 = Icges(6,3) * t154 - t153 * t226;
t339 = t154 * t60;
t227 = Icges(6,2) * t179 + t298;
t229 = Icges(6,1) * t178 + t297;
t223 = -t178 * t227 + t179 * t229;
t225 = Icges(6,5) * t178 + Icges(6,6) * t179;
t78 = t225 * t153;
t43 = t154 * t223 + t78;
t338 = qJD(1) * t43;
t136 = t153 * rSges(5,3);
t191 = sin(pkin(8));
t192 = cos(pkin(8));
t241 = -rSges(5,1) * t192 + rSges(5,2) * t191;
t337 = -t154 * t241 + t136;
t234 = -t178 * t64 - t179 * t61;
t77 = t225 * t154;
t274 = t318 * pkin(1) + t317 * qJ(2);
t159 = qJD(1) * t274;
t182 = qJD(2) * t318;
t123 = -t182 + t159;
t257 = qJD(1) * t318;
t276 = -pkin(2) * t257 + t182;
t260 = -t159 + t276;
t268 = t153 * qJD(4);
t94 = -t153 * pkin(3) + qJ(4) * t154;
t328 = qJD(1) * t94 - t260 - t268;
t184 = t318 * qJ(2);
t264 = t317 * pkin(1);
t160 = t264 - t184;
t181 = qJD(2) * t317;
t278 = qJ(2) * t257 + t181;
t326 = qJD(1) * t160 - t181 + t278;
t171 = pkin(4) * t192 + pkin(3);
t193 = -pkin(6) - qJ(4);
t288 = t153 * t179;
t289 = t153 * t178;
t69 = -rSges(6,1) * t288 + rSges(6,2) * t289 + t154 * rSges(6,3);
t333 = -t153 * t171 - t154 * t193 + t69;
t332 = t154 * rSges(4,1) - t153 * rSges(4,2);
t134 = t153 * qJ(4);
t250 = t154 * pkin(3) + t134;
t186 = t318 * rSges(3,3);
t262 = t317 * rSges(3,1);
t161 = t262 - t186;
t331 = -t160 - t161;
t263 = t317 * pkin(2);
t210 = -t264 - t263;
t330 = t210 + t263;
t74 = t154 * rSges(5,3) + t241 * t153;
t30 = qJD(1) * t74 + t328;
t275 = t318 * rSges(3,1) + t317 * rSges(3,3);
t327 = t274 + t275;
t132 = qJD(4) * t154;
t272 = qJD(5) * t154;
t126 = t153 * t193;
t135 = t153 * rSges(6,3);
t286 = t154 * t179;
t287 = t154 * t178;
t67 = rSges(6,1) * t286 - rSges(6,2) * t287 + t135;
t325 = t154 * t171 - t126 + t67;
t127 = t153 * qJD(1);
t247 = -t160 - t263;
t222 = t247 + t250;
t316 = pkin(3) - t171;
t152 = -rSges(6,1) * t179 + rSges(6,2) * t178;
t68 = t152 * t154 - t135;
t207 = -t154 * t316 - t126 - t134 + t222 - t68;
t240 = rSges(6,1) * t178 + rSges(6,2) * t179;
t273 = qJD(5) * t153;
t282 = t132 + t181;
t18 = qJD(1) * t207 + t240 * t273 + t282;
t128 = t154 * qJD(1);
t114 = t128 * qJ(4);
t124 = t152 * qJD(5);
t188 = t318 * pkin(2);
t194 = qJD(1) ^ 2;
t279 = qJD(1) * t182 + qJDD(2) * t317;
t212 = -t188 * t194 + t279;
t202 = qJD(4) * t127 + qJDD(4) * t154 + t212;
t291 = t128 * t193;
t301 = t127 * pkin(3) - t114 - t123 + t268;
t271 = qJD(5) * t178;
t215 = -t127 * t179 + t154 * t271;
t270 = qJD(5) * t179;
t216 = t127 * t178 + t154 * t270;
t304 = t128 * rSges(6,3);
t39 = rSges(6,1) * t215 + rSges(6,2) * t216 + t304;
t86 = qJD(5) * t128 - qJDD(5) * t153;
t5 = -t124 * t273 - t86 * t240 + (-t127 * t316 + t114 + t291 + t301 - t39) * qJD(1) + t207 * qJDD(1) + t202;
t324 = t127 * t18 + t154 * t5;
t251 = t128 * pkin(3) + t127 * qJ(4);
t323 = -qJD(1) * t250 + t326;
t310 = t227 * t154 - t64;
t312 = -t229 * t154 - t61;
t322 = t178 * t310 + t179 * t312;
t85 = qJD(5) * t127 + qJDD(5) * t154;
t321 = t85 / 0.2e1;
t320 = t86 / 0.2e1;
t319 = -m(5) - m(6);
t66 = Icges(6,5) * t154 - t153 * t230;
t315 = -t153 * t60 - t66 * t286;
t314 = -t66 * t288 + t339;
t313 = -t94 + t333;
t63 = Icges(6,6) * t154 - t153 * t228;
t311 = -t229 * t153 + t63;
t309 = t227 * t153 + t66;
t305 = t128 * rSges(5,3);
t303 = t178 * t63;
t302 = -t127 * t193 + t128 * t171;
t293 = t128 * t178;
t292 = t128 * t179;
t285 = rSges(6,1) * t292 + t127 * rSges(6,3);
t283 = t128 * rSges(4,1) - t127 * rSges(4,2);
t95 = -t153 * rSges(4,1) - t154 * rSges(4,2);
t281 = -t227 + t230;
t280 = -t228 - t229;
t269 = t226 * qJD(1);
t267 = m(4) - t319;
t266 = -t318 / 0.2e1;
t265 = t317 / 0.2e1;
t261 = t268 + t276;
t259 = t188 + t274;
t258 = t153 * t271;
t256 = qJD(1) * t317;
t255 = -t273 / 0.2e1;
t254 = t273 / 0.2e1;
t253 = -t272 / 0.2e1;
t252 = t272 / 0.2e1;
t232 = t178 * t66 + t179 * t63;
t213 = t258 + t292;
t214 = t153 * t270 - t293;
t36 = Icges(6,4) * t213 + Icges(6,2) * t214 + Icges(6,6) * t127;
t38 = Icges(6,1) * t213 + Icges(6,4) * t214 + Icges(6,5) * t127;
t198 = qJD(5) * t232 + t178 * t36 - t179 * t38;
t35 = Icges(6,4) * t215 + Icges(6,2) * t216 + Icges(6,6) * t128;
t37 = Icges(6,1) * t215 + Icges(6,4) * t216 + Icges(6,5) * t128;
t199 = qJD(5) * t234 + t178 * t35 - t179 * t37;
t231 = -t179 * t66 + t303;
t33 = Icges(6,5) * t215 + Icges(6,6) * t216 + Icges(6,3) * t128;
t34 = Icges(6,5) * t213 + Icges(6,6) * t214 + Icges(6,3) * t127;
t244 = -(-t127 * t335 - t128 * t58 - t153 * t33 + t154 * t199) * t153 + t154 * (t127 * t231 + t128 * t60 - t153 * t34 + t154 * t198);
t243 = -t153 * (-t127 * t58 + t128 * t335 + t153 * t199 + t154 * t33) + t154 * (t127 * t60 - t128 * t231 + t153 * t198 + t154 * t34);
t242 = t127 * rSges(4,1) + t128 * rSges(4,2);
t19 = t313 * qJD(1) + t240 * t272 + t328;
t239 = t153 * t18 + t154 * t19;
t21 = t287 * t63 + t315;
t238 = -t153 * t20 + t154 * t21;
t22 = -t154 * t58 + t288 * t64 - t289 * t61;
t23 = t289 * t63 + t314;
t237 = -t153 * t22 + t154 * t23;
t40 = rSges(6,1) * t258 + rSges(6,2) * t214 + t285;
t236 = -t153 * t40 - t154 * t39;
t235 = t153 * t69 + t154 * t68;
t224 = -t178 * t229 - t179 * t227;
t221 = t332 + t247;
t220 = t127 * rSges(5,3) - t241 * t128;
t219 = pkin(3) - t241;
t217 = t222 + t337;
t209 = -t262 - t264;
t208 = qJDD(1) * t274 - qJDD(2) * t318 + (-pkin(1) * t256 + t181 + t278) * qJD(1);
t165 = rSges(2,1) * t318 - rSges(2,2) * t317;
t162 = rSges(2,1) * t317 + rSges(2,2) * t318;
t204 = t178 * t309 + t179 * t311;
t203 = t184 + t210;
t201 = -t263 + t337;
t200 = (t178 * t280 + t179 * t281) * qJD(1);
t197 = qJDD(1) * t188 - t194 * t263 + t208;
t121 = t228 * qJD(5);
t122 = t230 * qJD(5);
t196 = qJD(5) * t224 - t121 * t178 + t122 * t179;
t195 = qJD(4) * t128 + qJD(1) * (t251 + t132) + qJDD(1) * t94 - qJDD(4) * t153 + t197;
t176 = rSges(3,3) * t257;
t120 = t226 * qJD(5);
t97 = qJD(1) * t331 + t181;
t84 = t240 * t153;
t83 = t240 * t154;
t70 = qJD(1) * t221 + t181;
t53 = qJDD(1) * t275 + qJD(1) * (-rSges(3,1) * t256 + t176) + t208;
t52 = -qJD(1) * t123 + qJDD(1) * t331 - t194 * t275 + t279;
t44 = t153 * t223 - t77;
t42 = t44 * qJD(1);
t29 = qJD(1) * t217 + t282;
t28 = qJD(1) * t283 + qJDD(1) * t95 + t197;
t27 = (-t123 + t242) * qJD(1) + t221 * qJDD(1) + t212;
t24 = qJD(5) * t235 - qJD(3);
t15 = -t120 * t154 - t127 * t225 - t128 * t223 + t153 * t196;
t14 = t120 * t153 + t127 * t223 - t128 * t225 + t154 * t196;
t13 = qJD(1) * t220 + qJDD(1) * t74 + t195;
t12 = (-t127 * t241 + t301 - t305) * qJD(1) + t217 * qJDD(1) + t202;
t11 = qJD(5) * t231 - t178 * t38 - t179 * t36;
t10 = -qJD(5) * t335 - t178 * t37 - t179 * t35;
t9 = qJD(5) * t236 - t68 * t85 + t69 * t86 + qJDD(3);
t8 = qJD(5) * t237 + t42;
t7 = qJD(5) * t238 + t338;
t6 = t313 * qJDD(1) + (-t251 + t40 + t302) * qJD(1) + t195 - t124 * t272 + t85 * t240;
t1 = [(qJD(5) * t223 + t121 * t179 + t122 * t178) * qJD(1) + (t42 + ((t21 - t22 - t315) * t153 + t314 * t154) * qJD(5)) * t254 - m(2) * (-g(1) * t162 + g(2) * t165) + (t44 - t232) * t321 + (t43 - t234) * t320 + (t11 + t15) * t252 + (t8 + t10 + t14) * t255 + (t7 + ((t22 + (-t231 + t58) * t154) * t154 + (t23 + (t58 - t303) * t153 + t339 - t314) * t153) * qJD(5) - t338) * t253 + (t5 * t203 - g(1) * (t203 + t325) + (t6 - g(2)) * (t259 + t333) + t324 * (-t152 + t171) + (-rSges(6,2) * t293 + t285 + t302 + t323) * t19 + (t5 * (rSges(6,3) - t193) - t24 * (t67 + t68) * qJD(5)) * t153 + (t261 + t291 - t304 + t328) * t18 + ((t250 - t325 + t330) * t19 + (-t274 + t313) * t18) * qJD(1)) * m(6) + (-g(1) * (-t160 + t201 + t250) + (t13 - g(2)) * (t74 + t94 + t259) + (t219 * t154 + t134 + t136 + t203) * t12 + (t220 + (-t201 + t210) * qJD(1) + t323 + t251) * t30 + (t219 * t127 - t114 - t159 + t261 + t30 - t305) * t29) * m(5) + ((t28 - g(2)) * (t259 + t95) + (t242 + t260) * t70 + (t70 + t283 + (t330 - t332) * qJD(1) + t326) * (qJD(1) * t95 - t260) + (t27 - g(1)) * (t203 + t332)) * m(4) + ((-qJD(1) * t327 + t182) * t97 + (t176 + t97 + (t161 + t209) * qJD(1) + t326) * (qJD(1) * t275 + t123) + (t53 - g(2)) * t327 + (t52 - g(1)) * (t184 + t186 + t209)) * m(3) + (-t224 + m(2) * (t162 ^ 2 + t165 ^ 2) + t192 ^ 2 * Icges(5,2) + (Icges(5,1) * t191 + 0.2e1 * Icges(5,4) * t192) * t191 + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJDD(1); (-m(3) - t267) * (g(1) * t317 - g(2) * t318) + 0.2e1 * (t265 * t5 + t266 * t6) * m(6) + 0.2e1 * (t12 * t265 + t13 * t266) * m(5) + 0.2e1 * (t265 * t27 + t266 * t28) * m(4) + 0.2e1 * (t265 * t52 + t266 * t53) * m(3); m(6) * t9 + (m(4) + m(5)) * qJDD(3) + t267 * g(3); t319 * (g(1) * t154 - g(2) * t153) + (t12 * t154 + t127 * t29 + t128 * t30 - t13 * t153 - (t153 * t29 + t154 * t30) * qJD(1)) * m(5) + (-qJD(1) * t239 + t128 * t19 - t153 * t6 + t324) * m(6); t127 * t8 / 0.2e1 + t154 * (qJD(1) * t15 + qJD(5) * t243 + qJDD(1) * t44 + t22 * t86 + t23 * t85) / 0.2e1 + t237 * t321 + (t127 * t23 + t128 * t22 + t243) * t252 + t128 * t7 / 0.2e1 - t153 * (qJD(1) * t14 + qJD(5) * t244 + qJDD(1) * t43 + t20 * t86 + t21 * t85) / 0.2e1 + t238 * t320 + (t127 * t21 + t128 * t20 + t244) * t255 + qJDD(1) * (t153 * t234 - t154 * t232) / 0.2e1 + qJD(1) * (-t10 * t153 + t11 * t154 - t127 * t232 - t128 * t234) / 0.2e1 + ((t272 * t78 - t269) * t154 + (t200 + (-t322 * t153 + (-t77 + t204) * t154) * qJD(5)) * t153) * t253 + ((t273 * t77 + t269) * t153 + (t200 + (t204 * t154 + (-t322 - t78) * t153) * qJD(5)) * t154) * t254 - qJD(1) * ((t281 * t178 - t280 * t179) * qJD(1) + ((t153 * t310 - t154 * t309) * t179 + (-t153 * t312 + t154 * t311) * t178) * qJD(5)) / 0.2e1 + (-t9 * t235 + t24 * (t127 * t68 - t128 * t69 - t236) - t239 * t124 - (-t127 * t19 + t128 * t18 - t153 * t5 - t154 * t6) * t240 - (-t18 * t83 + t19 * t84) * qJD(1) - (t24 * (t153 * t84 + t154 * t83) - t239 * t152) * qJD(5) - g(1) * t84 - g(2) * t83 - g(3) * t152) * m(6);];
tau = t1;
