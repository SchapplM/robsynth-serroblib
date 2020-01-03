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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:27
% DurationCPUTime: 4.00s
% Computational Cost: add. (11057->378), mult. (6222->471), div. (0->0), fcn. (4592->10), ass. (0->229)
t194 = qJD(1) + qJD(2);
t195 = qJ(1) + qJ(2);
t185 = pkin(9) + t195;
t172 = sin(t185);
t186 = sin(t195);
t174 = pkin(2) * t186;
t324 = pkin(3) * t172 + t174;
t226 = t324 * t194;
t197 = sin(qJ(1));
t301 = pkin(1) * qJD(1);
t259 = t197 * t301;
t210 = t259 + t226;
t184 = qJD(4) + t194;
t176 = qJ(4) + t185;
t166 = sin(t176);
t167 = cos(t176);
t325 = t166 * rSges(5,1) + t167 * rSges(5,2);
t340 = t325 * t184;
t59 = t210 + t340;
t173 = cos(t185);
t114 = t172 * rSges(4,1) + t173 * rSges(4,2);
t322 = t174 + t114;
t187 = cos(t195);
t323 = t186 * rSges(3,1) + t187 * rSges(3,2);
t104 = t323 * t194;
t98 = t259 + t104;
t196 = sin(qJ(5));
t198 = cos(qJ(5));
t188 = Icges(6,4) * t198;
t232 = -Icges(6,2) * t196 + t188;
t321 = Icges(6,1) * t196 + t188;
t271 = t321 + t232;
t294 = Icges(6,4) * t196;
t144 = Icges(6,2) * t198 + t294;
t147 = Icges(6,1) * t198 - t294;
t272 = t144 - t147;
t341 = (t196 * t271 + t198 * t272) * t184;
t109 = t167 * pkin(4) + t166 * pkin(8);
t282 = t167 * t196;
t257 = t184 * t282;
t283 = t167 * t184;
t286 = t166 * t184;
t274 = pkin(4) * t283 + pkin(8) * t286;
t281 = t167 * t198;
t261 = rSges(6,1) * t281;
t275 = rSges(6,3) * t286 + t184 * t261;
t78 = -rSges(6,2) * t282 + rSges(6,3) * t166 + t261;
t66 = t184 * t78;
t339 = -rSges(6,2) * t257 - t184 * t109 + t274 + t275 - t66;
t151 = rSges(6,1) * t196 + rSges(6,2) * t198;
t265 = qJD(5) * t167;
t159 = t166 * pkin(4);
t285 = t166 * t196;
t260 = rSges(6,2) * t285;
t284 = t166 * t198;
t245 = rSges(6,1) * t284 - t260;
t77 = -rSges(6,3) * t167 + t245;
t296 = -pkin(8) * t167 + t159 + t77;
t216 = -t151 * t265 - t184 * t296;
t30 = t210 - t216;
t199 = cos(qJ(1));
t182 = t199 * t301;
t266 = qJD(5) * t166;
t255 = t151 * t266;
t279 = t173 * t194;
t139 = pkin(3) * t279;
t277 = t187 * t194;
t150 = pkin(2) * t277;
t273 = t139 + t150;
t228 = -t255 + t273;
t62 = t109 + t78;
t31 = t184 * t62 + t182 + t228;
t95 = t151 * t166;
t96 = t151 * t167;
t338 = -t30 * t95 - t31 * t96;
t153 = rSges(6,1) * t198 - rSges(6,2) * t196;
t127 = t153 * qJD(5);
t193 = qJDD(1) + qJDD(2);
t183 = qJDD(4) + t193;
t165 = pkin(3) * t173;
t192 = t194 ^ 2;
t191 = t199 * pkin(1);
t200 = qJD(1) ^ 2;
t289 = pkin(1) * qJDD(1);
t267 = t200 * t191 + t197 * t289;
t308 = pkin(2) * t193;
t309 = pkin(2) * t192;
t248 = t186 * t308 + t187 * t309 + t267;
t307 = pkin(3) * t193;
t227 = t192 * t165 + t172 * t307 + t248;
t262 = qJD(5) * t198;
t263 = qJD(5) * t196;
t48 = -rSges(6,1) * t166 * t263 + (-t166 * t262 - t257) * rSges(6,2) + t275;
t264 = qJD(5) * t184;
t87 = -qJDD(5) * t167 + t166 * t264;
t11 = t127 * t265 - t87 * t151 + (t48 + t274) * t184 + t296 * t183 + t227;
t337 = t11 - g(3);
t190 = t197 * pkin(1);
t247 = -t190 * t200 + t199 * t289;
t229 = t187 * t308 + t247;
t207 = t173 * t307 - t192 * t324 + t229;
t246 = -pkin(4) * t286 + pkin(8) * t283;
t258 = t184 * t284;
t276 = -rSges(6,3) * t283 - t184 * t260;
t47 = rSges(6,2) * t167 * t262 + (t167 * t263 + t258) * rSges(6,1) + t276;
t86 = -qJDD(5) * t166 - t167 * t264;
t12 = -t127 * t266 + t86 * t151 + (t246 - t47) * t184 + t62 * t183 + t207;
t336 = t12 - g(2);
t83 = rSges(5,1) * t283 - rSges(5,2) * t286;
t335 = t183 * t325 + t184 * t83 - g(3) + t227;
t107 = rSges(5,1) * t167 - t166 * rSges(5,2);
t334 = t107 * t183 - t184 * t340 - g(2) + t207;
t137 = rSges(4,1) * t279;
t280 = t172 * t194;
t333 = t193 * t114 + t194 * (-rSges(4,2) * t280 + t137) + t248 - g(3);
t161 = t172 * rSges(4,2);
t115 = rSges(4,1) * t173 - t161;
t332 = -t114 * t192 + t193 * t115 - t186 * t309 - g(2) + t229;
t278 = t186 * t194;
t105 = rSges(3,1) * t277 - rSges(3,2) * t278;
t331 = t105 * t194 + t193 * t323 - g(3) + t267;
t123 = rSges(3,1) * t187 - t186 * rSges(3,2);
t330 = -t104 * t194 + t123 * t193 - g(2) + t247;
t270 = t150 + t182;
t94 = t184 * t107;
t60 = t139 + t270 + t94;
t327 = t194 * t322;
t102 = t194 * t115;
t326 = t137 - t102;
t320 = -t228 + t273 + t339;
t223 = t232 * t184;
t317 = -Icges(6,6) * t184 + qJD(5) * t144;
t42 = -t166 * t317 + t167 * t223;
t224 = t147 * t184;
t314 = -Icges(6,5) * t184 + qJD(5) * t321;
t44 = -t166 * t314 + t167 * t224;
t73 = -Icges(6,6) * t167 + t166 * t232;
t75 = -Icges(6,5) * t167 + t147 * t166;
t45 = t196 * t75 + t198 * t73;
t143 = Icges(6,5) * t198 - Icges(6,6) * t196;
t71 = -Icges(6,3) * t167 + t143 * t166;
t319 = qJD(5) * t45 - t184 * t71 + t196 * t42 - t198 * t44;
t142 = Icges(6,5) * t196 + Icges(6,6) * t198;
t318 = -Icges(6,3) * t184 + qJD(5) * t142;
t125 = t232 * qJD(5);
t126 = t147 * qJD(5);
t231 = t144 * t198 + t196 * t321;
t316 = qJD(5) * t231 + t125 * t196 - t126 * t198 - t142 * t184;
t74 = Icges(6,4) * t281 - Icges(6,2) * t282 + Icges(6,6) * t166;
t132 = Icges(6,4) * t282;
t76 = Icges(6,1) * t281 + Icges(6,5) * t166 - t132;
t235 = t196 * t76 + t198 * t74;
t41 = t166 * t223 + t167 * t317;
t43 = t166 * t224 + t167 * t314;
t72 = Icges(6,5) * t281 - Icges(6,6) * t282 + Icges(6,3) * t166;
t315 = qJD(5) * t235 - t184 * t72 - t196 * t41 + t198 * t43;
t312 = t86 / 0.2e1;
t311 = t87 / 0.2e1;
t310 = m(4) + m(5);
t305 = t166 * t321 + t73;
t304 = t167 * t321 + t74;
t303 = -t144 * t166 + t75;
t302 = -Icges(6,2) * t281 - t132 + t76;
t300 = t196 * t73;
t299 = t196 * t74;
t298 = t198 * t75;
t287 = t142 * t166;
t56 = -t144 * t282 + t281 * t321 + t287;
t297 = t56 * t184;
t89 = t142 * t167;
t222 = t143 * t184;
t175 = pkin(2) * t187;
t268 = t165 + t175;
t253 = -t266 / 0.2e1;
t252 = t266 / 0.2e1;
t251 = -t265 / 0.2e1;
t250 = t265 / 0.2e1;
t99 = t123 * t194 + t182;
t80 = t324 + t325;
t101 = t115 + t175;
t154 = rSges(2,1) * t199 - t197 * rSges(2,2);
t152 = rSges(2,1) * t197 + rSges(2,2) * t199;
t63 = t75 * t284;
t24 = -t167 * t71 - t285 * t73 + t63;
t64 = t76 * t284;
t25 = t167 * t72 + t285 * t74 - t64;
t240 = -t166 * t25 - t167 * t24;
t65 = t73 * t282;
t26 = -t166 * t71 - t281 * t75 + t65;
t234 = -t198 * t76 + t299;
t27 = t166 * t72 - t234 * t167;
t239 = -t166 * t27 - t167 * t26;
t238 = -t166 * t31 + t167 * t30;
t237 = t166 * t77 + t167 * t78;
t236 = t298 - t300;
t230 = -t144 * t196 + t198 * t321;
t81 = t107 + t268;
t225 = t83 + t273;
t218 = t196 * t303 + t198 * t305;
t217 = t196 * t302 + t198 * t304;
t214 = -t166 * t222 - t167 * t318 + t184 * t234;
t213 = t166 * t318 - t167 * t222 + t184 * t236;
t212 = -t143 * qJD(5) + t184 * t230;
t61 = t159 + (-rSges(6,3) - pkin(8)) * t167 + t245;
t208 = -rSges(6,1) * t258 + t246 - t276;
t54 = t62 + t268;
t53 = t61 + t324;
t10 = qJD(5) * t239 - t297;
t16 = qJD(5) * t236 + t196 * t44 + t198 * t42;
t17 = qJD(5) * t234 + t196 * t43 + t198 * t41;
t20 = t212 * t166 + t167 * t316;
t21 = -t166 * t316 + t212 * t167;
t55 = t166 * t230 - t89;
t52 = t55 * t184;
t9 = qJD(5) * t240 + t52;
t206 = (t52 + ((t26 + t64 - t65 + (t71 - t299) * t166) * t166 + (-t63 - t27 + (-t234 + t71) * t167 + (t298 + t300) * t166) * t167) * qJD(5)) * t252 - t235 * t312 - t86 * t56 / 0.2e1 + (qJD(5) * t230 + t125 * t198 + t126 * t196) * t184 + (t45 + t55) * t311 + (t297 + ((-t25 + t65 + (t72 - t298) * t167) * t167 + (t24 - t63 + (t72 + t300) * t166) * t166) * qJD(5) + t10) * t250 + (Icges(5,3) + t231) * t183 + (t17 + t20 + t9) * t253 + (t16 + t21) * t251;
t69 = t259 + t327;
t70 = t270 + t102;
t204 = (-t69 * t161 - t70 * t322) * t194;
t203 = t206 + (Icges(3,3) + Icges(4,3)) * t193;
t202 = -pkin(2) * t278 - pkin(3) * t280 + t208;
t201 = t338 * qJD(5);
t34 = qJD(5) * t237 + qJD(3);
t13 = -t77 * t86 - t78 * t87 + qJDD(3) + (t166 * t48 - t167 * t47) * qJD(5);
t6 = t166 * t315 + t214 * t167;
t5 = -t166 * t319 + t213 * t167;
t4 = t214 * t166 - t167 * t315;
t3 = t213 * t166 + t167 * t319;
t1 = [Icges(2,3) * qJDD(1) + t203 + (t330 * (t123 + t191) + t331 * (t190 + t323) + (-t99 + t105 + t182) * t98) * m(3) + ((t152 ^ 2 + t154 ^ 2) * qJDD(1) - g(2) * t154 - g(3) * t152) * m(2) + (t31 * (t202 - t259) + t201 + t336 * (t191 + t54) + t337 * (t190 + t53) + (t31 + t320) * t30) * m(6) + (t334 * (t191 + t81) + t335 * (t190 + t80) + (-t60 + t182 + t225) * t59) * m(5) + (-t70 * t259 + t204 + t332 * (t101 + t191) + t333 * (t190 + t322) + (t70 + t326) * t69) * m(4); t203 + (t201 + t336 * t54 + t337 * t53 + (t226 - t216 + t202) * t31 + t320 * t30) * m(6) + (t334 * t81 + t335 * t80 + (t225 - t94 - t273) * t59) * m(5) + (t332 * t101 + t322 * t333 + t326 * t69 + t70 * t327 + t204) * m(4) + (-t104 * t99 + t105 * t98 + (-t194 * t98 + t330) * t123 + (t194 * t99 + t331) * t323) * m(3); m(6) * t13 + t310 * qJDD(3) + (-m(6) - t310) * g(1); t206 + (t201 + t336 * t62 + t337 * t61 + (-t216 + t208) * t31 + (t255 + t339) * t30) * m(6) + (t59 * t83 - t60 * t340 + (-t184 * t59 + t334) * t107 + (t184 * t60 + t335) * t325) * m(5); t183 * (t166 * t235 - t167 * t45) / 0.2e1 + t184 * ((t184 * t235 - t16) * t167 + (t184 * t45 - t17) * t166) / 0.2e1 + t9 * t286 / 0.2e1 - t167 * (t55 * t183 + t21 * t184 + t24 * t87 + t25 * t86 + (-t166 * t6 - t167 * t5) * qJD(5)) / 0.2e1 + t240 * t311 + ((-t184 * t25 - t5) * t167 + (t184 * t24 - t6) * t166) * t251 - t10 * t283 / 0.2e1 - t166 * (-t56 * t183 + t20 * t184 + t26 * t87 + t27 * t86 + (-t166 * t4 - t167 * t3) * qJD(5)) / 0.2e1 + t239 * t312 + ((-t184 * t27 - t3) * t167 + (t184 * t26 - t4) * t166) * t253 - t184 * ((-t272 * t196 + t271 * t198) * t184 + ((t166 * t302 - t167 * t303) * t198 + (-t166 * t304 + t167 * t305) * t196) * qJD(5)) / 0.2e1 + ((-t265 * t287 - t222) * t167 + (-t341 + (-t217 * t166 + (t89 + t218) * t167) * qJD(5)) * t166) * t250 + ((t89 * t266 - t222) * t166 + (t341 + (-t218 * t167 + (-t287 + t217) * t166) * qJD(5)) * t167) * t252 + (t13 * t237 + t34 * ((t184 * t77 - t47) * t167 + (t48 - t66) * t166) + t238 * t127 + ((-t184 * t31 + t11) * t167 + (-t184 * t30 - t12) * t166) * t151 - t338 * t184 - (t34 * (-t166 * t95 - t167 * t96) + t238 * t153) * qJD(5) - g(1) * t153 + g(2) * t95 - g(3) * t96) * m(6);];
tau = t1;
