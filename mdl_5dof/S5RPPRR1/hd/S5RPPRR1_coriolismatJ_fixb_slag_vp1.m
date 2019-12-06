% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:00
% EndTime: 2019-12-05 17:38:08
% DurationCPUTime: 3.39s
% Computational Cost: add. (8525->301), mult. (12122->438), div. (0->0), fcn. (10840->6), ass. (0->202)
t235 = sin(qJ(4));
t237 = cos(qJ(4));
t215 = rSges(5,1) * t237 - rSges(5,2) * t235;
t236 = sin(qJ(1));
t232 = t236 ^ 2;
t238 = cos(qJ(1));
t233 = t238 ^ 2;
t220 = t232 + t233;
t134 = t220 * t215;
t234 = qJ(4) + qJ(5);
t226 = cos(t234);
t313 = rSges(6,1) * t226;
t225 = sin(t234);
t368 = rSges(6,2) * t225;
t201 = t313 - t368;
t318 = pkin(4) * t237;
t146 = (t201 + t318) * t236;
t283 = t237 * t238;
t224 = pkin(4) * t283;
t148 = t201 * t238 + t224;
t350 = m(6) / 0.2e1;
t366 = m(5) / 0.2e1;
t305 = t134 * t366 + (t146 * t236 + t148 * t238) * t350;
t284 = t236 * t237;
t286 = t235 * t236;
t193 = rSges(5,1) * t284 - rSges(5,2) * t286;
t194 = t215 * t238;
t117 = -t193 * t236 - t194 * t238;
t218 = t236 * t368;
t143 = t218 + (-t313 - t318) * t236;
t287 = t226 * t238;
t291 = t225 * t238;
t183 = rSges(6,1) * t287 - rSges(6,2) * t291;
t144 = t183 + t224;
t307 = t117 * t366 + (t143 * t236 - t144 * t238) * t350;
t24 = t307 - t305;
t369 = t24 * qJD(1);
t301 = Icges(6,4) * t226;
t197 = -Icges(6,2) * t225 + t301;
t256 = Icges(6,1) * t225 + t301;
t367 = t197 + t256;
t326 = -t236 / 0.2e1;
t364 = -t238 / 0.2e1;
t324 = t238 / 0.2e1;
t258 = rSges(6,1) * t225 + rSges(6,2) * t226;
t121 = t220 * t258;
t250 = Icges(6,5) * t225 + Icges(6,6) * t226;
t363 = t250 * t236;
t362 = t250 * t238;
t296 = t258 * t236;
t295 = t258 * t238;
t303 = Icges(5,4) * t237;
t210 = -Icges(5,2) * t235 + t303;
t257 = Icges(5,1) * t235 + t303;
t361 = (t257 / 0.2e1 + t210 / 0.2e1) * t237;
t182 = t236 * t313 - t218;
t106 = -t182 * t236 - t183 * t238;
t251 = Icges(6,5) * t226 - Icges(6,6) * t225;
t176 = t251 * t236;
t177 = t238 * t251;
t217 = Icges(6,4) * t287;
t154 = Icges(6,1) * t291 - Icges(6,5) * t236 + t217;
t273 = -Icges(6,2) * t291 + t154 + t217;
t153 = t238 * Icges(6,5) + t256 * t236;
t274 = t197 * t236 + t153;
t302 = Icges(6,4) * t225;
t254 = Icges(6,2) * t226 + t302;
t152 = -Icges(6,6) * t236 + t254 * t238;
t199 = Icges(6,1) * t226 - t302;
t275 = -t199 * t238 + t152;
t151 = t238 * Icges(6,6) + t254 * t236;
t276 = t199 * t236 - t151;
t355 = (-t273 * t236 + t274 * t238) * t226 + (t275 * t236 + t276 * t238) * t225;
t316 = (t232 * t177 + (-t236 * t176 + t355) * t238) * t326 + (t233 * t176 + (-t238 * t177 + t355) * t236) * t324;
t155 = t238 * rSges(6,3) + t296;
t228 = t236 * rSges(6,3);
t156 = -t228 + t295;
t223 = pkin(4) * t286;
t285 = t235 * t238;
t57 = (-t155 - t223) * t236 + (-pkin(4) * t285 - t156) * t238;
t6 = t316 + m(6) * (t57 * t106 - t146 * t296 - t148 * t295);
t359 = t6 * qJD(5);
t227 = t238 * qJ(2);
t259 = rSges(5,1) * t235 + rSges(5,2) * t237;
t317 = pkin(1) + qJ(3);
t245 = t259 + t317;
t119 = t227 + (-rSges(5,3) - pkin(6)) * t238 - t245 * t236;
t229 = t236 * rSges(5,3);
t120 = t245 * t238 - t229 + (-pkin(6) + qJ(2)) * t236;
t358 = t119 * t238 + t236 * t120;
t118 = -t193 * t238 + t194 * t236;
t88 = t143 * t238 + t144 * t236;
t306 = t118 * t366 + t88 * t350;
t262 = -t367 * t226 / 0.2e1 + (-t199 / 0.2e1 + t254 / 0.2e1) * t225;
t149 = Icges(6,3) * t238 + t363;
t150 = -Icges(6,3) * t236 + t362;
t248 = t151 * t226 + t153 * t225;
t279 = -t152 * t287 - t154 * t291;
t280 = t151 * t287 + t153 * t291;
t299 = t149 * t236;
t59 = t149 * t238 + t248 * t236;
t61 = t280 - t299;
t5 = (t280 * t238 + (-t59 + (t150 + t248) * t236 + t279) * t236) * t324 + (-t236 * (-t150 * t236 - t279) + t238 * t61) * t364 + ((-t61 - t299) * t236 - t149 * t233 + t238 * t59) * t326;
t221 = Icges(5,4) * t283;
t174 = Icges(5,1) * t285 - Icges(5,5) * t236 + t221;
t269 = -Icges(5,2) * t285 + t174 + t221;
t173 = t238 * Icges(5,5) + t257 * t236;
t270 = t210 * t236 + t173;
t304 = Icges(5,4) * t235;
t255 = Icges(5,2) * t237 + t304;
t172 = -Icges(5,6) * t236 + t255 * t238;
t212 = Icges(5,1) * t237 - t304;
t271 = -t212 * t238 + t172;
t171 = t238 * Icges(5,6) + t255 * t236;
t272 = t212 * t236 - t171;
t356 = (-t269 * t236 + t270 * t238) * t237 + (t271 * t236 + t272 * t238) * t235;
t354 = 0.2e1 * t220;
t353 = 0.4e1 * qJD(1);
t352 = 2 * qJD(4);
t351 = -m(6) / 0.2e1;
t268 = rSges(4,3) + t317;
t158 = rSges(4,2) * t238 - t268 * t236 + t227;
t159 = (rSges(4,2) + qJ(2)) * t236 + t268 * t238;
t349 = m(4) * (-t158 * t236 + t238 * t159);
t348 = m(4) * (t158 * t238 + t236 * t159);
t347 = m(5) * (-t119 * t193 + t120 * t194);
t346 = m(5) * (-t119 * t236 + t238 * t120);
t345 = m(5) * t358;
t239 = -pkin(7) - pkin(6);
t244 = -t258 - t317;
t101 = -t223 + t227 + (-rSges(6,3) + t239) * t238 + t244 * t236;
t319 = pkin(4) * t235;
t102 = -t228 + (qJ(2) + t239) * t236 + (-t244 + t319) * t238;
t97 = t236 * t102;
t315 = -t101 * t295 - t258 * t97;
t341 = m(6) * (t146 * t183 - t148 * t182 + t315);
t340 = m(6) * (t88 * t201 + t315);
t339 = m(6) * (t101 * t143 + t102 * t144);
t338 = m(6) * (-t101 * t182 + t102 * t183);
t337 = m(6) * (-t101 * t236 + t238 * t102);
t336 = m(6) * (t101 * t238 + t97);
t331 = m(6) * (t146 * t238 - t148 * t236);
t107 = -t182 * t238 + t183 * t236;
t329 = t107 / 0.2e1;
t328 = t220 / 0.2e1;
t325 = t236 / 0.2e1;
t323 = m(3) * ((rSges(3,3) * t238 + t227) * t238 + (rSges(3,3) + qJ(2)) * t232);
t314 = m(6) * qJD(5);
t252 = Icges(5,5) * t235 + Icges(5,6) * t237;
t169 = Icges(5,3) * t238 + t252 * t236;
t298 = t169 * t236;
t241 = t106 * t350;
t261 = m(6) * t220 * t201;
t64 = t241 - t261 / 0.2e1;
t282 = t64 * qJD(1);
t72 = (t351 - m(4) / 0.2e1 - m(5) / 0.2e1) * t354;
t281 = t72 * qJD(1);
t278 = t171 * t283 + t173 * t285;
t277 = -t172 * t283 - t174 * t285;
t267 = t107 * t351;
t170 = -Icges(5,3) * t236 + t252 * t238;
t68 = t238 * t170 + t172 * t284 + t174 * t286;
t253 = Icges(5,5) * t237 - Icges(5,6) * t235;
t145 = -t223 - t296;
t147 = (-t258 - t319) * t238;
t249 = t145 * t236 + t147 * t238;
t246 = t171 * t237 + t173 * t235;
t240 = (-t254 + t199) * t226 - t367 * t225;
t243 = -t5 + (-t273 * t225 - t275 * t226 + t240 * t238 + t363) * t326 + (-t274 * t225 + t276 * t226 + t240 * t236 - t362) * t324;
t242 = -t262 + (t325 + t326) * (t151 * t225 - t153 * t226);
t188 = t238 * t253;
t187 = t253 * t236;
t113 = t121 * t314;
t99 = t314 * t329;
t93 = -t236 * t155 - t238 * t156;
t86 = -t220 * t318 + t106;
t84 = t331 / 0.2e1;
t71 = (-m(6) / 0.4e1 - m(4) / 0.4e1 - m(5) / 0.4e1) * t354 + (m(4) + m(5) + m(6)) * t328;
t69 = t278 - t298;
t67 = t169 * t238 + t246 * t236;
t63 = t241 + t261 / 0.2e1;
t45 = -t236 * (-t170 * t236 - t277) + t238 * t69;
t44 = -t236 * t68 + t238 * t67;
t34 = t262 + t338;
t33 = t340 / 0.2e1;
t31 = t341 / 0.2e1;
t30 = t84 - t306;
t29 = t84 + t306;
t28 = -t331 / 0.2e1 + t306;
t27 = t337 + t346 + t349;
t25 = t305 + t307;
t20 = t323 + t336 + t345 + t348;
t17 = -t361 + (-t212 / 0.2e1 + t255 / 0.2e1) * t235 + t347 + t339 + t262;
t16 = t278 * t238 + (-t67 + (t170 + t246) * t236 + t277) * t236;
t15 = (t68 - t69 - t298) * t236 - t169 * t233;
t8 = m(6) * (t93 * t106 - t201 * t121) + t316;
t7 = t8 * qJD(5);
t4 = t31 - t340 / 0.2e1 + t5;
t3 = t33 - t341 / 0.2e1 + t5;
t2 = t31 + t33 + t243;
t1 = (-t45 / 0.2e1 + t16 / 0.2e1) * t238 + (-t15 / 0.2e1 - t44 / 0.2e1) * t236 + t5;
t9 = [qJD(2) * t20 + qJD(3) * t27 + qJD(4) * t17 + qJD(5) * t34, qJD(1) * t20 + qJD(3) * t71 + qJD(4) * t25 + qJD(5) * t63, qJD(1) * t27 + qJD(2) * t71 + qJD(4) * t29 + t99, t17 * qJD(1) + t25 * qJD(2) + t29 * qJD(3) + t2 * qJD(5) + ((t101 * t147 + t102 * t145 + t143 * t148 + t144 * t146) * t350 + (t118 * t215 - t259 * t358) * t366) * t352 + (t16 * t364 + (-t269 * t235 - t271 * t237) * t326 + t243 - (t232 / 0.2e1 + t233 / 0.2e1) * t252 + (t15 + t44) * t325 + (-t270 * t235 + t272 * t237 + t45) * t324) * qJD(4), t34 * qJD(1) + t63 * qJD(2) + t2 * qJD(4) + t243 * qJD(5) + (qJD(3) * t329 + (t107 * t201 + t315) * qJD(5)) * m(6); t72 * qJD(3) + t24 * qJD(4) + t64 * qJD(5) + (-t336 / 0.4e1 - t323 / 0.4e1 - t348 / 0.4e1 - t345 / 0.4e1) * t353, 0, t281, t369 + m(6) * (-t145 * t238 + t147 * t236) * qJD(4), t282; -t72 * qJD(2) + t28 * qJD(4) + t99 + (-t337 / 0.4e1 - t349 / 0.4e1 - t346 / 0.4e1) * t353, -t281, 0, t28 * qJD(1) + (-m(5) * t259 * t328 + t249 * t350) * t352 - t113, -t113 + 0.2e1 * (t107 * qJD(1) / 0.4e1 - t121 * qJD(4) / 0.2e1) * m(6); (t242 + (-t255 + t212) * t235 / 0.2e1 + t361) * qJD(1) - t24 * qJD(2) + t30 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t347 / 0.4e1 - t339 / 0.4e1) * t353, -t369, t30 * qJD(1), t1 * qJD(1) + (m(5) * ((-t236 * (rSges(5,3) * t238 + t259 * t236) + (-t259 * t238 + t229) * t238) * t117 - t259 * t134) + (t232 * t188 + (-t236 * t187 + t356) * t238) * t326 + (t233 * t187 + (-t238 * t188 + t356) * t236) * t324 + m(6) * (t145 * t146 + t147 * t148 + t57 * t86) + t316) * qJD(4) + t359, t4 * qJD(1) + t6 * qJD(4) + t359; (t242 - t338) * qJD(1) - t64 * qJD(2) + qJD(3) * t267 + t3 * qJD(4) + t5 * qJD(5), -t282, qJD(1) * t267, t3 * qJD(1) + ((t249 * t201 + t86 * t93) * m(6) + t316) * qJD(4) + t7, qJD(1) * t5 + qJD(4) * t8 + t7;];
Cq = t9;
