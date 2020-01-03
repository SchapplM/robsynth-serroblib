% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:59:01
% DurationCPUTime: 4.77s
% Computational Cost: add. (4003->323), mult. (9311->438), div. (0->0), fcn. (8490->4), ass. (0->202)
t236 = sin(qJ(1));
t238 = cos(qJ(1));
t235 = sin(qJ(2));
t224 = Icges(4,5) * t235;
t237 = cos(qJ(2));
t317 = Icges(4,1) * t237;
t253 = t224 + t317;
t144 = Icges(4,4) * t236 + t238 * t253;
t315 = Icges(3,4) * t235;
t195 = Icges(3,1) * t237 - t315;
t146 = Icges(3,5) * t236 + t195 * t238;
t384 = t144 + t146;
t196 = pkin(2) * t235 - qJ(3) * t237;
t321 = rSges(4,1) * t235;
t198 = -rSges(4,3) * t237 + t321;
t282 = t196 + t198;
t104 = t282 * t236;
t106 = t282 * t238;
t304 = t235 * t238;
t306 = t235 * t236;
t301 = t237 * t238;
t302 = t236 * t237;
t231 = t238 * rSges(4,2);
t232 = t238 * pkin(5);
t319 = rSges(4,3) + qJ(3);
t333 = rSges(4,1) + pkin(2);
t359 = t319 * t235 + t333 * t237 + pkin(1);
t83 = -t236 * t359 + t231 + t232;
t84 = (rSges(4,2) + pkin(5)) * t236 + t359 * t238;
t325 = t83 * t301 + t84 * t302;
t332 = rSges(5,1) + pkin(3);
t276 = pkin(2) + t332;
t320 = rSges(5,2) + qJ(3);
t360 = t320 * t235 + t276 * t237 + pkin(1);
t370 = rSges(5,3) + qJ(4);
t73 = -t236 * t360 - t238 * t370 + t232;
t74 = t238 * t360 + (pkin(5) - t370) * t236;
t326 = t73 * t301 + t74 * t302;
t354 = m(5) / 0.2e1;
t355 = m(4) / 0.2e1;
t197 = rSges(5,1) * t235 - rSges(5,2) * t237;
t331 = pkin(3) * t235;
t258 = t196 + t197 + t331;
t93 = t258 * t236;
t95 = t258 * t238;
t328 = (-t93 * t304 + t95 * t306 + t326) * t354 + (-t104 * t304 + t106 * t306 + t325) * t355;
t218 = pkin(2) * t306;
t100 = t218 + (-t319 * t237 + t321) * t236;
t208 = qJ(3) * t301;
t216 = rSges(4,3) * t301;
t101 = -t333 * t304 + t208 + t216;
t88 = t218 + (t332 * t235 - t320 * t237) * t236;
t217 = rSges(5,2) * t301;
t89 = -t276 * t304 + t208 + t217;
t329 = ((t236 * t89 + t238 * t88) * t235 + t326) * t354 + ((t100 * t238 + t101 * t236) * t235 + t325) * t355;
t2 = t329 - t328;
t383 = t2 * qJD(1);
t234 = t238 ^ 2;
t211 = Icges(4,5) * t301;
t132 = Icges(4,6) * t236 + Icges(4,3) * t304 + t211;
t184 = Icges(3,5) * t237 - Icges(3,6) * t235;
t134 = Icges(3,3) * t236 + t184 * t238;
t187 = Icges(4,4) * t237 + Icges(4,6) * t235;
t138 = Icges(4,2) * t236 + t187 * t238;
t382 = t132 * t304 + t384 * t301 + (t134 + t138) * t236;
t233 = t236 ^ 2;
t277 = t233 + t234;
t112 = t146 * t302;
t181 = Icges(5,5) * t237 + Icges(5,6) * t235;
t129 = Icges(5,3) * t238 + t236 * t181;
t381 = -t236 * t129 - t134 * t238 + t112;
t133 = Icges(3,5) * t302 - Icges(3,6) * t306 - Icges(3,3) * t238;
t314 = Icges(5,4) * t237;
t186 = Icges(5,2) * t235 + t314;
t135 = Icges(5,6) * t238 + t186 * t236;
t226 = Icges(5,4) * t235;
t316 = Icges(5,1) * t237;
t252 = t226 + t316;
t141 = Icges(5,5) * t238 + t236 * t252;
t213 = Icges(3,4) * t306;
t145 = Icges(3,1) * t302 - Icges(3,5) * t238 - t213;
t380 = -t236 * t133 - t135 * t304 + (-t141 - t145) * t301;
t379 = (Icges(3,6) - Icges(4,6) + Icges(5,6)) * t237 + (Icges(4,4) + Icges(3,5) - Icges(5,5)) * t235;
t142 = -Icges(5,5) * t236 + t238 * t252;
t188 = Icges(3,2) * t237 + t315;
t310 = Icges(4,3) * t237;
t248 = t310 - t224;
t311 = Icges(5,2) * t237;
t250 = t311 - t226;
t378 = t142 + (-t188 - t248 - t250) * t238 + t384;
t245 = t135 * t235 + t141 * t237;
t308 = (-Icges(4,2) * t238 + t236 * t187) * t238;
t313 = Icges(4,5) * t237;
t183 = Icges(4,3) * t235 + t313;
t131 = -Icges(4,6) * t238 + t183 * t236;
t143 = -Icges(4,4) * t238 + t236 * t253;
t369 = (t131 * t235 + t143 * t237) * t236;
t377 = t129 * t238 + t236 * t245 - t308 + t369;
t228 = Icges(3,4) * t237;
t312 = Icges(3,2) * t235;
t140 = Icges(3,6) * t236 + (t228 - t312) * t238;
t376 = -t140 * t304 + t382;
t375 = t308 + t382;
t139 = Icges(3,4) * t302 - Icges(3,2) * t306 - Icges(3,6) * t238;
t374 = -t139 * t304 - t140 * t306 - t380 + t381;
t373 = -t236 / 0.2e1;
t335 = t236 / 0.2e1;
t372 = -t238 / 0.2e1;
t305 = t235 * t237;
t280 = t277 * t305;
t368 = (m(4) / 0.4e1 + m(5) / 0.4e1) * (t280 - t305);
t202 = rSges(5,1) * t237 + rSges(5,2) * t235;
t367 = pkin(3) * t237 + t202;
t366 = qJD(2) * t238;
t365 = t379 * t236;
t364 = t379 * t238;
t363 = t378 * t236;
t318 = Icges(3,1) * t235;
t254 = -t228 - t318;
t267 = (t254 * t238 - t140) * t236;
t212 = Icges(5,4) * t301;
t136 = Icges(5,2) * t304 - Icges(5,6) * t236 + t212;
t269 = (-Icges(5,1) * t304 + t136 + t212) * t236;
t271 = (-Icges(4,1) * t304 + t132 + t211) * t236;
t362 = t267 + t269 + t271;
t190 = Icges(5,1) * t235 - t314;
t192 = Icges(4,1) * t235 - t313;
t361 = -t235 * (t195 / 0.2e1 - t188 / 0.2e1 + t224 + t317 / 0.2e1 - t310 / 0.2e1 + t226 + t316 / 0.2e1 - t311 / 0.2e1) - t237 * (t228 + t318 / 0.2e1 - t312 / 0.2e1 + t192 / 0.2e1 - t183 / 0.2e1 + t190 / 0.2e1 - t186 / 0.2e1);
t358 = 0.4e1 * qJD(1);
t357 = 0.2e1 * qJD(2);
t322 = rSges(3,1) * t237;
t273 = pkin(1) + t322;
t279 = rSges(3,2) * t306 + t238 * rSges(3,3);
t102 = -t236 * t273 + t232 + t279;
t215 = rSges(3,2) * t304;
t103 = -t215 + t273 * t238 + (rSges(3,3) + pkin(5)) * t236;
t199 = rSges(3,1) * t235 + rSges(3,2) * t237;
t176 = t199 * t236;
t178 = t199 * t238;
t353 = m(3) * (t102 * t176 - t103 * t178);
t179 = t277 * t235;
t323 = -t104 * t302 - t106 * t301;
t203 = rSges(4,1) * t237 + rSges(4,3) * t235;
t201 = pkin(2) * t237 + qJ(3) * t235;
t284 = t277 * t201;
t49 = t236 * (t203 * t236 - t231) + (t236 * rSges(4,2) + t203 * t238) * t238 + t284;
t349 = m(4) * (t179 * t49 + t323);
t347 = m(4) * (t100 * t83 + t101 * t84);
t346 = m(4) * (t84 * t304 - t83 * t306);
t324 = -t95 * t301 - t93 * t302;
t43 = (pkin(3) * t302 + t202 * t236) * t236 + t284 + t234 * t367;
t344 = m(5) * (t179 * t43 + t324);
t340 = m(5) * (t73 * t88 + t74 * t89);
t339 = m(5) * (t74 * t304 - t73 * t306);
t338 = m(5) * (-t74 * t236 - t238 * t73);
t337 = m(5) * (-t236 * t88 + t238 * t89);
t336 = m(5) * (t236 * t93 + t238 * t95);
t307 = t139 * t235;
t90 = m(5) * t179;
t299 = t90 * qJD(1);
t297 = t136 * t304 + t142 * t301;
t294 = -t192 * t236 + t131;
t293 = -t190 * t236 + t135;
t292 = -t254 * t236 + t139;
t291 = -t250 * t236 + t141;
t289 = -t248 * t236 + t143;
t287 = -Icges(3,2) * t302 + t145 - t213;
t285 = t236 * (qJ(3) * t302 - t218) + t238 * (-pkin(2) * t304 + t208);
t281 = -t201 - t203;
t274 = qJD(1) * t338;
t272 = t294 * t238;
t270 = t293 * t238;
t268 = t292 * t238;
t266 = t291 * t238;
t264 = t289 * t238;
t262 = t287 * t238;
t259 = t140 * t235 - t133;
t257 = -t201 - t367;
t256 = -t132 * t306 + t138 * t238 - t144 * t302;
t255 = -t181 / 0.2e1 + t187 / 0.2e1 + t184 / 0.2e1;
t242 = t256 * t373 + t129 * t234 / 0.2e1 + (t238 * t259 - t375 + t376) * t238 / 0.2e1 + (-t236 * (-t145 * t237 + t307) - t133 * t238 + t377) * t372 + (t131 * t304 + t143 * t301 + t236 * t259 + t256 + t374 - t381) * t335;
t130 = -Icges(5,3) * t236 + t181 * t238;
t241 = ((-t130 - t245) * t236 + t297 - t369 + t375 + t377) * t373 + (-t236 * t130 + t297 + t376) * t335 + (-t112 + (t134 + t307) * t238 + t374 + t380) * t372;
t204 = -rSges(3,2) * t235 + t322;
t107 = t281 * t238;
t105 = t281 * t236;
t96 = t257 * t238;
t94 = t257 * t236;
t72 = 0.4e1 * t368;
t68 = t238 * (-rSges(4,1) * t304 + t216) - t198 * t233 + t285;
t53 = t336 / 0.2e1;
t48 = t337 / 0.2e1;
t46 = t238 * (-rSges(5,1) * t304 + t217) - t197 * t233 - t277 * t331 + t285;
t18 = t339 + t346;
t17 = t53 - t337 / 0.2e1;
t16 = t53 + t48;
t15 = t48 - t336 / 0.2e1;
t12 = t344 + t349;
t5 = t340 + t347 + t353 - t361;
t3 = t328 + t329;
t1 = t236 * t242 + t238 * t241;
t4 = [t5 * qJD(2) + t18 * qJD(3) + qJD(4) * t338, t5 * qJD(1) + t3 * qJD(3) + t16 * qJD(4) + ((t73 * t96 + t74 * t94 - t88 * t95 - t89 * t93) * t354 + (-t100 * t106 - t101 * t104 + t105 * t84 + t107 * t83) * t355) * t357 + (m(3) * (-t102 * t204 - t176 * t199) + t255 * t238 - t241) * t366 + ((m(3) * (-t103 * t204 + t178 * t199) + t255 * t236 - t242) * t236 + (t269 / 0.2e1 + t271 / 0.2e1 + t267 / 0.2e1 - t270 / 0.2e1 - t272 / 0.2e1 + t268 / 0.2e1) * t235 + (-t266 / 0.2e1 - t264 / 0.2e1 - t262 / 0.2e1 + t378 * t335) * t237) * qJD(2), qJD(1) * t18 + qJD(2) * t3, t16 * qJD(2) + t274; t1 * qJD(2) - t2 * qJD(3) + t17 * qJD(4) + (-t353 / 0.4e1 - t340 / 0.4e1 - t347 / 0.4e1) * t358 + t361 * qJD(1), t1 * qJD(1) + t12 * qJD(3) - ((t364 * t238 + (-t270 - t272 + t268 + t362) * t237 + ((t287 + t289 + t291) * t238 - t363) * t235) * t236 - t365 * t234) * t366 / 0.2e1 + (m(4) * (-t104 * t105 - t106 * t107 + t49 * t68) + m(5) * (t43 * t46 - t93 * t94 - t95 * t96) + m(3) * ((t236 * (rSges(3,1) * t302 - t279) + t238 * (rSges(3,1) * t301 + t236 * rSges(3,3) - t215)) * (-t236 * t176 - t178 * t238) + t277 * t204 * t199) + ((t365 * t236 + (t264 + t266 + t262 - t363) * t235 + ((t292 - t293 - t294) * t238 + t362) * t237) * t238 - t364 * t233) * t335) * qJD(2), -t383 + t12 * qJD(2) + (-0.4e1 * t368 + 0.2e1 * (t355 + t354) * (-t179 * t237 + t280)) * qJD(3), t17 * qJD(1); t2 * qJD(2) - t90 * qJD(4) + (-t346 / 0.4e1 - t339 / 0.4e1) * t358, t383 + t72 * qJD(3) + 0.4e1 * (-t349 / 0.4e1 - t344 / 0.4e1) * qJD(2) + ((-t237 * t68 + t323) * t355 + (-t237 * t46 + t324) * t354 + ((t105 * t236 + t107 * t238 + t49) * t355 + (t236 * t94 + t238 * t96 + t43) * t354) * t235) * t357, t72 * qJD(2), -t299; t15 * qJD(2) + t90 * qJD(3) - t274, t15 * qJD(1) + m(5) * (-t236 * t96 + t238 * t94) * qJD(2), t299, 0;];
Cq = t4;
