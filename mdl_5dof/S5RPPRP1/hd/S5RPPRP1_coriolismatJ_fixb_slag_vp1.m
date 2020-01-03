% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:10
% EndTime: 2020-01-03 11:25:26
% DurationCPUTime: 6.18s
% Computational Cost: add. (13270->297), mult. (15985->424), div. (0->0), fcn. (17069->8), ass. (0->169)
t230 = qJ(1) + pkin(7);
t226 = sin(t230);
t227 = cos(t230);
t236 = cos(qJ(4));
t232 = cos(pkin(8));
t234 = sin(qJ(4));
t288 = t232 * t234;
t196 = -t226 * t288 - t227 * t236;
t287 = t232 * t236;
t290 = t227 * t234;
t197 = t226 * t287 - t290;
t375 = rSges(6,1) + pkin(4);
t119 = -rSges(6,2) * t197 + t196 * t375;
t166 = rSges(5,1) * t196 - rSges(5,2) * t197;
t198 = -t226 * t236 + t227 * t288;
t292 = t226 * t234;
t199 = t227 * t287 + t292;
t168 = rSges(5,1) * t198 + rSges(5,2) * t199;
t271 = rSges(6,2) * t199 + t198 * t375;
t339 = m(6) / 0.2e1;
t340 = m(5) / 0.2e1;
t312 = (-t119 * t226 + t227 * t271) * t339 + (-t166 * t226 + t168 * t227) * t340;
t252 = rSges(5,1) * t197 + rSges(5,2) * t196;
t231 = sin(pkin(8));
t293 = t226 * t231;
t146 = rSges(5,3) * t293 + t252;
t195 = -rSges(5,3) * t232 + (rSges(5,1) * t236 - rSges(5,2) * t234) * t231;
t102 = t146 * t232 + t195 * t293;
t291 = t227 * t231;
t148 = rSges(5,1) * t199 - rSges(5,2) * t198 + rSges(5,3) * t291;
t380 = -t148 * t232 - t195 * t291;
t233 = -qJ(5) - pkin(6);
t316 = pkin(6) + t233;
t225 = pkin(4) * t236 + pkin(3);
t317 = -pkin(3) + t225;
t270 = (t316 - rSges(6,3)) * t232 + (rSges(6,1) * t236 - rSges(6,2) * t234 + t317) * t231;
t319 = pkin(3) * t232;
t294 = t225 * t232;
t388 = -rSges(6,1) * t199 + rSges(6,2) * t198 - rSges(6,3) * t291 - pkin(4) * t292 - t227 * t294;
t272 = (t231 * t316 + t319) * t227 + t388;
t389 = t232 * t272 - t270 * t291;
t318 = pkin(6) * t231;
t289 = t231 * t233;
t354 = rSges(6,1) * t197 + rSges(6,2) * t196 - pkin(4) * t290 - t226 * t289;
t273 = (t232 * t317 - t318) * t226 + rSges(6,3) * t293 + t354;
t67 = t232 * t273 + t270 * t293;
t314 = (t226 * t389 + t227 * t67) * t339 + (t102 * t227 + t226 * t380) * t340;
t4 = t314 - t312;
t404 = t4 * qJD(1);
t402 = Icges(5,5) + Icges(6,5);
t401 = Icges(5,6) + Icges(6,6);
t400 = Icges(5,3) + Icges(6,3);
t391 = -t198 * t401 + t199 * t402 + t291 * t400;
t403 = t391 * t291;
t202 = (-Icges(6,5) * t234 - Icges(6,6) * t236) * t231;
t203 = (-Icges(5,5) * t234 - Icges(5,6) * t236) * t231;
t304 = Icges(5,4) * t234;
t193 = -Icges(5,5) * t232 + (Icges(5,1) * t236 - t304) * t231;
t205 = (-Icges(5,2) * t236 - t304) * t231;
t266 = t193 + t205;
t301 = Icges(6,4) * t234;
t192 = -Icges(6,5) * t232 + (Icges(6,1) * t236 - t301) * t231;
t204 = (-Icges(6,2) * t236 - t301) * t231;
t267 = t192 + t204;
t303 = Icges(5,4) * t236;
t191 = -Icges(5,6) * t232 + (-Icges(5,2) * t234 + t303) * t231;
t207 = (-Icges(5,1) * t234 - t303) * t231;
t268 = t191 - t207;
t300 = Icges(6,4) * t236;
t190 = -Icges(6,6) * t232 + (-Icges(6,2) * t234 + t300) * t231;
t206 = (-Icges(6,1) * t234 - t300) * t231;
t269 = t190 - t206;
t399 = ((t203 + t202) * t232 + ((t268 + t269) * t236 + (t266 + t267) * t234) * t231) * t232;
t176 = Icges(6,4) * t199;
t133 = Icges(6,2) * t198 - Icges(6,6) * t291 - t176;
t179 = Icges(5,4) * t199;
t136 = Icges(5,2) * t198 - Icges(5,6) * t291 - t179;
t175 = Icges(6,4) * t198;
t138 = Icges(6,1) * t199 + Icges(6,5) * t291 - t175;
t178 = Icges(5,4) * t198;
t141 = Icges(5,1) * t199 + Icges(5,5) * t291 - t178;
t396 = -t196 * (t133 + t136) + (t141 + t138) * t197;
t395 = -t231 * (-(t193 / 0.2e1 + t205 / 0.2e1 + t192 / 0.2e1 + t204 / 0.2e1) * t234 + (t207 / 0.2e1 - t191 / 0.2e1 + t206 / 0.2e1 - t190 / 0.2e1) * t236) + (t203 / 0.2e1 + t202 / 0.2e1) * t232;
t392 = t196 * t401 + t197 * t402 + t293 * t400;
t373 = t293 * t391 + t396;
t382 = t391 * t226;
t381 = t392 * t227;
t379 = t119 * t227 + t226 * t271;
t378 = t392 * t291;
t302 = Icges(6,4) * t197;
t131 = Icges(6,2) * t196 + Icges(6,6) * t293 + t302;
t305 = Icges(5,4) * t197;
t134 = Icges(5,2) * t196 + Icges(5,6) * t293 + t305;
t174 = Icges(6,4) * t196;
t137 = Icges(6,1) * t197 + Icges(6,5) * t293 + t174;
t177 = Icges(5,4) * t196;
t140 = Icges(5,1) * t197 + Icges(5,5) * t293 + t177;
t374 = t392 * t293 + (t137 + t140) * t197 + (t131 + t134) * t196;
t153 = Icges(6,5) * t196 - Icges(6,6) * t197;
t154 = Icges(6,5) * t198 + Icges(6,6) * t199;
t155 = Icges(5,5) * t196 - Icges(5,6) * t197;
t156 = Icges(5,5) * t198 + Icges(5,6) * t199;
t368 = -(t156 + t154) * t291 + (t155 + t153) * t293;
t365 = t231 / 0.2e1;
t364 = -t232 / 0.2e1;
t361 = m(6) * t231;
t274 = Icges(5,2) * t199 - t141 + t178;
t276 = Icges(6,2) * t199 - t138 + t175;
t353 = -t274 - t276;
t275 = -Icges(5,2) * t197 + t140 + t177;
t277 = -Icges(6,2) * t197 + t137 + t174;
t352 = t275 + t277;
t278 = -Icges(5,1) * t198 + t136 - t179;
t280 = -Icges(6,1) * t198 + t133 - t176;
t351 = t278 + t280;
t279 = -Icges(5,1) * t196 + t134 + t305;
t281 = -Icges(6,1) * t196 + t131 + t302;
t350 = t279 + t281;
t260 = t227 * pkin(2) + t226 * qJ(3) + cos(qJ(1)) * pkin(1);
t349 = -t227 * t289 + t260 - t388;
t348 = (t318 + t319) * t227 + t260 + t148;
t345 = t227 ^ 2;
t343 = 0.4e1 * qJD(1);
t342 = 2 * qJD(4);
t255 = sin(qJ(1)) * pkin(1) - t227 * qJ(3);
t337 = m(4) * (-(-rSges(4,2) * t293 - t227 * rSges(4,3) + pkin(2) * t226 + t255) * t227 + (-rSges(4,2) * t291 + rSges(4,3) * t226 + t260) * t226);
t100 = (t319 + pkin(2) + (rSges(5,3) + pkin(6)) * t231) * t226 + t252 + t255;
t334 = m(5) * (t100 * t166 - t168 * t348);
t333 = m(5) * (-t100 * t227 + t226 * t348);
t29 = t291 * t389 - t293 * t67;
t329 = m(6) * ((t226 * t67 - t227 * t389) * t231 + t29);
t327 = m(6) * t29;
t92 = (rSges(6,3) * t231 + pkin(2) + t294) * t226 + t255 + t354;
t325 = m(6) * (t119 * t92 - t271 * t349);
t324 = m(6) * (t226 * t349 - t227 * t92);
t323 = t379 * t361;
t43 = t291 * t349 + t293 * t92;
t310 = m(6) * qJD(1);
t151 = (t226 ^ 2 + t345) * t361;
t263 = t151 * qJD(1);
t248 = -(((t381 - t382) * t231 + t373 - t378) * t227 + (t374 + t403) * t226) * t231 / 0.2e1 + (t226 * t374 + t227 * t373) * t365;
t247 = (t391 * t231 * t345 - t403 * t227 + ((t381 + t382) * t231 - t378 - t373 + t396) * t226) * t365 + (t364 + t232 / 0.2e1) * ((-t400 * t232 + (-t234 * t401 + t236 * t402) * t231) * t291 + (t193 + t192) * t199 + (-t191 - t190) * t198);
t246 = (rSges(6,2) * t236 + t234 * t375) * t231 ^ 2;
t209 = (-rSges(5,1) * t234 - rSges(5,2) * t236) * t231;
t118 = t168 * t232 - t209 * t291;
t117 = -t166 * t232 - t209 * t293;
t90 = (t166 * t227 + t168 * t226) * t231;
t89 = t227 * t246 + t232 * t271;
t88 = -t119 * t232 + t226 * t246;
t71 = -t323 / 0.2e1;
t70 = t379 * t231;
t62 = t198 * t266 + t199 * t268 - t203 * t291;
t61 = t198 * t267 + t199 * t269 - t202 * t291;
t60 = t196 * t266 - t197 * t268 + t203 * t293;
t59 = t196 * t267 - t197 * t269 + t202 * t293;
t38 = -t156 * t232 + (-t234 * t274 - t236 * t278) * t231;
t37 = -t155 * t232 + (-t234 * t275 - t236 * t279) * t231;
t36 = -t154 * t232 + (-t234 * t276 - t236 * t280) * t231;
t35 = -t153 * t232 + (-t234 * t277 - t236 * t281) * t231;
t26 = t327 / 0.2e1;
t23 = t324 + t333 + t337;
t15 = t325 + t334 - t395;
t13 = t329 / 0.2e1;
t7 = t26 + t71 - t329 / 0.2e1;
t6 = t26 + t13 + t323 / 0.2e1;
t5 = t71 + t13 - t327 / 0.2e1;
t3 = t312 + t314;
t1 = (t226 * t247 + t227 * t248) * t231;
t2 = [m(6) * t43 * qJD(5) + qJD(3) * t23 + qJD(4) * t15, 0, qJD(1) * t23 + qJD(4) * t3, t15 * qJD(1) + t3 * qJD(3) + t7 * qJD(5) + ((t100 * t117 - t102 * t166 + t118 * t348 - t168 * t380) * t340 + (-t119 * t67 - t271 * t389 + t349 * t89 + t88 * t92) * t339) * t342 + (t399 + ((-t61 / 0.2e1 - t38 / 0.2e1 - t36 / 0.2e1 - t62 / 0.2e1 - t248) * t227 + (t37 / 0.2e1 + t35 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1 - t247) * t226) * t231) * qJD(4), qJD(4) * t7 + t310 * t43; 0, 0, 0, (t339 * t70 + t340 * t90) * t342, 0; -t4 * qJD(4) - t151 * qJD(5) + (-t337 / 0.4e1 - t333 / 0.4e1 - t324 / 0.4e1) * t343, 0, 0, -t404 + ((-t117 * t226 - t118 * t227) * t340 + (-t226 * t88 - t227 * t89) * t339) * t342, -t263; t4 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t334 / 0.4e1 - t325 / 0.4e1) * t343 + t395 * qJD(1), 0, t404, t1 * qJD(1) + ((t399 + ((-t36 - t38) * t227 + (t35 + t37) * t226) * t231) * t364 + ((t196 * t353 + t197 * t351) * t291 + (-t60 - t59) * t232 + (t352 * t196 - t350 * t197 + t368) * t293) * t293 / 0.2e1 - ((t198 * t352 + t199 * t350) * t293 + (-t62 - t61) * t232 + (t353 * t198 - t351 * t199 - t368) * t291) * t291 / 0.2e1 + m(5) * (-t102 * t117 + t380 * t118 + (t146 * t227 - t148 * t226) * t231 * t90) + ((t226 * t272 + t227 * t273) * t231 * t70 - t67 * t88 + t389 * t89) * m(6)) * qJD(4), t6 * qJD(1); (-t226 * t92 - t227 * t349) * t231 * t310 + t151 * qJD(3) + t5 * qJD(4), 0, t263, t5 * qJD(1) + m(6) * (-t232 * t70 + (t226 * t89 - t227 * t88) * t231) * qJD(4), 0;];
Cq = t2;
