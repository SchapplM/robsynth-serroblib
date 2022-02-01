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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:14
% EndTime: 2022-01-23 09:12:23
% DurationCPUTime: 5.76s
% Computational Cost: add. (13270->287), mult. (15985->405), div. (0->0), fcn. (17069->8), ass. (0->163)
t230 = qJ(1) + pkin(7);
t227 = sin(t230);
t228 = cos(t230);
t236 = cos(qJ(4));
t232 = cos(pkin(8));
t234 = sin(qJ(4));
t293 = t232 * t234;
t204 = t227 * t293 + t228 * t236;
t292 = t232 * t236;
t294 = t228 * t234;
t205 = t227 * t292 - t294;
t404 = rSges(6,1) + pkin(4);
t121 = rSges(6,2) * t205 + t204 * t404;
t171 = -rSges(5,1) * t204 - rSges(5,2) * t205;
t297 = t227 * t234;
t207 = t228 * t292 + t297;
t375 = -t227 * t236 + t228 * t293;
t173 = -rSges(5,1) * t375 - rSges(5,2) * t207;
t276 = t207 * rSges(6,2) + t375 * t404;
t341 = m(6) / 0.2e1;
t343 = m(5) / 0.2e1;
t312 = (t121 * t227 + t228 * t276) * t341 + (-t171 * t227 - t173 * t228) * t343;
t251 = t207 * rSges(5,1) - rSges(5,2) * t375;
t231 = sin(pkin(8));
t295 = t228 * t231;
t154 = rSges(5,3) * t295 + t251;
t203 = -rSges(5,3) * t232 + (rSges(5,1) * t236 - rSges(5,2) * t234) * t231;
t106 = t154 * t232 + t203 * t295;
t298 = t227 * t231;
t360 = -t205 * rSges(5,1) + t204 * rSges(5,2);
t150 = rSges(5,3) * t298 - t360;
t393 = t150 * t232 + t203 * t298;
t233 = -qJ(5) - pkin(6);
t316 = pkin(6) + t233;
t319 = pkin(3) * t232;
t226 = pkin(4) * t236 + pkin(3);
t299 = t226 * t232;
t401 = -t205 * rSges(6,1) + t204 * rSges(6,2) + pkin(4) * t294 - t227 * t299;
t279 = -(t316 * t231 + t319) * t227 + rSges(6,3) * t298 - t401;
t317 = -pkin(3) + t226;
t378 = (t316 - rSges(6,3)) * t232 + (rSges(6,1) * t236 - rSges(6,2) * t234 + t317) * t231;
t400 = t279 * t232 + t378 * t298;
t359 = t207 * rSges(6,1) - rSges(6,2) * t375 - t233 * t295;
t278 = pkin(4) * t297 + (-pkin(6) * t231 + t317 * t232) * t228 + rSges(6,3) * t295 + t359;
t69 = t278 * t232 + t378 * t295;
t314 = (-t69 * t227 + t400 * t228) * t341 + (-t106 * t227 + t393 * t228) * t343;
t4 = t314 - t312;
t409 = t4 * qJD(1);
t397 = Icges(5,5) + Icges(6,5);
t396 = -Icges(5,6) - Icges(6,6);
t395 = Icges(6,3) + Icges(5,3);
t210 = (-Icges(6,5) * t234 - Icges(6,6) * t236) * t231;
t211 = (-Icges(5,5) * t234 - Icges(5,6) * t236) * t231;
t304 = Icges(5,4) * t234;
t201 = -Icges(5,5) * t232 + (Icges(5,1) * t236 - t304) * t231;
t213 = (-Icges(5,2) * t236 - t304) * t231;
t269 = t201 + t213;
t301 = Icges(6,4) * t234;
t200 = -Icges(6,5) * t232 + (Icges(6,1) * t236 - t301) * t231;
t212 = (-Icges(6,2) * t236 - t301) * t231;
t270 = t200 + t212;
t303 = Icges(5,4) * t236;
t199 = -Icges(5,6) * t232 + (-Icges(5,2) * t234 + t303) * t231;
t215 = (-Icges(5,1) * t234 - t303) * t231;
t271 = -t199 + t215;
t300 = Icges(6,4) * t236;
t198 = -Icges(6,6) * t232 + (-Icges(6,2) * t234 + t300) * t231;
t214 = (-Icges(6,1) * t234 - t300) * t231;
t272 = -t198 + t214;
t408 = ((t211 + t210) * t232 + ((-t271 - t272) * t236 + (t269 + t270) * t234) * t231) * t232;
t407 = t396 * t204 + t397 * t205 + t395 * t298;
t406 = -t231 * (-(t201 / 0.2e1 + t213 / 0.2e1 + t200 / 0.2e1 + t212 / 0.2e1) * t234 + (t215 / 0.2e1 - t199 / 0.2e1 + t214 / 0.2e1 - t198 / 0.2e1) * t236) + (t211 / 0.2e1 + t210 / 0.2e1) * t232;
t183 = Icges(6,4) * t205;
t135 = -Icges(6,2) * t204 + Icges(6,6) * t298 + t183;
t186 = Icges(5,4) * t205;
t138 = -Icges(5,2) * t204 + Icges(5,6) * t298 + t186;
t182 = Icges(6,4) * t204;
t142 = -Icges(6,1) * t205 - Icges(6,5) * t298 + t182;
t185 = Icges(5,4) * t204;
t145 = -Icges(5,1) * t205 - Icges(5,5) * t298 + t185;
t405 = (t142 + t145) * t207 + t375 * (t135 + t138);
t255 = -sin(qJ(1)) * pkin(1) + t228 * qJ(3);
t377 = (-pkin(2) + (-rSges(6,3) + t233) * t231) * t227 + t255 + t401;
t318 = pkin(4) * t234;
t320 = cos(qJ(1)) * pkin(1);
t93 = t320 + (qJ(3) + t318) * t227 + (rSges(6,3) * t231 + pkin(2) + t299) * t228 + t359;
t399 = m(6) * (t93 * t295 - t298 * t377);
t394 = t407 * t295;
t392 = t394 - t405;
t391 = qJD(1) * t399;
t390 = t121 * t228 - t276 * t227;
t302 = Icges(6,4) * t207;
t137 = -Icges(6,2) * t375 + Icges(6,6) * t295 + t302;
t305 = Icges(5,4) * t207;
t140 = -Icges(5,2) * t375 + Icges(5,6) * t295 + t305;
t184 = Icges(6,4) * t375;
t143 = Icges(6,1) * t207 + Icges(6,5) * t295 - t184;
t187 = Icges(5,4) * t375;
t146 = Icges(5,1) * t207 + Icges(5,5) * t295 - t187;
t387 = (t143 + t146) * t207 + (-t137 - t140) * t375 + (t207 * t397 + t295 * t395 + t375 * t396) * t295;
t352 = (rSges(5,3) + pkin(6)) * t231 + pkin(2) + t319;
t376 = -t227 * t352 + t255 + t360;
t373 = t392 - t394;
t370 = t231 / 0.2e1;
t369 = -t232 / 0.2e1;
t366 = m(6) * t231;
t280 = -Icges(5,2) * t207 + t146 - t187;
t282 = -Icges(6,2) * t207 + t143 - t184;
t358 = t280 + t282;
t281 = -Icges(5,2) * t205 - t145 - t185;
t283 = -Icges(6,2) * t205 - t142 - t182;
t357 = t281 + t283;
t284 = -Icges(5,1) * t375 - t140 - t305;
t286 = -Icges(6,1) * t375 - t137 - t302;
t356 = t284 + t286;
t285 = -Icges(5,1) * t204 - t138 - t186;
t287 = -Icges(6,1) * t204 - t135 - t183;
t355 = t285 + t287;
t158 = -Icges(6,5) * t204 - Icges(6,6) * t205;
t159 = -Icges(6,5) * t375 - Icges(6,6) * t207;
t160 = -Icges(5,5) * t204 - Icges(5,6) * t205;
t161 = -Icges(5,5) * t375 - Icges(5,6) * t207;
t350 = (t160 + t158) * t298 + (t159 + t161) * t295;
t346 = 0.4e1 * qJD(1);
t345 = 2 * qJD(4);
t339 = m(4) * ((rSges(4,2) * t298 + rSges(4,3) * t228 + t255) * t228 + (t320 - rSges(4,2) * t295 + (rSges(4,3) + qJ(3)) * t227) * t227);
t103 = t227 * qJ(3) + t228 * t352 + t251 + t320;
t336 = m(5) * (t103 * t173 - t171 * t376);
t335 = m(5) * (t103 * t227 + t376 * t228);
t329 = m(6) * (-t69 * t295 - t298 * t400);
t327 = m(6) * (t121 * t377 - t276 * t93);
t326 = m(6) * (t93 * t227 + t377 * t228);
t325 = t390 * t366;
t157 = (-t227 ^ 2 - t228 ^ 2) * t366;
t268 = t157 * qJD(1);
t248 = (t373 + t405) * t370 * t228 + (t369 + t232 / 0.2e1) * ((-t395 * t232 + (t234 * t396 + t236 * t397) * t231) * t298 + (t201 + t200) * t205 + (-t199 - t198) * t204);
t247 = -(t392 * t227 + t387 * t228) * t231 / 0.2e1 + ((t407 * t298 + t387) * t228 + t373 * t227) * t370;
t246 = (rSges(6,1) * t234 + rSges(6,2) * t236 + t318) * t231 ^ 2;
t217 = (-rSges(5,1) * t234 - rSges(5,2) * t236) * t231;
t120 = -t173 * t232 - t217 * t295;
t119 = t171 * t232 + t217 * t298;
t90 = (t171 * t228 - t173 * t227) * t231;
t89 = t246 * t228 + t276 * t232;
t88 = -t121 * t232 - t246 * t227;
t71 = t325 / 0.2e1;
t70 = t390 * t231;
t62 = t271 * t207 + t211 * t295 - t269 * t375;
t61 = t272 * t207 + t210 * t295 - t270 * t375;
t60 = -t269 * t204 + t271 * t205 + t211 * t298;
t59 = -t270 * t204 + t272 * t205 + t210 * t298;
t38 = -t161 * t232 + (-t280 * t234 + t284 * t236) * t231;
t37 = -t160 * t232 + (-t281 * t234 + t285 * t236) * t231;
t36 = -t159 * t232 + (-t282 * t234 + t286 * t236) * t231;
t35 = -t158 * t232 + (-t283 * t234 + t287 * t236) * t231;
t26 = t329 / 0.2e1;
t23 = t326 + t335 + t339;
t15 = t327 + t336 - t406;
t7 = t26 + t71;
t6 = t26 - t325 / 0.2e1;
t5 = t71 - t329 / 0.2e1;
t2 = t312 + t314;
t1 = (t247 * t227 + t248 * t228) * t231;
t3 = [t23 * qJD(3) + t15 * qJD(4) + qJD(5) * t399, 0, qJD(1) * t23 + qJD(4) * t2, t15 * qJD(1) + t2 * qJD(3) + t7 * qJD(5) + ((t121 * t400 + t276 * t69 + t377 * t88 + t89 * t93) * t341 + (t103 * t120 - t106 * t173 + t119 * t376 - t171 * t393) * t343) * t345 + (t408 + ((t62 / 0.2e1 + t61 / 0.2e1 + t38 / 0.2e1 + t36 / 0.2e1 - t248) * t228 + (t37 / 0.2e1 + t35 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1 - t247) * t227) * t231) * qJD(4), t7 * qJD(4) + t391; 0, 0, 0, (-t70 * t341 + t90 * t343) * t345, 0; -t4 * qJD(4) + t157 * qJD(5) + (-t326 / 0.4e1 - t335 / 0.4e1 - t339 / 0.4e1) * t346, 0, 0, -t409 + ((t119 * t227 - t120 * t228) * t343 + (t227 * t88 - t228 * t89) * t341) * t345, t268; t4 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t327 / 0.4e1 - t336 / 0.4e1) * t346 + t406 * qJD(1), 0, t409, t1 * qJD(1) + ((t408 + ((t36 + t38) * t228 + (t35 + t37) * t227) * t231) * t369 + ((-t204 * t358 + t205 * t356) * t295 + (-t60 - t59) * t232 + (-t357 * t204 + t355 * t205 + t350) * t298) * t298 / 0.2e1 + ((t207 * t355 - t357 * t375) * t298 + (-t62 - t61) * t232 + (t356 * t207 - t358 * t375 + t350) * t295) * t295 / 0.2e1 + (t393 * t119 - t106 * t120 + (t150 * t228 - t154 * t227) * t231 * t90) * m(5) + (-(-t278 * t227 + t279 * t228) * t231 * t70 + t400 * t88 - t69 * t89) * m(6)) * qJD(4), t6 * qJD(1); -t157 * qJD(3) + t5 * qJD(4) - t391, 0, -t268, t5 * qJD(1) + m(6) * (t232 * t70 + (t227 * t89 + t228 * t88) * t231) * qJD(4), 0;];
Cq = t3;
