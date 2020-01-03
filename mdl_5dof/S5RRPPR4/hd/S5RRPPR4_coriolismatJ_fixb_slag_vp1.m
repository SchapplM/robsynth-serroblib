% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:41
% EndTime: 2019-12-31 19:27:45
% DurationCPUTime: 2.93s
% Computational Cost: add. (17558->231), mult. (16431->320), div. (0->0), fcn. (18550->8), ass. (0->171)
t339 = m(6) / 0.2e1;
t273 = qJ(1) + qJ(2);
t258 = sin(t273);
t259 = cos(t273);
t298 = sin(pkin(8));
t299 = cos(pkin(8));
t196 = -t258 * t299 + t259 * t298;
t225 = sin(qJ(5));
t227 = cos(qJ(5));
t249 = rSges(6,1) * t225 + rSges(6,2) * t227;
t359 = t249 * t196;
t195 = -t258 * t298 - t259 * t299;
t360 = t249 * t195;
t373 = m(6) * (t258 * t359 + t259 * t360);
t118 = -t373 / 0.2e1;
t131 = t373 / 0.2e1;
t372 = 0.2e1 * t339 * qJD(1);
t296 = Icges(6,4) * t227;
t246 = -Icges(6,2) * t225 + t296;
t247 = Icges(6,1) * t225 + t296;
t371 = t246 + t247;
t297 = Icges(6,4) * t225;
t245 = Icges(6,2) * t227 + t297;
t248 = Icges(6,1) * t227 - t297;
t255 = t371 * t227 / 0.2e1 + (-t245 / 0.2e1 + t248 / 0.2e1) * t225;
t232 = t259 * pkin(2) + t258 * qJ(3);
t229 = t259 * pkin(3) + t232;
t285 = t195 * t225;
t256 = rSges(6,2) * t285 + t196 * rSges(6,3);
t302 = rSges(6,1) * t227;
t260 = pkin(4) + t302;
t108 = t196 * pkin(7) - t260 * t195 + t229 + t256;
t306 = cos(qJ(1)) * pkin(1);
t103 = t108 + t306;
t241 = -t258 * pkin(2) + t259 * qJ(3);
t231 = -t258 * pkin(3) + t241;
t282 = t196 * t225;
t301 = t195 * rSges(6,3);
t257 = -rSges(6,2) * t282 + t301;
t107 = t195 * pkin(7) + t260 * t196 + t231 + t257;
t307 = sin(qJ(1)) * pkin(1);
t102 = t107 - t307;
t87 = t102 * t359;
t40 = t103 * t360 - t87;
t88 = t108 * t360;
t41 = -t107 * t359 + t88;
t251 = (t41 + t40) * t339 + t255;
t244 = Icges(6,5) * t227 - Icges(6,6) * t225;
t362 = t244 * t195;
t140 = Icges(6,3) * t196 - t362;
t368 = t196 * t140;
t361 = t244 * t196;
t138 = Icges(6,3) * t195 + t361;
t143 = Icges(6,6) * t196 - t246 * t195;
t367 = -t143 * t225 + t138;
t158 = -t195 * rSges(5,1) - t196 * rSges(5,2) + t229;
t155 = t158 + t306;
t350 = t196 * rSges(5,1) - t195 * rSges(5,2) + t231;
t353 = t350 - t307;
t65 = t155 * t350 - t158 * t353;
t141 = Icges(6,6) * t195 + t246 * t196;
t144 = Icges(6,5) * t195 + t248 * t196;
t366 = t141 * t225 - t144 * t227;
t364 = m(4) / 0.2e1;
t340 = m(5) / 0.2e1;
t316 = m(3) * (t306 * (-t258 * rSges(3,1) - t259 * rSges(3,2)) + (t259 * rSges(3,1) - t258 * rSges(3,2)) * t307);
t184 = -t258 * rSges(4,1) + t259 * rSges(4,3) + t241;
t180 = t184 - t307;
t185 = t259 * rSges(4,1) + t258 * rSges(4,3) + t232;
t181 = t185 + t306;
t314 = m(4) * (-t185 * t180 + t181 * t184);
t250 = -rSges(6,2) * t225 + t302;
t238 = t250 * t196;
t109 = t196 * pkin(4) + t238 + (rSges(6,3) + pkin(7)) * t195 + t231;
t358 = t107 - t109;
t128 = t180 * t259 + t181 * t258;
t148 = t184 * t259 + t185 * t258;
t254 = t108 * t258;
t240 = -t259 * t109 - t254;
t95 = t102 * t259;
t46 = t103 * t258 + t95;
t91 = t155 * t258 + t353 * t259;
t97 = t158 * t258 + t259 * t350;
t265 = (t240 + t46) * t339 + (-t97 + t91) * t340 + (t128 - t148) * t364;
t53 = t107 * t259 + t254;
t351 = (t240 + t53) * t339;
t355 = t265 - t351;
t352 = qJD(1) + qJD(2);
t281 = t196 * t227;
t149 = -rSges(6,1) * t281 - t257;
t320 = m(6) * (t149 + t238 + t301) * t195;
t264 = (t53 + t46) * t339 + (t97 + t91) * t340 + (t148 + t128) * t364;
t348 = t195 ^ 2;
t347 = t196 ^ 2;
t345 = 0.4e1 * qJD(1);
t344 = 0.2e1 * qJD(2);
t343 = 0.4e1 * qJD(2);
t334 = m(5) * t65;
t332 = m(5) * t91;
t331 = m(5) * t97;
t329 = m(6) * (-t87 - t88 - (-t103 * t195 - t109 * t196) * t249);
t21 = -t88 - (-t108 * t195 + t358 * t196) * t249;
t328 = m(6) * t21;
t32 = -t102 * t108 + t103 * t109;
t325 = m(6) * t32;
t33 = t358 * t108;
t324 = m(6) * t33;
t284 = t195 * t227;
t82 = t196 * t149 + t195 * (-rSges(6,1) * t284 + t256);
t323 = t82 * t320;
t322 = m(6) * t46;
t321 = m(6) * t53;
t319 = -t195 / 0.2e1;
t317 = t196 / 0.2e1;
t312 = m(4) * t128;
t311 = m(4) * t148;
t262 = t320 / 0.2e1;
t64 = t118 + t131;
t305 = t64 * qJD(3) + qJD(4) * t262;
t263 = -t320 / 0.2e1;
t304 = t352 * t263;
t303 = m(6) * qJD(5);
t291 = t360 * t249;
t286 = t195 * t250;
t146 = Icges(6,5) * t196 - t248 * t195;
t272 = -t195 * t140 - t146 * t281;
t271 = -t146 * t284 + t368;
t270 = -t247 * t196 - t141;
t269 = -t247 * t195 + t143;
t268 = -t245 * t196 + t144;
t267 = t245 * t195 + t146;
t51 = -t138 * t196 - t195 * t366;
t52 = t143 * t285 + t271;
t11 = (t51 + (t146 * t227 + t367) * t196) * t196 + (t367 * t195 - t271 + t368 + t52) * t195;
t50 = t143 * t282 + t272;
t12 = (t50 - t51 - t272) * t195 + t271 * t196;
t24 = -t195 * (t195 * t138 - t366 * t196) + t196 * t50;
t25 = -t195 * t51 + t196 * t52;
t2 = t323 + (t11 / 0.2e1 + t24 / 0.2e1) * t196 + (t25 / 0.2e1 - t12 / 0.2e1) * t195;
t62 = 0.2e1 * t131;
t266 = t62 * qJD(3) + qJD(4) * t263 + t2 * qJD(5);
t243 = Icges(6,5) * t225 + Icges(6,6) * t227;
t234 = (-t245 + t248) * t227 - t371 * t225;
t237 = t195 * t12 / 0.2e1 - t323 - (t11 + t24) * t196 / 0.2e1 + (t234 * t195 + t269 * t225 - t267 * t227 - t361) * t317 + (t234 * t196 + t270 * t225 + t268 * t227 + t25 + t362) * t319;
t236 = t268 * t225 - t270 * t227;
t235 = t267 * t225 + t269 * t227;
t165 = t243 * t195;
t164 = t243 * t196;
t129 = t359 * t360;
t104 = t109 - t307;
t101 = t195 * t360 + t196 * t359;
t72 = qJD(5) * t262;
t71 = qJD(5) * t263;
t63 = 0.2e1 * t118;
t58 = t64 * qJD(5);
t57 = t63 * qJD(5);
t39 = m(6) * t41 + t255;
t38 = m(6) * t40 + t255;
t35 = t311 + t321 + t331;
t34 = t312 + t322 + t332;
t20 = t328 / 0.2e1;
t18 = t329 / 0.2e1;
t14 = t324 * qJD(2);
t13 = t314 + t316 + t325 + t334;
t8 = t20 - t329 / 0.2e1 + t251;
t7 = t18 - t328 / 0.2e1 + t251;
t6 = t264 - t355;
t5 = t264 + t355;
t4 = t351 + t265 - t264;
t3 = t18 + t20 - t251;
t1 = [m(6) * (-t102 + t104) * t103 * t345 / 0.4e1 + t13 * qJD(2) + t34 * qJD(3) + t38 * qJD(5), t13 * qJD(1) + t5 * qJD(3) + t7 * qJD(5) + t324 * t343 / 0.4e1 + ((-t33 + t32) * t339 + t314 / 0.2e1 + t65 * t340 + t316 / 0.2e1) * t344, qJD(1) * t34 + qJD(2) * t5 + t58, t72, t38 * qJD(1) + t7 * qJD(2) + (m(6) * (t102 * t286 - t129 + (t103 * t250 + t291) * t196) + t237) * qJD(5) + t305; t6 * qJD(3) + t8 * qJD(5) - t14 + (-t325 / 0.4e1 - t314 / 0.4e1 - t334 / 0.4e1 - t316 / 0.4e1) * t345 + (-t103 * t107 + t104 * t108 + t32) * t372, -qJD(1) * t324 + qJD(3) * t35 + qJD(5) * t39 - t14, qJD(1) * t6 + qJD(2) * t35 + t58, t72, t8 * qJD(1) + t39 * qJD(2) + (m(6) * (t107 * t286 - t129 + (t108 * t250 + t291) * t196) + t237) * qJD(5) + t305; t4 * qJD(2) + t57 + (-t322 / 0.4e1 - t312 / 0.4e1 - t332 / 0.4e1) * t345 + (-t259 * t104 + t95) * t372, t4 * qJD(1) + t57 + (-t321 / 0.4e1 - t311 / 0.4e1 - t331 / 0.4e1) * t343 + t351 * t344, 0, 0, -(-t195 * t258 + t196 * t259) * t250 * t303 + t352 * t63; t71, t71, 0, 0, -t101 * t303 + t304; (-t255 + (t104 * t359 - t40 - t87) * m(6)) * qJD(1) + t3 * qJD(2) + t266, t3 * qJD(1) + (-t255 + (t21 - t41) * m(6)) * qJD(2) + t266, t352 * t62, t304, (m(6) * (t82 * t101 + (t347 + t348) * t250 * t249) + (t347 * t165 + (t236 * t195 + (-t164 + t235) * t196) * t195) * t317 + (t348 * t164 + (t235 * t196 + (-t165 + t236) * t195) * t196) * t319) * qJD(5) + t352 * t2;];
Cq = t1;
