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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:35:51
% DurationCPUTime: 6.49s
% Computational Cost: add. (13270->287), mult. (15985->410), div. (0->0), fcn. (17069->8), ass. (0->164)
t228 = qJ(1) + pkin(7);
t225 = sin(t228);
t226 = cos(t228);
t233 = cos(qJ(4));
t230 = cos(pkin(8));
t231 = sin(qJ(4));
t300 = t230 * t231;
t200 = -t225 * t233 + t226 * t300;
t299 = t230 * t233;
t304 = t225 * t231;
t201 = t226 * t299 + t304;
t396 = rSges(6,1) + pkin(4);
t116 = rSges(6,2) * t201 + t200 * t396;
t198 = t225 * t300 + t226 * t233;
t302 = t226 * t231;
t199 = t225 * t299 - t302;
t161 = rSges(5,1) * t198 + rSges(5,2) * t199;
t163 = -rSges(5,1) * t200 - rSges(5,2) * t201;
t282 = -t199 * rSges(6,2) - t198 * t396;
t345 = m(6) / 0.2e1;
t346 = m(5) / 0.2e1;
t318 = (t116 * t226 - t225 * t282) * t345 + (t161 * t225 - t163 * t226) * t346;
t229 = sin(pkin(8));
t305 = t225 * t229;
t322 = -qJ(5) - pkin(6);
t361 = -rSges(6,1) * t199 + rSges(6,2) * t198 + pkin(4) * t302 + t305 * t322;
t284 = (-pkin(4) * t299 + pkin(6) * t229) * t225 - rSges(6,3) * t305 + t361;
t325 = pkin(4) * t233;
t379 = ((rSges(6,1) * t233 - rSges(6,2) * t231 + t325) * t229 + (-rSges(6,3) - qJ(5)) * t230) * t229;
t65 = -t225 * t379 + t230 * t284;
t303 = t226 * t229;
t324 = t230 * pkin(3);
t257 = t230 * (pkin(3) + t325);
t402 = -t201 * rSges(6,1) + t200 * rSges(6,2) - pkin(4) * t304 - t226 * t257;
t283 = (-qJ(5) * t229 + t324) * t226 - rSges(6,3) * t303 + t402;
t67 = t226 * t379 - t230 * t283;
t278 = rSges(5,1) * t199 - rSges(5,2) * t198;
t140 = -rSges(5,3) * t305 - t278;
t238 = t229 * (-rSges(5,3) * t230 + (rSges(5,1) * t233 - rSges(5,2) * t231) * t229);
t98 = t140 * t230 - t225 * t238;
t362 = -rSges(5,1) * t201 + rSges(5,2) * t200;
t142 = rSges(5,3) * t303 - t362;
t99 = t142 * t230 + t226 * t238;
t320 = (-t225 * t67 - t226 * t65) * t345 + (-t225 * t99 - t226 * t98) * t346;
t4 = t320 - t318;
t414 = t4 * qJD(1);
t411 = Icges(5,5) + Icges(6,5);
t410 = Icges(5,6) + Icges(6,6);
t409 = -Icges(5,3) - Icges(6,3);
t397 = t200 * t410 - t201 * t411 + t303 * t409;
t413 = t397 * t303;
t204 = (-Icges(6,5) * t231 - Icges(6,6) * t233) * t229;
t205 = (-Icges(5,5) * t231 - Icges(5,6) * t233) * t229;
t310 = Icges(5,4) * t231;
t197 = -Icges(5,5) * t230 + (Icges(5,1) * t233 - t310) * t229;
t207 = (-Icges(5,2) * t233 - t310) * t229;
t272 = t197 + t207;
t307 = Icges(6,4) * t231;
t196 = -Icges(6,5) * t230 + (Icges(6,1) * t233 - t307) * t229;
t206 = (-Icges(6,2) * t233 - t307) * t229;
t273 = t196 + t206;
t309 = Icges(5,4) * t233;
t195 = -Icges(5,6) * t230 + (-Icges(5,2) * t231 + t309) * t229;
t209 = (-Icges(5,1) * t231 - t309) * t229;
t274 = t195 - t209;
t306 = Icges(6,4) * t233;
t194 = -Icges(6,6) * t230 + (-Icges(6,2) * t231 + t306) * t229;
t208 = (-Icges(6,1) * t231 - t306) * t229;
t275 = t194 - t208;
t412 = ((t205 + t204) * t230 + ((t274 + t275) * t233 + (t272 + t273) * t231) * t229) * t230;
t176 = Icges(6,4) * t201;
t126 = -Icges(6,2) * t200 + Icges(6,6) * t303 + t176;
t179 = Icges(5,4) * t201;
t129 = -Icges(5,2) * t200 + Icges(5,6) * t303 + t179;
t175 = Icges(6,4) * t200;
t133 = -Icges(6,1) * t201 - Icges(6,5) * t303 + t175;
t178 = Icges(5,4) * t200;
t136 = -Icges(5,1) * t201 - Icges(5,5) * t303 + t178;
t391 = (t129 + t126) * t198 - (-t136 - t133) * t199;
t408 = -(-(t197 / 0.2e1 + t207 / 0.2e1 + t196 / 0.2e1 + t206 / 0.2e1) * t231 + (t209 / 0.2e1 - t195 / 0.2e1 + t208 / 0.2e1 - t194 / 0.2e1) * t233) * t229 + (t205 / 0.2e1 + t204 / 0.2e1) * t230;
t368 = m(6) * t229;
t258 = sin(qJ(1)) * pkin(1) - t226 * qJ(3);
t88 = (rSges(6,3) * t229 + pkin(2) + t257) * t225 + t258 - t361;
t323 = cos(qJ(1)) * pkin(1);
t253 = -t225 * qJ(3) - t323;
t89 = (-pkin(2) + (-rSges(6,3) + t322) * t229) * t226 + t253 + t402;
t401 = (t225 * t88 - t226 * t89) * t368;
t398 = t198 * t410 - t199 * t411 + t305 * t409;
t395 = qJD(1) * t401;
t308 = Icges(6,4) * t199;
t125 = Icges(6,2) * t198 - Icges(6,6) * t305 - t308;
t311 = Icges(5,4) * t199;
t128 = Icges(5,2) * t198 - Icges(5,6) * t305 - t311;
t174 = Icges(6,4) * t198;
t131 = -Icges(6,1) * t199 - Icges(6,5) * t305 + t174;
t177 = Icges(5,4) * t198;
t134 = -Icges(5,1) * t199 - Icges(5,5) * t305 + t177;
t392 = (t131 + t134) * t199 + (-t125 - t128) * t198;
t390 = t116 * t225 + t226 * t282;
t387 = t305 * t397 + t391;
t354 = (rSges(5,3) + pkin(6)) * t229 + pkin(2) + t324;
t97 = -t226 * t354 + t253 + t362;
t378 = t305 * t398 + t392;
t376 = t229 * t397;
t148 = Icges(6,5) * t198 + Icges(6,6) * t199;
t149 = -Icges(6,5) * t200 - Icges(6,6) * t201;
t150 = Icges(5,5) * t198 + Icges(5,6) * t199;
t151 = -Icges(5,5) * t200 - Icges(5,6) * t201;
t375 = -(t151 + t149) * t303 + (t150 + t148) * t305;
t373 = -t229 / 0.2e1;
t371 = -t230 / 0.2e1;
t285 = -Icges(5,2) * t201 - t136 - t178;
t287 = -Icges(6,2) * t201 - t133 - t175;
t360 = t285 + t287;
t286 = Icges(5,2) * t199 + t134 + t177;
t288 = Icges(6,2) * t199 + t131 + t174;
t359 = t286 + t288;
t289 = Icges(5,1) * t200 + t129 + t179;
t291 = Icges(6,1) * t200 + t126 + t176;
t358 = t289 + t291;
t290 = -Icges(5,1) * t198 + t128 - t311;
t292 = -Icges(6,1) * t198 + t125 - t308;
t357 = t290 + t292;
t352 = t225 ^ 2;
t349 = 0.4e1 * qJD(1);
t348 = 2 * qJD(4);
t344 = m(4) * (-(-t323 + (-rSges(4,3) - qJ(3)) * t225 + rSges(4,2) * t303) * t225 - t226 * (-rSges(4,2) * t305 - t226 * rSges(4,3) + t258));
t96 = t225 * t354 + t258 + t278;
t341 = m(5) * (-t161 * t96 - t163 * t97);
t340 = m(5) * (-t225 * t97 - t226 * t96);
t334 = (t225 * t65 - t226 * t67) * t368;
t332 = m(6) * (t116 * t89 + t282 * t88);
t331 = m(6) * (-t225 * t89 - t226 * t88);
t330 = t390 * t368;
t146 = (t226 ^ 2 + t352) * t368;
t269 = t146 * qJD(1);
t250 = (t225 * t378 + t226 * t387) * t373 + (t391 * t226 + (t378 + t413) * t225) * t229 / 0.2e1;
t247 = (t352 * t376 + (t226 * t376 + t378 - t392 - t413) * t226 + (-t303 * t398 - t387 + t391) * t225) * t373 + (t371 + t230 / 0.2e1) * ((t409 * t230 + (-t231 * t410 + t233 * t411) * t229) * t303 + (t197 + t196) * t201 + (-t195 - t194) * t200);
t246 = (-rSges(6,2) * t233 - t231 * t396) * t229 ^ 2;
t211 = (-rSges(5,1) * t231 - rSges(5,2) * t233) * t229;
t114 = t163 * t230 + t211 * t303;
t113 = -t161 * t230 + t211 * t305;
t86 = (-t161 * t226 - t163 * t225) * t229;
t85 = -t116 * t230 + t226 * t246;
t84 = t225 * t246 + t230 * t282;
t69 = -t330 / 0.2e1;
t68 = t390 * t229;
t60 = -t200 * t272 - t201 * t274 + t205 * t303;
t59 = -t200 * t273 - t201 * t275 + t204 * t303;
t58 = t198 * t272 + t199 * t274 - t205 * t305;
t57 = t198 * t273 + t199 * t275 - t204 * t305;
t38 = -t151 * t230 + (-t231 * t285 - t233 * t289) * t229;
t37 = -t150 * t230 + (-t231 * t286 - t233 * t290) * t229;
t36 = -t149 * t230 + (-t231 * t287 - t233 * t291) * t229;
t35 = -t148 * t230 + (-t231 * t288 - t233 * t292) * t229;
t26 = t334 / 0.2e1;
t23 = t331 + t340 + t344;
t15 = t332 + t341 - t408;
t7 = t26 + t69;
t6 = t26 + t330 / 0.2e1;
t5 = t69 - t334 / 0.2e1;
t3 = t318 + t320;
t1 = (t225 * t247 + t226 * t250) * t229;
t2 = [qJD(3) * t23 + qJD(4) * t15 + qJD(5) * t401, 0, qJD(1) * t23 + qJD(4) * t3, t15 * qJD(1) + t3 * qJD(3) + t7 * qJD(5) + ((-t113 * t96 + t114 * t97 - t161 * t98 - t163 * t99) * t346 + (t116 * t67 + t282 * t65 - t84 * t88 + t85 * t89) * t345) * t348 + (t412 + ((t38 / 0.2e1 + t36 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1 - t250) * t226 + (-t37 / 0.2e1 - t35 / 0.2e1 - t58 / 0.2e1 - t57 / 0.2e1 - t247) * t225) * t229) * qJD(4), qJD(4) * t7 + t395; 0, 0, 0, (t345 * t68 + t346 * t86) * t348, 0; -t4 * qJD(4) - t146 * qJD(5) + (-t344 / 0.4e1 - t340 / 0.4e1 - t331 / 0.4e1) * t349, 0, 0, -t414 + ((t113 * t225 + t114 * t226) * t346 + (t225 * t84 + t226 * t85) * t345) * t348, -t269; t4 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t341 / 0.4e1 - t332 / 0.4e1) * t349 + t408 * qJD(1), 0, t414, t1 * qJD(1) + ((t412 + ((t36 + t38) * t226 + (-t35 - t37) * t225) * t229) * t371 - ((t198 * t360 + t199 * t358) * t303 + (-t58 - t57) * t230 + (-t359 * t198 - t357 * t199 + t375) * t305) * t305 / 0.2e1 + ((t200 * t359 + t201 * t357) * t305 + (-t60 - t59) * t230 + (-t360 * t200 - t358 * t201 - t375) * t303) * t303 / 0.2e1 + m(5) * (-t113 * t98 + t114 * t99 + (-t140 * t226 - t142 * t225) * t229 * t86) + m(6) * ((t225 * t283 - t226 * t284) * t229 * t68 - t65 * t84 + t67 * t85)) * qJD(4), t6 * qJD(1); t146 * qJD(3) + t5 * qJD(4) - t395, 0, t269, t5 * qJD(1) + m(6) * (-t230 * t68 + (-t225 * t85 + t226 * t84) * t229) * qJD(4), 0;];
Cq = t2;
