% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:46
% EndTime: 2019-12-31 16:27:55
% DurationCPUTime: 7.89s
% Computational Cost: add. (5449->352), mult. (6801->458), div. (0->0), fcn. (5325->4), ass. (0->203)
t187 = pkin(6) + qJ(2);
t183 = sin(t187);
t184 = cos(t187);
t188 = sin(qJ(3));
t185 = Icges(5,5) * t188;
t189 = cos(qJ(3));
t225 = Icges(5,1) * t189 + t185;
t90 = -Icges(5,4) * t184 + t183 * t225;
t282 = t183 * t188;
t161 = Icges(4,4) * t282;
t281 = t183 * t189;
t297 = Icges(4,5) * t184;
t92 = Icges(4,1) * t281 - t161 - t297;
t391 = t90 + t92;
t299 = Icges(4,4) * t188;
t144 = Icges(4,2) * t189 + t299;
t223 = Icges(5,3) * t189 - t185;
t390 = t144 + t223;
t394 = t390 * t189;
t296 = Icges(5,5) * t189;
t146 = Icges(5,1) * t188 - t296;
t186 = Icges(4,4) * t189;
t148 = Icges(4,1) * t188 + t186;
t389 = t146 + t148;
t293 = Icges(4,6) * t184;
t88 = Icges(4,4) * t281 - Icges(4,2) * t282 - t293;
t313 = t188 * t88;
t139 = Icges(5,3) * t188 + t296;
t82 = -Icges(5,6) * t184 + t139 * t183;
t386 = t188 * t82 + t391 * t189 - t313;
t290 = Icges(4,3) * t184;
t84 = Icges(4,5) * t281 - Icges(4,6) * t282 - t290;
t393 = t386 * t183 - t184 * t84;
t141 = Icges(4,5) * t189 - Icges(4,6) * t188;
t85 = Icges(4,3) * t183 + t141 * t184;
t143 = Icges(5,4) * t189 + Icges(5,6) * t188;
t87 = Icges(5,2) * t183 + t143 * t184;
t392 = t85 + t87;
t91 = Icges(5,4) * t183 + t184 * t225;
t149 = Icges(4,1) * t189 - t299;
t93 = Icges(4,5) * t183 + t149 * t184;
t384 = t91 + t93;
t257 = qJD(3) * t188;
t388 = t384 * t188;
t140 = Icges(4,5) * t188 + Icges(4,6) * t189;
t142 = Icges(5,4) * t188 - Icges(5,6) * t189;
t387 = t142 + t140;
t377 = -t390 * t188 + t389 * t189;
t385 = qJD(3) * t394 + t389 * t257;
t279 = t184 * t189;
t280 = t184 * t188;
t160 = Icges(5,5) * t279;
t292 = Icges(5,6) * t183;
t83 = Icges(5,3) * t280 + t160 + t292;
t383 = t392 * t183 + t384 * t279 + t83 * t280;
t86 = -Icges(5,2) * t184 + t143 * t183;
t80 = t183 * t86;
t382 = -t183 * t84 - t391 * t279 - t82 * t280 - t80;
t224 = -Icges(4,2) * t188 + t186;
t381 = (t139 - t224) * qJD(3);
t380 = (t149 + t225) * qJD(3);
t314 = t184 * t86;
t367 = -t314 + t393;
t366 = -t280 * t88 - t382;
t89 = Icges(4,6) * t183 + t184 * t224;
t365 = -t280 * t89 + t383;
t283 = t142 * t184;
t285 = t140 * t184;
t376 = t377 * t183 - t283 - t285;
t284 = t142 * t183;
t286 = t140 * t183;
t375 = t377 * t184 + t284 + t286;
t239 = t184 * t87 - t91 * t281 - t83 * t282;
t71 = t93 * t281;
t243 = t184 * t85 - t71;
t33 = -t282 * t89 - t243;
t374 = -t239 + t33;
t373 = t387 * qJD(3);
t312 = t188 * t89;
t372 = t188 * t83 + t384 * t189 - t312;
t363 = rSges(5,3) + qJ(4);
t371 = (-t373 * t183 + (-t386 + t392) * qJD(2)) * t184;
t370 = t380 * t189 + t381 * t188 + (-t389 * t188 - t394) * qJD(3) + t387 * qJD(2);
t368 = (-t141 - t143) * qJD(3) + t377 * qJD(2);
t260 = qJD(2) * t184;
t328 = rSges(5,1) + pkin(3);
t364 = t375 * qJD(2);
t151 = pkin(3) * t188 - qJ(4) * t189;
t152 = rSges(5,1) * t188 - rSges(5,3) * t189;
t266 = t151 + t152;
t154 = pkin(3) * t189 + qJ(4) * t188;
t155 = rSges(5,1) * t189 + rSges(5,3) * t188;
t362 = t154 + t155;
t359 = (t365 * t183 - t366 * t184) * qJD(3);
t358 = (t374 * t183 - t367 * t184) * qJD(3);
t357 = t376 * qJD(2);
t355 = -t373 * t184 + (-t141 * t183 + t290 - t372 - t86) * qJD(2);
t354 = 0.2e1 * qJD(3);
t353 = t357 + t358;
t352 = t359 + t364;
t351 = -t386 * qJD(3) + t385 * t183 + ((t139 * t184 + t292 - t89) * t189 - t388) * qJD(2);
t350 = t372 * qJD(3) - t385 * t184 + ((-t183 * t224 + t293 + t82) * t189 + (-t149 * t183 + t297 - t90) * t188) * qJD(2);
t349 = -t183 * t368 + t184 * t370;
t348 = t183 * t370 + t184 * t368;
t347 = (-t82 + t88) * t189 + t391 * t188;
t346 = (-t83 + t89) * t189 + t388;
t255 = qJD(4) * t189;
t177 = t183 * rSges(5,2);
t303 = t279 * t328 + t363 * t280 + t177;
t179 = t184 * rSges(5,2);
t304 = t183 * t362 - t179;
t25 = -t255 + qJD(1) + (t183 * t304 + t184 * t303) * qJD(3);
t344 = qJD(3) * t25;
t240 = t266 * qJD(3);
t256 = qJD(4) * t188;
t207 = -t240 + t256;
t343 = t207 * t183;
t126 = t184 * pkin(2) + t183 * pkin(5);
t342 = t126 + t303;
t341 = t314 + t383;
t150 = t184 * t256;
t258 = qJD(3) * t184;
t251 = t189 * t258;
t340 = rSges(5,2) * t260 + t251 * t363 + t150;
t306 = -t148 * t183 - t88;
t310 = -Icges(4,2) * t281 - t161 + t92;
t333 = -t188 * t310 + t189 * t306;
t308 = -t146 * t183 + t82;
t321 = -t223 * t183 + t90;
t332 = -t188 * t321 + t189 * t308;
t111 = t152 * t183;
t324 = t154 * t260 + (-qJD(3) * t151 + t256) * t183 - qJD(3) * t111 + (t155 * t184 + t177) * qJD(2);
t261 = qJD(2) * t183;
t206 = -t184 * t257 - t189 * t261;
t254 = t188 * t261;
t325 = t206 * t328 - t363 * t254 + t340;
t1 = (t256 + t325 * t184 + t324 * t183 + (-t183 * t303 + t184 * t304) * qJD(2)) * qJD(3);
t331 = m(5) * t1;
t326 = qJD(2) / 0.2e1;
t320 = -t223 * t184 + t91;
t319 = rSges(4,1) * t189;
t153 = rSges(4,1) * t188 + rSges(4,2) * t189;
t116 = t153 * t184;
t259 = qJD(3) * t183;
t253 = t153 * t259;
t176 = t183 * rSges(4,3);
t97 = rSges(4,1) * t279 - rSges(4,2) * t280 + t176;
t302 = t126 + t97;
t40 = qJD(2) * t302 - t253;
t318 = t116 * t40;
t181 = t184 * pkin(5);
t125 = pkin(2) * t183 - t181;
t252 = t153 * t258;
t264 = rSges(4,2) * t282 + t184 * rSges(4,3);
t95 = rSges(4,1) * t281 - t264;
t39 = -t252 + (-t125 - t95) * qJD(2);
t317 = t183 * t39;
t316 = t184 * t39;
t309 = -t144 * t184 + t93;
t307 = -Icges(5,1) * t280 + t160 + t83;
t305 = -t148 * t184 - t89;
t276 = -t151 * t183 - t111;
t275 = t266 * t184;
t274 = -qJD(3) * t362 + t255;
t271 = rSges(4,2) * t254 + rSges(4,3) * t260;
t270 = -t223 + t225;
t269 = t139 - t146;
t268 = -t144 + t149;
t267 = t148 + t224;
t263 = qJD(2) * t141;
t262 = qJD(2) * t143;
t250 = -pkin(2) - t319;
t247 = -t259 / 0.2e1;
t244 = t258 / 0.2e1;
t242 = -t84 + t312;
t236 = -t189 * t328 - pkin(2);
t234 = -rSges(4,2) * t188 + t319;
t233 = -t183 * t40 - t316;
t232 = t183 * t95 + t184 * t97;
t220 = -t184 * t240 + t150;
t219 = qJD(3) * (t255 + t274);
t112 = t153 * t183;
t205 = -t188 * t320 + t189 * t307;
t204 = -t188 * t309 + t189 * t305;
t203 = (t188 * t269 + t189 * t270) * qJD(2);
t202 = (-t188 * t267 + t189 * t268) * qJD(2);
t66 = rSges(4,1) * t206 - rSges(4,2) * t251 + t271;
t68 = -qJD(3) * t112 + (t184 * t234 + t176) * qJD(2);
t201 = t183 * t68 + t184 * t66 + (-t183 * t97 + t184 * t95) * qJD(2);
t169 = pkin(5) * t260;
t134 = t234 * qJD(3);
t122 = qJD(2) * t125;
t121 = t126 * qJD(2);
t118 = qJD(2) * (-pkin(2) * t261 + t169);
t38 = qJD(3) * t232 + qJD(1);
t29 = qJD(2) * t342 + t343;
t28 = (-t125 - t304) * qJD(2) + t220;
t27 = -t134 * t258 + (-t121 - t68 + t253) * qJD(2);
t26 = -t134 * t259 + t118 + (t66 - t252) * qJD(2);
t16 = t201 * qJD(3);
t15 = t184 * t219 + (-t121 - t324 - t343) * qJD(2);
t14 = t118 + t183 * t219 + (t184 * t207 + t325) * qJD(2);
t2 = [m(4) * t16 + t331; (((t33 - t71 + (t85 + t313) * t184 + t382) * t184 + (t341 + t367 - t393) * t183) * qJD(3) + t364) * t244 + (t377 * qJD(3) + t380 * t188 - t381 * t189) * qJD(2) + (t15 * (t179 + t181) + t14 * t342 + t29 * (t169 + t340) + (-t29 * t328 * t257 + t28 * (-t188 * t363 + t236) * qJD(2)) * t184 + (t15 * t236 + (-t28 * qJD(4) - t15 * t363) * t188 + t28 * (t188 * t328 - t189 * t363) * qJD(3) + (t28 * (-rSges(5,2) - pkin(5)) + t29 * (-pkin(2) - t362)) * qJD(2)) * t183 - (-qJD(2) * t304 - t122 + t220 - t28) * t29) * m(5) + (t27 * (t183 * t250 + t181 + t264) + t26 * t302 + t40 * (t169 + t271) + (t153 * t317 - t318) * qJD(3) + ((-pkin(2) - t234) * t316 + (t39 * (-rSges(4,3) - pkin(5)) + t40 * t250) * t183) * qJD(2) - (-qJD(2) * t95 - t122 - t252 - t39) * t40) * m(4) + (((t184 * t242 - t341 + t365) * t184 + (t183 * t242 + t239 + t243 + t366 - t80) * t183) * qJD(3) + t353 - t357) * t247 + (t349 + t350) * t259 / 0.2e1 - (t348 - t351 + t352) * t258 / 0.2e1 + ((t347 + t376) * t183 + (t346 + t375) * t184) * qJD(3) * t326; -((((-t310 - t321) * t184 + (t309 + t320) * t183) * t189 + ((-t306 - t308) * t184 + (t305 + t307) * t183) * t188) * qJD(3) + ((t267 - t269) * t189 + (t268 + t270) * t188) * qJD(2)) * qJD(2) / 0.2e1 + (t351 * t184 + t350 * t183 + (t183 * t347 + t184 * t346) * qJD(2)) * t326 + ((-t259 * t283 + t262) * t183 + (t203 + (-t332 * t184 + (t284 + t205) * t183) * qJD(3)) * t184 + (-t259 * t285 + t263) * t183 + (t202 + (-t333 * t184 + (t286 + t204) * t183) * qJD(3)) * t184) * t247 + ((-t258 * t284 - t262) * t184 + (t203 + (t205 * t183 + (t283 - t332) * t184) * qJD(3)) * t183 + (-t258 * t286 - t263) * t184 + (t202 + (t204 * t183 + (t285 - t333) * t184) * qJD(3)) * t183) * t244 + (-(t188 * t25 + (t183 * t29 + t184 * t28) * t189) * qJD(4) - (-t275 * t29 - t276 * t28) * qJD(2) - ((-t25 * t275 - t28 * t362) * t184 + (t25 * t276 - t29 * t362) * t183) * qJD(3) + (-t15 * t266 + t28 * t274 + t1 * t303 + t25 * t325 + (t25 * t304 - t266 * t29) * qJD(2)) * t184 + (-t14 * t266 + t29 * t274 + t1 * t304 + t25 * t324 + (-t25 * t303 + t266 * t28) * qJD(2)) * t183) * m(5) + (t16 * t232 + t38 * t201 + t233 * t134 + (-t26 * t183 - t27 * t184 + (-t184 * t40 + t317) * qJD(2)) * t153 - (t112 * t39 - t318) * qJD(2) - (t38 * (-t112 * t183 - t116 * t184) + t233 * t234) * qJD(3)) * m(4) + (t349 * qJD(2) + (t365 * t260 + (t366 * qJD(2) + t355 * t183 - t371) * t183) * t354) * t183 / 0.2e1 - (t348 * qJD(2) + ((t374 * qJD(2) + t371) * t184 + (t367 * qJD(2) - t184 * t355) * t183) * t354) * t184 / 0.2e1 + (t353 + t358) * t261 / 0.2e1 + (t352 + t359) * t260 / 0.2e1; -t189 * t331 + 0.2e1 * (m(5) * (t14 * t183 + t15 * t184 + t344) / 0.2e1 - m(5) * (t183 ^ 2 + t184 ^ 2) * t344 / 0.2e1) * t188;];
tauc = t2(:);
