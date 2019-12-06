% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:05
% EndTime: 2019-12-05 18:20:13
% DurationCPUTime: 3.79s
% Computational Cost: add. (23412->246), mult. (17995->350), div. (0->0), fcn. (18859->10), ass. (0->155)
t232 = qJ(1) + qJ(2);
t229 = pkin(8) + t232;
t227 = sin(t229);
t228 = cos(t229);
t230 = sin(t232);
t309 = pkin(2) * t230;
t251 = -t228 * qJ(4) + t309;
t237 = cos(qJ(5));
t234 = cos(pkin(9));
t235 = sin(qJ(5));
t280 = t234 * t235;
t201 = t227 * t280 + t228 * t237;
t279 = t234 * t237;
t202 = t227 * t279 - t228 * t235;
t269 = t202 * rSges(6,1) - t201 * rSges(6,2);
t233 = sin(pkin(9));
t342 = pkin(4) * t234 + pkin(3) + (rSges(6,3) + pkin(7)) * t233;
t126 = t342 * t227 + t251 + t269;
t307 = sin(qJ(1)) * pkin(1);
t121 = t126 + t307;
t113 = t228 * t121;
t252 = rSges(5,1) * t234 + pkin(3);
t287 = t227 * t233;
t162 = -rSges(5,2) * t287 - t228 * rSges(5,3) + t252 * t227 + t251;
t159 = t162 + t307;
t144 = t228 * t159;
t284 = t228 * t233;
t231 = cos(t232);
t308 = pkin(2) * t231;
t163 = -t308 + rSges(5,2) * t284 - t252 * t228 + (-rSges(5,3) - qJ(4)) * t227;
t306 = cos(qJ(1)) * pkin(1);
t160 = t163 - t306;
t285 = t228 * t162;
t286 = t228 * t126;
t335 = m(6) / 0.2e1;
t203 = -t227 * t237 + t228 * t280;
t204 = t227 * t235 + t228 * t279;
t345 = -t204 * rSges(6,1) + t203 * rSges(6,2);
t349 = -t227 * qJ(4) - t342 * t228 - t308 + t345;
t353 = m(5) / 0.2e1;
t355 = t349 - t306;
t303 = (-t113 - t286 + (-t355 - t349) * t227) * t335 + (-t144 - t285 + (-t160 - t163) * t227) * t353;
t278 = t227 * t349 + t286;
t65 = -t227 * t355 - t113;
t304 = (t65 + t278) * t335 + (-t144 + t285 + (-t160 + t163) * t227) * t353;
t6 = t304 - t303;
t364 = t6 * qJD(1);
t172 = t201 * rSges(6,1) + t202 * rSges(6,2);
t173 = -t203 * rSges(6,1) - t204 * rSges(6,2);
t54 = -t126 * t172 - t173 * t349;
t363 = m(6) * t54;
t51 = -t121 * t172 - t173 * t355;
t362 = m(6) * t51;
t43 = -t121 * t349 + t126 * t355;
t155 = rSges(6,3) * t284 - t345;
t241 = t233 * (-t234 * rSges(6,3) + (rSges(6,1) * t237 - rSges(6,2) * t235) * t233);
t360 = t234 * t155 + t228 * t241;
t294 = Icges(6,4) * t235;
t209 = -Icges(6,5) * t234 + (Icges(6,1) * t237 - t294) * t233;
t266 = t209 + (-Icges(6,2) * t237 - t294) * t233;
t356 = t266 * t235;
t190 = Icges(6,4) * t204;
t149 = -Icges(6,2) * t203 + Icges(6,6) * t284 + t190;
t189 = Icges(6,4) * t203;
t153 = -Icges(6,1) * t204 - Icges(6,5) * t284 + t189;
t275 = t201 * t149 + t202 * t153;
t248 = -t203 * t149 - t204 * t153;
t332 = m(5) * (-t163 * t159 + t160 * t162);
t317 = m(3) * (-t306 * (t230 * rSges(3,1) + t231 * rSges(3,2)) - (-t231 * rSges(3,1) + t230 * rSges(3,2)) * t307);
t315 = m(4) * (-t306 * (t227 * rSges(4,1) + t228 * rSges(4,2) + t309) - (-t228 * rSges(4,1) + t227 * rSges(4,2) - t308) * t307);
t207 = -Icges(6,3) * t234 + (Icges(6,5) * t237 - Icges(6,6) * t235) * t233;
t293 = Icges(6,4) * t237;
t208 = -Icges(6,6) * t234 + (-Icges(6,2) * t235 + t293) * t233;
t351 = (-t203 * t208 + t204 * t209 + t207 * t284) * t234;
t265 = qJD(1) + qJD(2);
t212 = (-Icges(6,5) * t235 - Icges(6,6) * t237) * t233;
t281 = t234 * t212;
t347 = -t281 / 0.2e1 - t233 * t356 / 0.2e1;
t145 = -Icges(6,5) * t202 + Icges(6,6) * t201 - Icges(6,3) * t287;
t295 = Icges(6,4) * t202;
t148 = Icges(6,2) * t201 - Icges(6,6) * t287 - t295;
t188 = Icges(6,4) * t201;
t151 = -Icges(6,1) * t202 - Icges(6,5) * t287 + t188;
t343 = t145 * t284 - t203 * t148 + t204 * t151;
t341 = 0.4e1 * qJD(1);
t154 = -rSges(6,3) * t287 - t269;
t130 = t234 * t154 - t227 * t241;
t330 = m(6) * ((-t121 + t126) * t360 + (-t349 + t355) * t130);
t327 = m(6) * (t54 + t51);
t324 = m(6) * t43;
t321 = m(6) * t65;
t320 = m(6) * t278;
t319 = m(6) * (-t130 * t228 - t227 * t360);
t313 = m(5) * (-t160 * t227 - t144);
t312 = m(5) * (-t163 * t227 - t285);
t310 = m(6) * (t227 * t172 - t228 * t173);
t112 = t310 / 0.2e1;
t73 = t319 / 0.2e1;
t25 = t112 + t73;
t305 = t25 * qJD(4);
t299 = -t130 * t172 - t173 * t360;
t297 = m(6) * qJD(5);
t292 = (t201 * t208 - t202 * t209 - t207 * t287) * t234;
t214 = (-Icges(6,1) * t235 - t293) * t233;
t267 = t208 - t214;
t291 = (-t281 + (-t267 * t237 - t356) * t233) * t234;
t282 = t233 * t237;
t274 = -t201 * t148 + t202 * t151;
t273 = -Icges(6,1) * t201 + t148 - t295;
t272 = Icges(6,1) * t203 + t149 + t190;
t271 = Icges(6,2) * t202 + t151 + t188;
t270 = -Icges(6,2) * t204 - t153 - t189;
t246 = t145 * t287 + t274;
t146 = Icges(6,5) * t204 - Icges(6,6) * t203 + Icges(6,3) * t284;
t64 = t146 * t284 + t248;
t11 = -t292 + (t275 * t228 + (t246 + t248 - t64) * t227) * t233;
t62 = -t146 * t287 + t275;
t12 = t351 + (-(-t275 - t343) * t227 + t246 * t228 + (-t248 - t274) * t228 - t227 * t62 + (-t146 * t227 ^ 2 + (-t145 * t227 - t146 * t228) * t228) * t233) * t233;
t33 = -t292 + (t227 * t246 + t228 * t62) * t233;
t34 = -t351 + (-t227 * t343 + t228 * t64) * t233;
t2 = ((-t33 / 0.2e1 + t11 / 0.2e1) * t228 + (-t12 / 0.2e1 - t34 / 0.2e1) * t227) * t233;
t26 = t73 - t310 / 0.2e1;
t264 = t26 * qJD(4) + t2 * qJD(5);
t259 = -t287 / 0.2e1;
t257 = t284 / 0.2e1;
t254 = -t282 / 0.2e1;
t253 = t282 / 0.2e1;
t249 = t208 * t254 + t214 * t253 + t347;
t245 = t327 / 0.2e1 + t249;
t242 = t208 * t253 + t214 * t254 - t347;
t166 = Icges(6,5) * t201 + Icges(6,6) * t202;
t55 = -t234 * t166 + (-t271 * t235 - t273 * t237) * t233;
t167 = -Icges(6,5) * t203 - Icges(6,6) * t204;
t56 = -t234 * t167 + (-t270 * t235 - t272 * t237) * t233;
t85 = t266 * t201 + t267 * t202 - t212 * t287;
t86 = -t266 * t203 - t267 * t204 + t212 * t284;
t240 = -t11 * t284 / 0.2e1 - t291 + (t55 + t85) * t259 + (t12 + t34) * t287 / 0.2e1 + (t33 + t56 + t86) * t257;
t215 = (-rSges(6,1) * t235 - rSges(6,2) * t237) * t233;
t141 = t234 * t173 + t215 * t284;
t140 = -t234 * t172 + t215 * t287;
t108 = (-t172 * t228 - t173 * t227) * t233;
t50 = t312 - t320;
t49 = t313 + t321;
t42 = t249 + t363;
t40 = t249 + t362;
t24 = t112 - t319 / 0.2e1;
t21 = t25 * qJD(5);
t20 = t24 * qJD(5);
t14 = t330 / 0.2e1;
t13 = t315 + t317 + t324 + t332;
t8 = t303 + t304;
t5 = -t330 / 0.2e1 + t245;
t4 = t14 + t245;
t3 = t14 - t327 / 0.2e1 + t242;
t1 = [t13 * qJD(2) + t49 * qJD(4) + t40 * qJD(5), t13 * qJD(1) + t8 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t317 / 0.2e1 + t315 / 0.2e1 + t332 / 0.2e1 + t43 * t335) * qJD(2), 0, t49 * qJD(1) + t8 * qJD(2) + t21, t40 * qJD(1) + t4 * qJD(2) + (m(6) * (-t140 * t121 + t141 * t355 + t299) + t240) * qJD(5) + t305; -t6 * qJD(4) + t5 * qJD(5) + (-t317 / 0.4e1 - t315 / 0.4e1 - t332 / 0.4e1 - t324 / 0.4e1) * t341, t50 * qJD(4) + t42 * qJD(5), 0, t50 * qJD(2) + t21 - t364, t5 * qJD(1) + t42 * qJD(2) + (m(6) * (-t140 * t126 + t141 * t349 + t299) + t240) * qJD(5) + t305; 0, 0, 0, 0, t108 * t297; t6 * qJD(2) + t20 + (-t313 / 0.4e1 - t321 / 0.4e1) * t341, t364 + t20 + 0.4e1 * (t320 / 0.4e1 - t312 / 0.4e1) * qJD(2), 0, 0, (t140 * t227 + t141 * t228) * t297 + t265 * t24; (t242 - t362) * qJD(1) + t3 * qJD(2) + t264, t3 * qJD(1) + (t242 - t363) * qJD(2) + t264, 0, t265 * t26, (m(6) * ((-t154 * t228 - t155 * t227) * t233 * t108 - t130 * t140 + t141 * t360) - t234 * (-t291 + (-t227 * t55 + t228 * t56) * t233) / 0.2e1 + (-t85 * t234 - (-t166 * t287 + t271 * t201 + t273 * t202) * t287 + (-t167 * t287 + t270 * t201 + t272 * t202) * t284) * t259 + (-t86 * t234 - (t166 * t284 - t271 * t203 - t273 * t204) * t287 + (t167 * t284 - t270 * t203 - t272 * t204) * t284) * t257) * qJD(5) + t265 * t2;];
Cq = t1;
