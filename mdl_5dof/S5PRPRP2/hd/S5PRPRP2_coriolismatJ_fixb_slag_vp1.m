% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:32
% EndTime: 2019-12-05 15:30:47
% DurationCPUTime: 6.01s
% Computational Cost: add. (13186->283), mult. (15893->403), div. (0->0), fcn. (16983->6), ass. (0->162)
t228 = pkin(7) + qJ(2);
t225 = sin(t228);
t226 = cos(t228);
t233 = cos(qJ(4));
t230 = cos(pkin(8));
t232 = sin(qJ(4));
t288 = t230 * t232;
t204 = t225 * t288 + t226 * t233;
t287 = t230 * t233;
t289 = t226 * t232;
t205 = t225 * t287 - t289;
t398 = rSges(6,1) + pkin(4);
t121 = t205 * rSges(6,2) + t398 * t204;
t171 = -rSges(5,1) * t204 - rSges(5,2) * t205;
t292 = t225 * t232;
t207 = t226 * t287 + t292;
t369 = -t225 * t233 + t226 * t288;
t173 = -rSges(5,1) * t369 - rSges(5,2) * t207;
t271 = t207 * rSges(6,2) + t398 * t369;
t335 = m(6) / 0.2e1;
t337 = m(5) / 0.2e1;
t307 = (t121 * t225 + t226 * t271) * t335 + (-t171 * t225 - t173 * t226) * t337;
t247 = rSges(5,1) * t207 - rSges(5,2) * t369;
t229 = sin(pkin(8));
t290 = t226 * t229;
t154 = rSges(5,3) * t290 + t247;
t203 = -t230 * rSges(5,3) + (rSges(5,1) * t233 - rSges(5,2) * t232) * t229;
t103 = t154 * t230 + t203 * t290;
t293 = t225 * t229;
t354 = -t205 * rSges(5,1) + t204 * rSges(5,2);
t150 = rSges(5,3) * t293 - t354;
t387 = t150 * t230 + t203 * t293;
t231 = -qJ(5) - pkin(6);
t311 = pkin(6) + t231;
t314 = pkin(3) * t230;
t224 = pkin(4) * t233 + pkin(3);
t294 = t224 * t230;
t395 = -t205 * rSges(6,1) + t204 * rSges(6,2) + pkin(4) * t289 - t225 * t294;
t274 = -(t311 * t229 + t314) * t225 + rSges(6,3) * t293 - t395;
t312 = -pkin(3) + t224;
t372 = (t311 - rSges(6,3)) * t230 + (rSges(6,1) * t233 - rSges(6,2) * t232 + t312) * t229;
t394 = t274 * t230 + t372 * t293;
t353 = rSges(6,1) * t207 - rSges(6,2) * t369 - t231 * t290;
t273 = pkin(4) * t292 + (-pkin(6) * t229 + t312 * t230) * t226 + rSges(6,3) * t290 + t353;
t69 = t273 * t230 + t372 * t290;
t309 = (-t69 * t225 + t394 * t226) * t335 + (-t103 * t225 + t387 * t226) * t337;
t4 = t309 - t307;
t403 = t4 * qJD(2);
t391 = Icges(5,5) + Icges(6,5);
t390 = -Icges(5,6) - Icges(6,6);
t389 = Icges(6,3) + Icges(5,3);
t208 = (-Icges(6,5) * t232 - Icges(6,6) * t233) * t229;
t209 = (-Icges(5,5) * t232 - Icges(5,6) * t233) * t229;
t299 = Icges(5,4) * t232;
t201 = -Icges(5,5) * t230 + (Icges(5,1) * t233 - t299) * t229;
t211 = (-Icges(5,2) * t233 - t299) * t229;
t264 = t201 + t211;
t296 = Icges(6,4) * t232;
t200 = -Icges(6,5) * t230 + (Icges(6,1) * t233 - t296) * t229;
t210 = (-Icges(6,2) * t233 - t296) * t229;
t265 = t200 + t210;
t298 = Icges(5,4) * t233;
t199 = -Icges(5,6) * t230 + (-Icges(5,2) * t232 + t298) * t229;
t213 = (-Icges(5,1) * t232 - t298) * t229;
t266 = -t199 + t213;
t295 = Icges(6,4) * t233;
t198 = -Icges(6,6) * t230 + (-Icges(6,2) * t232 + t295) * t229;
t212 = (-Icges(6,1) * t232 - t295) * t229;
t267 = -t198 + t212;
t402 = ((t209 + t208) * t230 + ((-t266 - t267) * t233 + (t264 + t265) * t232) * t229) * t230;
t401 = t390 * t204 + t391 * t205 + t389 * t293;
t400 = -t229 * (-(t201 / 0.2e1 + t211 / 0.2e1 + t200 / 0.2e1 + t210 / 0.2e1) * t232 + (t213 / 0.2e1 - t199 / 0.2e1 + t212 / 0.2e1 - t198 / 0.2e1) * t233) + (t209 / 0.2e1 + t208 / 0.2e1) * t230;
t183 = Icges(6,4) * t205;
t135 = -Icges(6,2) * t204 + Icges(6,6) * t293 + t183;
t186 = Icges(5,4) * t205;
t138 = -Icges(5,2) * t204 + Icges(5,6) * t293 + t186;
t182 = Icges(6,4) * t204;
t142 = -Icges(6,1) * t205 - Icges(6,5) * t293 + t182;
t185 = Icges(5,4) * t204;
t145 = -Icges(5,1) * t205 - Icges(5,5) * t293 + t185;
t399 = (t142 + t145) * t207 + t369 * (t135 + t138);
t223 = t226 * qJ(3);
t371 = t223 + (-pkin(2) + (-rSges(6,3) + t231) * t229) * t225 + t395;
t313 = pkin(4) * t232;
t94 = (qJ(3) + t313) * t225 + (rSges(6,3) * t229 + pkin(2) + t294) * t226 + t353;
t393 = m(6) * (t94 * t290 - t293 * t371);
t388 = t401 * t290;
t386 = t388 - t399;
t385 = qJD(2) * t393;
t384 = t121 * t226 - t271 * t225;
t297 = Icges(6,4) * t207;
t137 = -Icges(6,2) * t369 + Icges(6,6) * t290 + t297;
t300 = Icges(5,4) * t207;
t140 = -Icges(5,2) * t369 + Icges(5,6) * t290 + t300;
t184 = Icges(6,4) * t369;
t143 = Icges(6,1) * t207 + Icges(6,5) * t290 - t184;
t187 = Icges(5,4) * t369;
t146 = Icges(5,1) * t207 + Icges(5,5) * t290 - t187;
t381 = (t143 + t146) * t207 + (-t137 - t140) * t369 + (t391 * t207 + t389 * t290 + t390 * t369) * t290;
t346 = (rSges(5,3) + pkin(6)) * t229 + pkin(2) + t314;
t370 = -t225 * t346 + t223 + t354;
t367 = t386 - t388;
t364 = t229 / 0.2e1;
t363 = -t230 / 0.2e1;
t360 = m(6) * t229;
t275 = -Icges(5,2) * t207 + t146 - t187;
t277 = -Icges(6,2) * t207 + t143 - t184;
t352 = t275 + t277;
t276 = -Icges(5,2) * t205 - t145 - t185;
t278 = -Icges(6,2) * t205 - t142 - t182;
t351 = t276 + t278;
t279 = -Icges(5,1) * t369 - t140 - t300;
t281 = -Icges(6,1) * t369 - t137 - t297;
t350 = t279 + t281;
t280 = -Icges(5,1) * t204 - t138 - t186;
t282 = -Icges(6,1) * t204 - t135 - t183;
t349 = t280 + t282;
t158 = -Icges(6,5) * t204 - Icges(6,6) * t205;
t159 = -Icges(6,5) * t369 - Icges(6,6) * t207;
t160 = -Icges(5,5) * t204 - Icges(5,6) * t205;
t161 = -Icges(5,5) * t369 - Icges(5,6) * t207;
t344 = (t160 + t158) * t293 + (t159 + t161) * t290;
t340 = 0.4e1 * qJD(2);
t339 = 2 * qJD(4);
t106 = qJ(3) * t225 + t226 * t346 + t247;
t331 = m(5) * (t106 * t173 - t171 * t370);
t329 = m(5) * (t106 * t225 + t370 * t226);
t324 = m(6) * (-t69 * t290 - t293 * t394);
t322 = m(6) * (t121 * t371 - t271 * t94);
t321 = m(6) * (t94 * t225 + t371 * t226);
t320 = t384 * t360;
t316 = m(4) * ((rSges(4,2) * t293 + rSges(4,3) * t226 + t223) * t226 + (-rSges(4,2) * t290 + (rSges(4,3) + qJ(3)) * t225) * t225);
t157 = (-t225 ^ 2 - t226 ^ 2) * t360;
t263 = t157 * qJD(2);
t244 = (t367 + t399) * t364 * t226 + (t363 + t230 / 0.2e1) * ((-t389 * t230 + (t390 * t232 + t391 * t233) * t229) * t293 + (t201 + t200) * t205 + (-t199 - t198) * t204);
t243 = -(t386 * t225 + t381 * t226) * t229 / 0.2e1 + ((t401 * t293 + t381) * t226 + t367 * t225) * t364;
t242 = (rSges(6,1) * t232 + rSges(6,2) * t233 + t313) * t229 ^ 2;
t215 = (-rSges(5,1) * t232 - rSges(5,2) * t233) * t229;
t120 = -t173 * t230 - t215 * t290;
t119 = t171 * t230 + t215 * t293;
t90 = (t171 * t226 - t173 * t225) * t229;
t87 = t242 * t226 + t271 * t230;
t86 = -t121 * t230 - t242 * t225;
t71 = t320 / 0.2e1;
t70 = t384 * t229;
t62 = t266 * t207 + t209 * t290 - t264 * t369;
t61 = t267 * t207 + t208 * t290 - t265 * t369;
t60 = -t264 * t204 + t266 * t205 + t209 * t293;
t59 = -t265 * t204 + t267 * t205 + t208 * t293;
t38 = -t230 * t161 + (-t275 * t232 + t279 * t233) * t229;
t37 = -t230 * t160 + (-t276 * t232 + t280 * t233) * t229;
t36 = -t230 * t159 + (-t277 * t232 + t281 * t233) * t229;
t35 = -t230 * t158 + (-t278 * t232 + t282 * t233) * t229;
t26 = t324 / 0.2e1;
t23 = t316 + t321 + t329;
t15 = t322 + t331 - t400;
t7 = t26 + t71;
t6 = t26 - t320 / 0.2e1;
t5 = t71 - t324 / 0.2e1;
t3 = t307 + t309;
t1 = (t243 * t225 + t244 * t226) * t229;
t2 = [0, 0, 0, (-t70 * t335 + t90 * t337) * t339, 0; 0, t23 * qJD(3) + t15 * qJD(4) + qJD(5) * t393, qJD(2) * t23 + qJD(4) * t3, t15 * qJD(2) + t3 * qJD(3) + t7 * qJD(5) + ((-t103 * t173 + t106 * t120 + t119 * t370 - t171 * t387) * t337 + (t121 * t394 + t271 * t69 + t371 * t86 + t87 * t94) * t335) * t339 + (t402 + ((t38 / 0.2e1 + t36 / 0.2e1 + t62 / 0.2e1 + t61 / 0.2e1 - t244) * t226 + (t37 / 0.2e1 + t35 / 0.2e1 + t60 / 0.2e1 + t59 / 0.2e1 - t243) * t225) * t229) * qJD(4), t7 * qJD(4) + t385; 0, -t4 * qJD(4) + t157 * qJD(5) + (-t316 / 0.4e1 - t329 / 0.4e1 - t321 / 0.4e1) * t340, 0, -t403 + ((t119 * t225 - t120 * t226) * t337 + (t225 * t86 - t226 * t87) * t335) * t339, t263; 0, t4 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + (-t331 / 0.4e1 - t322 / 0.4e1) * t340 + t400 * qJD(2), t403, t1 * qJD(2) + ((t402 + ((t36 + t38) * t226 + (t35 + t37) * t225) * t229) * t363 + ((-t204 * t352 + t205 * t350) * t290 + (-t60 - t59) * t230 + (-t351 * t204 + t349 * t205 + t344) * t293) * t293 / 0.2e1 + ((t207 * t349 - t351 * t369) * t293 + (-t62 - t61) * t230 + (t350 * t207 - t352 * t369 + t344) * t290) * t290 / 0.2e1 + m(5) * (t387 * t119 - t103 * t120 + (t150 * t226 - t154 * t225) * t229 * t90) + m(6) * (-(-t273 * t225 + t274 * t226) * t229 * t70 + t394 * t86 - t69 * t87)) * qJD(4), t6 * qJD(2); 0, -t157 * qJD(3) + t5 * qJD(4) - t385, -t263, t5 * qJD(2) + m(6) * (t230 * t70 + (t225 * t87 + t226 * t86) * t229) * qJD(4), 0;];
Cq = t2;
