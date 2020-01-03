% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:21
% EndTime: 2020-01-03 12:00:26
% DurationCPUTime: 2.38s
% Computational Cost: add. (25408->225), mult. (12475->300), div. (0->0), fcn. (10650->10), ass. (0->151)
t225 = qJ(1) + qJ(2);
t220 = pkin(9) + t225;
t219 = qJ(4) + t220;
t213 = sin(t219);
t214 = cos(t219);
t176 = rSges(5,1) * t213 + rSges(5,2) * t214;
t216 = sin(t220);
t221 = sin(t225);
t288 = pkin(2) * t221;
t242 = pkin(3) * t216 + t288;
t158 = t176 + t242;
t177 = t214 * rSges(5,1) - rSges(5,2) * t213;
t217 = cos(t220);
t222 = cos(t225);
t218 = pkin(2) * t222;
t247 = pkin(3) * t217 + t218;
t159 = t177 + t247;
t255 = -t158 * t177 + t176 * t159;
t226 = sin(qJ(5));
t264 = t213 * t226;
t250 = -rSges(6,2) * t264 - t214 * rSges(6,3);
t228 = cos(qJ(5));
t283 = rSges(6,1) * t228;
t137 = -t214 * pkin(8) + (pkin(4) + t283) * t213 + t250;
t126 = t137 + t242;
t260 = t214 * t228;
t261 = t214 * t226;
t245 = rSges(6,1) * t260 - rSges(6,2) * t261 + t213 * rSges(6,3);
t138 = t214 * pkin(4) + t213 * pkin(8) + t245;
t127 = t138 + t247;
t284 = -t126 * t138 + t137 * t127;
t320 = m(6) / 0.2e1;
t321 = m(5) / 0.2e1;
t289 = sin(qJ(1)) * pkin(1);
t118 = t126 + t289;
t224 = cos(qJ(1)) * pkin(1);
t119 = t224 + t127;
t51 = t118 * t138 - t137 * t119;
t149 = t158 + t289;
t150 = t159 + t224;
t89 = t149 * t177 - t176 * t150;
t285 = (t284 + t51) * t320 + (t255 + t89) * t321;
t286 = (-t284 + t51) * t320 + (-t255 + t89) * t321;
t5 = t285 - t286;
t341 = t5 * qJD(1);
t223 = Icges(6,4) * t228;
t198 = -Icges(6,2) * t226 + t223;
t334 = Icges(6,1) * t226 + t223;
t340 = t198 + t334;
t319 = m(4) * (-t224 * (rSges(4,1) * t216 + rSges(4,2) * t217 + t288) + t289 * (t217 * rSges(4,1) - rSges(4,2) * t216 + t218));
t292 = m(3) * (-t224 * (rSges(3,1) * t221 + rSges(3,2) * t222) + t289 * (t222 * rSges(3,1) - rSges(3,2) * t221));
t203 = -rSges(6,2) * t226 + t283;
t338 = m(6) * t203;
t281 = Icges(6,4) * t226;
t197 = Icges(6,2) * t228 + t281;
t200 = Icges(6,1) * t228 - t281;
t243 = t340 * t228 / 0.2e1 + (-t197 / 0.2e1 + t200 / 0.2e1) * t226;
t201 = rSges(6,1) * t226 + rSges(6,2) * t228;
t172 = t201 * t213;
t173 = t201 * t214;
t63 = -t118 * t172 - t119 * t173;
t275 = t138 * t173;
t67 = -t137 * t172 - t275;
t308 = m(6) * (t67 + t63);
t336 = t243 + t308 / 0.2e1;
t64 = -t126 * t172 - t127 * t173;
t306 = m(6) * (t67 + t64);
t335 = t243 + t306 / 0.2e1;
t70 = t149 * t159 - t158 * t150;
t332 = (t197 - t200) * t228 + t340 * t226;
t192 = Icges(6,4) * t261;
t156 = Icges(6,1) * t260 + Icges(6,5) * t213 - t192;
t251 = -Icges(6,2) * t260 + t156 - t192;
t154 = Icges(6,4) * t260 - Icges(6,2) * t261 + Icges(6,6) * t213;
t253 = t214 * t334 + t154;
t331 = -t226 * t251 - t228 * t253;
t155 = -Icges(6,5) * t214 + t200 * t213;
t252 = -t197 * t213 + t155;
t153 = -Icges(6,6) * t214 + t198 * t213;
t254 = t213 * t334 + t153;
t330 = -t226 * t252 - t228 * t254;
t329 = t213 ^ 2;
t328 = t214 ^ 2;
t327 = 0.4e1 * qJD(1);
t324 = 2 * qJD(4);
t315 = m(5) * t70;
t313 = m(5) * t89;
t312 = m(5) * t255;
t309 = m(6) * (t64 + t63);
t307 = m(6) * ((-t119 + t127) * t214 + (-t118 + t126) * t213) * t201;
t305 = m(6) * (t275 + (-t119 * t214 + (-t118 + t137) * t213) * t201);
t304 = m(6) * (t275 + (-t127 * t214 + (-t126 + t137) * t213) * t201);
t48 = t118 * t127 - t126 * t119;
t303 = m(6) * t48;
t301 = m(6) * t51;
t300 = m(6) * t284;
t298 = m(6) * t63;
t297 = m(6) * t64;
t296 = m(6) * t67;
t295 = -t213 / 0.2e1;
t294 = -t214 / 0.2e1;
t293 = t214 / 0.2e1;
t263 = t213 * t228;
t143 = t155 * t263;
t144 = t156 * t263;
t145 = t153 * t261;
t196 = Icges(6,5) * t228 - Icges(6,6) * t226;
t265 = t196 * t213;
t151 = -Icges(6,3) * t214 + t265;
t273 = t154 * t226;
t234 = t156 * t228 - t273;
t272 = t155 * t228;
t274 = t153 * t226;
t79 = -t151 * t213 - t155 * t260 + t145;
t152 = Icges(6,5) * t260 - Icges(6,6) * t261 + Icges(6,3) * t213;
t80 = t152 * t213 + t234 * t214;
t17 = (t79 + t144 - t145 + (t151 - t273) * t213) * t213 + (-t143 - t80 + (t151 + t234) * t214 + (t272 + t274) * t213) * t214;
t77 = -t151 * t214 - t153 * t264 + t143;
t78 = t152 * t214 + t154 * t264 - t144;
t18 = (t145 - t78 + (t152 - t272) * t214) * t214 + (-t143 + t77 + (t152 + t274) * t213) * t213;
t46 = -t213 * t78 - t214 * t77;
t47 = -t213 * t80 - t214 * t79;
t2 = (-t18 / 0.2e1 - t47 / 0.2e1) * t214 + (t46 / 0.2e1 - t17 / 0.2e1) * t213;
t287 = t2 * qJD(5);
t240 = t309 / 0.2e1 + t243;
t236 = Icges(6,5) * t226 + Icges(6,6) * t228;
t232 = t213 * t17 / 0.2e1 + (-t196 * t214 - t332 * t213 - t254 * t226 + t252 * t228) * t294 + (t18 + t47) * t293 + (t332 * t214 + t253 * t226 - t251 * t228 - t265 + t46) * t295;
t231 = -t243 + (t293 + t294) * (t228 * t154 + t226 * t156);
t167 = t236 * t214;
t166 = t213 * t236;
t125 = -t172 * t213 - t173 * t214;
t59 = t243 + t296;
t52 = t243 + t297;
t50 = t243 + t298;
t40 = t304 / 0.2e1;
t38 = t305 / 0.2e1;
t33 = t307 / 0.2e1;
t31 = -t300 - t312;
t30 = t301 + t313;
t19 = t292 + t303 + t315 + t319;
t14 = -t304 / 0.2e1 + t335;
t13 = t40 + t335;
t12 = -t305 / 0.2e1 + t336;
t11 = t38 + t336;
t10 = -t307 / 0.2e1 + t240;
t9 = t33 + t240;
t8 = t40 - t306 / 0.2e1 + t231;
t7 = t38 - t308 / 0.2e1 + t231;
t6 = t33 - t309 / 0.2e1 + t231;
t4 = t285 + t286;
t1 = [qJD(2) * t19 + qJD(4) * t30 + qJD(5) * t50, t19 * qJD(1) + t4 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t319 / 0.2e1 + t292 / 0.2e1 + t48 * t320 + t70 * t321) * qJD(2), 0, t30 * qJD(1) + t4 * qJD(2) + t11 * qJD(5) + (t51 * t320 + t89 * t321) * t324, t50 * qJD(1) + t9 * qJD(2) + t11 * qJD(4) + ((t118 * t214 - t119 * t213) * t338 + t232) * qJD(5); -t5 * qJD(4) + t10 * qJD(5) + (-t292 / 0.4e1 - t319 / 0.4e1 - t315 / 0.4e1 - t303 / 0.4e1) * t327, qJD(4) * t31 + qJD(5) * t52, 0, -t341 + t31 * qJD(2) + t13 * qJD(5) + (-t255 * t321 - t284 * t320) * t324, t10 * qJD(1) + t52 * qJD(2) + t13 * qJD(4) + ((t126 * t214 - t127 * t213) * t338 + t232) * qJD(5); 0, 0, 0, 0, m(6) * t125 * qJD(5); t5 * qJD(2) + t12 * qJD(5) + (-t313 / 0.4e1 - t301 / 0.4e1) * t327, t341 + t14 * qJD(5) + 0.4e1 * (t300 / 0.4e1 + t312 / 0.4e1) * qJD(2), 0, qJD(5) * t59, t12 * qJD(1) + t14 * qJD(2) + t59 * qJD(4) + ((t137 * t214 - t138 * t213) * t338 + t232) * qJD(5); (t231 - t298) * qJD(1) + t6 * qJD(2) + t7 * qJD(4) + t287, t6 * qJD(1) + (t231 - t297) * qJD(2) + t8 * qJD(4) + t287, 0, t7 * qJD(1) + t8 * qJD(2) + (t231 - t296) * qJD(4) + t287, (m(6) * ((t214 * t245 + t213 * (rSges(6,1) * t263 + t250)) * t125 + (t328 + t329) * t203 * t201) + (-t328 * t166 + (t331 * t213 + (t167 - t330) * t214) * t213) * t294 + (t329 * t167 + (t330 * t214 + (-t166 - t331) * t213) * t214) * t295) * qJD(5) + (qJD(1) + qJD(2) + qJD(4)) * t2;];
Cq = t1;
