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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:33:59
% EndTime: 2022-01-20 10:34:05
% DurationCPUTime: 2.55s
% Computational Cost: add. (25408->224), mult. (12475->302), div. (0->0), fcn. (10650->10), ass. (0->150)
t222 = qJ(1) + qJ(2);
t218 = pkin(9) + t222;
t217 = qJ(4) + t218;
t212 = sin(t217);
t213 = cos(t217);
t177 = -rSges(5,1) * t212 - rSges(5,2) * t213;
t215 = sin(t218);
t219 = sin(t222);
t289 = pkin(2) * t219;
t243 = -pkin(3) * t215 - t289;
t159 = t177 + t243;
t178 = t213 * rSges(5,1) - t212 * rSges(5,2);
t216 = cos(t218);
t220 = cos(t222);
t288 = pkin(2) * t220;
t242 = pkin(3) * t216 + t288;
t160 = t178 + t242;
t257 = t178 * t159 - t160 * t177;
t225 = cos(qJ(5));
t283 = rSges(6,1) * t225;
t248 = pkin(4) + t283;
t223 = sin(qJ(5));
t266 = t212 * t223;
t250 = rSges(6,2) * t266 + t213 * rSges(6,3);
t137 = t213 * pkin(8) - t248 * t212 + t250;
t126 = t137 + t243;
t263 = t213 * t223;
t195 = rSges(6,2) * t263;
t138 = -t195 + t248 * t213 + (rSges(6,3) + pkin(8)) * t212;
t127 = t138 + t242;
t284 = t138 * t126 - t127 * t137;
t322 = m(6) / 0.2e1;
t323 = m(5) / 0.2e1;
t291 = sin(qJ(1)) * pkin(1);
t119 = t126 - t291;
t290 = cos(qJ(1)) * pkin(1);
t120 = t127 + t290;
t51 = -t138 * t119 + t120 * t137;
t150 = t159 - t291;
t151 = t160 + t290;
t89 = -t178 * t150 + t151 * t177;
t285 = (t284 + t51) * t322 + (t257 + t89) * t323;
t286 = (-t284 + t51) * t322 + (-t257 + t89) * t323;
t5 = t285 - t286;
t338 = t5 * qJD(1);
t221 = Icges(6,4) * t225;
t198 = -Icges(6,2) * t223 + t221;
t199 = Icges(6,1) * t223 + t221;
t337 = t198 + t199;
t321 = m(4) * (t290 * (-rSges(4,1) * t215 - rSges(4,2) * t216 - t289) + (t216 * rSges(4,1) - t215 * rSges(4,2) + t288) * t291);
t294 = m(3) * (t290 * (-rSges(3,1) * t219 - rSges(3,2) * t220) + (t220 * rSges(3,1) - t219 * rSges(3,2)) * t291);
t281 = Icges(6,4) * t223;
t197 = Icges(6,2) * t225 + t281;
t200 = Icges(6,1) * t225 - t281;
t244 = t337 * t225 / 0.2e1 + (-t197 / 0.2e1 + t200 / 0.2e1) * t223;
t201 = rSges(6,1) * t223 + rSges(6,2) * t225;
t173 = t201 * t212;
t108 = t119 * t173;
t174 = t201 * t213;
t63 = -t120 * t174 + t108;
t67 = t137 * t173 - t138 * t174;
t310 = m(6) * (t67 + t63);
t334 = t244 + t310 / 0.2e1;
t64 = t126 * t173 - t127 * t174;
t308 = m(6) * (t67 + t64);
t333 = t244 + t308 / 0.2e1;
t70 = -t160 * t150 + t151 * t159;
t331 = t212 ^ 2;
t330 = t213 ^ 2;
t329 = 0.4e1 * qJD(1);
t326 = 2 * qJD(4);
t317 = m(5) * t70;
t315 = m(5) * t89;
t314 = m(5) * t257;
t311 = m(6) * (t64 + t63);
t309 = m(6) * (t108 + (-t126 * t212 + (-t120 + t127) * t213) * t201);
t307 = m(6) * (t108 + (-t137 * t212 + (-t120 + t138) * t213) * t201);
t306 = m(6) * ((-t127 + t138) * t213 + (t126 - t137) * t212) * t201;
t48 = -t127 * t119 + t120 * t126;
t305 = m(6) * t48;
t303 = m(6) * t51;
t302 = m(6) * t284;
t300 = m(6) * t63;
t299 = m(6) * t64;
t298 = m(6) * t67;
t297 = -t212 / 0.2e1;
t296 = t212 / 0.2e1;
t295 = -t213 / 0.2e1;
t265 = t212 * t225;
t152 = Icges(6,5) * t265 - Icges(6,6) * t266 - Icges(6,3) * t213;
t155 = Icges(6,6) * t212 + t198 * t213;
t245 = t155 * t223 - t152;
t157 = Icges(6,5) * t212 + t200 * t213;
t143 = t157 * t265;
t196 = Icges(6,5) * t225 - Icges(6,6) * t223;
t268 = t196 * t213;
t153 = Icges(6,3) * t212 + t268;
t246 = t213 * t153 - t143;
t262 = t213 * t225;
t255 = t212 * t153 + t157 * t262;
t154 = Icges(6,4) * t265 - Icges(6,2) * t266 - Icges(6,6) * t213;
t193 = Icges(6,4) * t266;
t156 = Icges(6,1) * t265 - Icges(6,5) * t213 - t193;
t256 = -t212 * t152 - t156 * t262;
t79 = -t154 * t263 - t256;
t80 = -t155 * t263 + t255;
t17 = (t245 * t213 - t255 + t80) * t213 + (t245 * t212 + t246 + t79) * t212;
t275 = t154 * t223;
t78 = -t155 * t266 - t246;
t18 = (t78 - t143 + (t153 + t275) * t213 + t256) * t213 + t255 * t212;
t46 = t212 * t78 - t213 * (-(-t156 * t225 + t275) * t212 - t213 * t152);
t47 = t212 * t80 - t213 * t79;
t2 = (t47 / 0.2e1 - t18 / 0.2e1) * t213 + (t17 / 0.2e1 + t46 / 0.2e1) * t212;
t287 = t2 * qJD(5);
t254 = t199 * t212 + t154;
t253 = -t199 * t213 - t155;
t252 = -Icges(6,2) * t265 + t156 - t193;
t251 = -t197 * t213 + t157;
t241 = t311 / 0.2e1 + t244;
t237 = Icges(6,5) * t223 + Icges(6,6) * t225;
t233 = (-t173 * t213 + t174 * t212) * t201;
t228 = (-t197 + t200) * t225 - t337 * t223;
t232 = t213 * t18 / 0.2e1 + (t17 + t46) * t297 + (t196 * t212 + t228 * t213 + t253 * t223 + t251 * t225) * t296 + (t228 * t212 - t254 * t223 + t252 * t225 - t268 + t47) * t295;
t231 = -t244 + (t296 + t297) * (t225 * t154 + t223 * t156);
t230 = t252 * t223 + t254 * t225;
t229 = -t251 * t223 + t253 * t225;
t203 = -rSges(6,2) * t223 + t283;
t168 = t213 * t237;
t167 = t237 * t212;
t125 = -t173 * t212 - t174 * t213;
t59 = t244 + t298;
t52 = t244 + t299;
t50 = t244 + t300;
t40 = t306 / 0.2e1;
t38 = t307 / 0.2e1;
t33 = t309 / 0.2e1;
t31 = -t302 - t314;
t30 = t303 + t315;
t19 = t294 + t305 + t317 + t321;
t14 = -t306 / 0.2e1 + t333;
t13 = t40 + t333;
t12 = -t307 / 0.2e1 + t334;
t11 = t38 + t334;
t10 = -t309 / 0.2e1 + t241;
t9 = t33 + t241;
t8 = t40 - t308 / 0.2e1 + t231;
t7 = t38 - t310 / 0.2e1 + t231;
t6 = t33 - t311 / 0.2e1 + t231;
t3 = t285 + t286;
t1 = [qJD(2) * t19 + qJD(4) * t30 + qJD(5) * t50, t19 * qJD(1) + t3 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t321 / 0.2e1 + t294 / 0.2e1 + t48 * t322 + t70 * t323) * qJD(2), 0, t30 * qJD(1) + t3 * qJD(2) + t11 * qJD(5) + (t51 * t322 + t89 * t323) * t326, t50 * qJD(1) + t9 * qJD(2) + t11 * qJD(4) + (((-t119 * t213 - t120 * t212) * t203 + t233) * m(6) + t232) * qJD(5); -t5 * qJD(4) + t10 * qJD(5) + (-t305 / 0.4e1 - t317 / 0.4e1 - t294 / 0.4e1 - t321 / 0.4e1) * t329, qJD(4) * t31 + qJD(5) * t52, 0, -t338 + t31 * qJD(2) + t13 * qJD(5) + (-t257 * t323 - t284 * t322) * t326, t10 * qJD(1) + t52 * qJD(2) + t13 * qJD(4) + (((-t126 * t213 - t127 * t212) * t203 + t233) * m(6) + t232) * qJD(5); 0, 0, 0, 0, m(6) * t125 * qJD(5); t5 * qJD(2) + t12 * qJD(5) + (-t303 / 0.4e1 - t315 / 0.4e1) * t329, t338 + t14 * qJD(5) + 0.4e1 * (t302 / 0.4e1 + t314 / 0.4e1) * qJD(2), 0, qJD(5) * t59, t12 * qJD(1) + t14 * qJD(2) + t59 * qJD(4) + (((-t137 * t213 - t138 * t212) * t203 + t233) * m(6) + t232) * qJD(5); (t231 - t300) * qJD(1) + t6 * qJD(2) + t7 * qJD(4) + t287, t6 * qJD(1) + (t231 - t299) * qJD(2) + t8 * qJD(4) + t287, 0, t7 * qJD(1) + t8 * qJD(2) + (t231 - t298) * qJD(4) + t287, (m(6) * ((t212 * (rSges(6,1) * t265 - t250) + t213 * (rSges(6,1) * t262 + t212 * rSges(6,3) - t195)) * t125 + (t330 + t331) * t203 * t201) + (-t331 * t168 + (t230 * t213 + (t167 + t229) * t212) * t213) * t296 + (-t330 * t167 + (t229 * t212 + (t168 + t230) * t213) * t212) * t295) * qJD(5) + (qJD(1) + qJD(2) + qJD(4)) * t2;];
Cq = t1;
