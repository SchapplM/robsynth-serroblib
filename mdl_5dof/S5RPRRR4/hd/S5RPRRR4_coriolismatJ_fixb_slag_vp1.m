% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:33
% DurationCPUTime: 2.32s
% Computational Cost: add. (24670->217), mult. (12027->292), div. (0->0), fcn. (10264->10), ass. (0->147)
t213 = qJ(1) + pkin(9);
t211 = qJ(3) + t213;
t208 = qJ(4) + t211;
t204 = sin(t208);
t205 = cos(t208);
t169 = -rSges(5,1) * t204 - rSges(5,2) * t205;
t206 = sin(t211);
t281 = pkin(3) * t206;
t159 = t169 - t281;
t170 = t205 * rSges(5,1) - t204 * rSges(5,2);
t207 = cos(t211);
t280 = pkin(3) * t207;
t160 = t170 + t280;
t248 = t170 * t159 - t160 * t169;
t216 = cos(qJ(5));
t275 = rSges(6,1) * t216;
t239 = pkin(4) + t275;
t214 = sin(qJ(5));
t259 = t204 * t214;
t241 = rSges(6,2) * t259 + t205 * rSges(6,3);
t131 = t205 * pkin(8) - t239 * t204 + t241;
t125 = t131 - t281;
t274 = rSges(6,2) * t214;
t187 = t205 * t274;
t132 = -t187 + t239 * t205 + (rSges(6,3) + pkin(8)) * t204;
t126 = t132 + t280;
t276 = t132 * t125 - t126 * t131;
t310 = m(6) / 0.2e1;
t311 = m(5) / 0.2e1;
t234 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t213);
t116 = t125 + t234;
t233 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t213);
t117 = t126 + t233;
t50 = -t132 * t116 + t117 * t131;
t146 = t159 + t234;
t147 = t160 + t233;
t84 = -t170 * t146 + t147 * t169;
t277 = (t276 + t50) * t310 + (t248 + t84) * t311;
t278 = (-t276 + t50) * t310 + (-t248 + t84) * t311;
t5 = t277 - t278;
t325 = t5 * qJD(1);
t212 = Icges(6,4) * t216;
t190 = -Icges(6,2) * t214 + t212;
t191 = Icges(6,1) * t214 + t212;
t324 = t190 + t191;
t284 = m(4) * (t233 * (-rSges(4,1) * t206 - rSges(4,2) * t207) - (t207 * rSges(4,1) - t206 * rSges(4,2)) * t234);
t272 = Icges(6,4) * t214;
t189 = Icges(6,2) * t216 + t272;
t192 = Icges(6,1) * t216 - t272;
t235 = t324 * t216 / 0.2e1 + (-t189 / 0.2e1 + t192 / 0.2e1) * t214;
t193 = rSges(6,1) * t214 + rSges(6,2) * t216;
t167 = t193 * t204;
t105 = t116 * t167;
t168 = t193 * t205;
t62 = -t117 * t168 + t105;
t67 = t131 * t167 - t132 * t168;
t300 = m(6) * (t67 + t62);
t321 = t235 + t300 / 0.2e1;
t64 = t125 * t167 - t126 * t168;
t298 = m(6) * (t67 + t64);
t320 = t235 + t298 / 0.2e1;
t77 = -t160 * t146 + t147 * t159;
t318 = t204 ^ 2;
t317 = t205 ^ 2;
t316 = 0.4e1 * qJD(1);
t313 = 2 * qJD(4);
t307 = m(5) * t77;
t306 = m(5) * t84;
t304 = m(5) * t248;
t301 = m(6) * (t64 + t62);
t299 = m(6) * (t105 + (-t125 * t204 + (-t117 + t126) * t205) * t193);
t297 = m(6) * (t105 + (-t131 * t204 + (-t117 + t132) * t205) * t193);
t296 = m(6) * ((-t126 + t132) * t205 + (t125 - t131) * t204) * t193;
t48 = -t126 * t116 + t117 * t125;
t295 = m(6) * t48;
t294 = m(6) * t50;
t292 = m(6) * t276;
t290 = m(6) * t62;
t289 = m(6) * t64;
t288 = m(6) * t67;
t287 = -t204 / 0.2e1;
t286 = t204 / 0.2e1;
t285 = -t205 / 0.2e1;
t258 = t204 * t216;
t148 = Icges(6,5) * t258 - Icges(6,6) * t259 - Icges(6,3) * t205;
t151 = Icges(6,6) * t204 + t190 * t205;
t255 = t214 * t151;
t236 = -t148 + t255;
t153 = Icges(6,5) * t204 + t192 * t205;
t251 = t216 * t153;
t137 = t204 * t251;
t188 = Icges(6,5) * t216 - Icges(6,6) * t214;
t261 = t188 * t205;
t149 = Icges(6,3) * t204 + t261;
t237 = t205 * t149 - t137;
t246 = t204 * t149 + t205 * t251;
t185 = Icges(6,4) * t259;
t152 = Icges(6,1) * t258 - Icges(6,5) * t205 - t185;
t252 = t216 * t152;
t247 = -t204 * t148 - t205 * t252;
t150 = Icges(6,4) * t258 - Icges(6,2) * t259 - Icges(6,6) * t205;
t256 = t214 * t150;
t75 = -t205 * t256 - t247;
t76 = -t205 * t255 + t246;
t15 = (t236 * t205 - t246 + t76) * t205 + (t236 * t204 + t237 + t75) * t204;
t74 = -t204 * t255 - t237;
t16 = (t74 - t137 + (t149 + t256) * t205 + t247) * t205 + t246 * t204;
t46 = t204 * t74 - t205 * (-(-t252 + t256) * t204 - t205 * t148);
t47 = t204 * t76 - t205 * t75;
t2 = (t47 / 0.2e1 - t16 / 0.2e1) * t205 + (t15 / 0.2e1 + t46 / 0.2e1) * t204;
t279 = t2 * qJD(5);
t245 = t191 * t204 + t150;
t244 = -t191 * t205 - t151;
t243 = -Icges(6,2) * t258 + t152 - t185;
t242 = -t189 * t205 + t153;
t232 = t301 / 0.2e1 + t235;
t228 = Icges(6,5) * t214 + Icges(6,6) * t216;
t224 = (-t167 * t205 + t168 * t204) * t193;
t219 = (-t189 + t192) * t216 - t324 * t214;
t223 = t205 * t16 / 0.2e1 + (t15 + t46) * t287 + (t188 * t204 + t219 * t205 + t244 * t214 + t242 * t216) * t286 + (t219 * t204 - t245 * t214 + t243 * t216 - t261 + t47) * t285;
t222 = -t235 + (t286 + t287) * (t216 * t150 + t214 * t152);
t221 = t243 * t214 + t245 * t216;
t220 = -t242 * t214 + t244 * t216;
t195 = -t274 + t275;
t162 = t205 * t228;
t161 = t228 * t204;
t121 = -t167 * t204 - t168 * t205;
t54 = t235 + t288;
t52 = t235 + t289;
t49 = t235 + t290;
t40 = t296 / 0.2e1;
t38 = t297 / 0.2e1;
t36 = -t292 - t304;
t33 = t299 / 0.2e1;
t28 = t294 + t306;
t19 = t284 + t295 + t307;
t18 = -t296 / 0.2e1 + t320;
t17 = t40 + t320;
t12 = -t297 / 0.2e1 + t321;
t11 = t38 + t321;
t10 = -t299 / 0.2e1 + t232;
t9 = t33 + t232;
t8 = t40 - t298 / 0.2e1 + t222;
t7 = t38 - t300 / 0.2e1 + t222;
t6 = t33 - t301 / 0.2e1 + t222;
t3 = t277 + t278;
t1 = [qJD(3) * t19 + qJD(4) * t28 + qJD(5) * t49, 0, t19 * qJD(1) + t3 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t284 / 0.2e1 + t48 * t310 + t77 * t311) * qJD(3), t28 * qJD(1) + t3 * qJD(3) + t11 * qJD(5) + (t50 * t310 + t84 * t311) * t313, t49 * qJD(1) + t9 * qJD(3) + t11 * qJD(4) + (((-t116 * t205 - t117 * t204) * t195 + t224) * m(6) + t223) * qJD(5); 0, 0, 0, 0, m(6) * t121 * qJD(5); -t5 * qJD(4) + t10 * qJD(5) + (-t295 / 0.4e1 - t307 / 0.4e1 - t284 / 0.4e1) * t316, 0, qJD(4) * t36 + qJD(5) * t52, -t325 + t36 * qJD(3) + t17 * qJD(5) + (-t248 * t311 - t276 * t310) * t313, t10 * qJD(1) + t52 * qJD(3) + t17 * qJD(4) + (((-t125 * t205 - t126 * t204) * t195 + t224) * m(6) + t223) * qJD(5); t5 * qJD(3) + t12 * qJD(5) + (-t294 / 0.4e1 - t306 / 0.4e1) * t316, 0, t325 + t18 * qJD(5) + 0.4e1 * (t304 / 0.4e1 + t292 / 0.4e1) * qJD(3), t54 * qJD(5), t12 * qJD(1) + t18 * qJD(3) + t54 * qJD(4) + (((-t131 * t205 - t132 * t204) * t195 + t224) * m(6) + t223) * qJD(5); (t222 - t290) * qJD(1) + t6 * qJD(3) + t7 * qJD(4) + t279, 0, t6 * qJD(1) + (t222 - t289) * qJD(3) + t8 * qJD(4) + t279, t7 * qJD(1) + t8 * qJD(3) + (t222 - t288) * qJD(4) + t279, (m(6) * ((t204 * (rSges(6,1) * t258 - t241) + t205 * (t204 * rSges(6,3) + t205 * t275 - t187)) * t121 + (t317 + t318) * t195 * t193) + (-t318 * t162 + (t221 * t205 + (t161 + t220) * t204) * t205) * t286 + (-t317 * t161 + (t220 * t204 + (t162 + t221) * t205) * t204) * t285) * qJD(5) + (qJD(1) + qJD(3) + qJD(4)) * t2;];
Cq = t1;
