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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:52:03
% EndTime: 2020-01-03 11:52:08
% DurationCPUTime: 2.42s
% Computational Cost: add. (24670->218), mult. (12027->291), div. (0->0), fcn. (10264->10), ass. (0->146)
t216 = qJ(1) + pkin(9);
t213 = qJ(3) + t216;
t210 = qJ(4) + t213;
t205 = sin(t210);
t206 = cos(t210);
t168 = rSges(5,1) * t205 + rSges(5,2) * t206;
t208 = sin(t213);
t278 = pkin(3) * t208;
t158 = t168 + t278;
t169 = t206 * rSges(5,1) - rSges(5,2) * t205;
t209 = cos(t213);
t204 = pkin(3) * t209;
t159 = t169 + t204;
t246 = -t158 * t169 + t168 * t159;
t217 = sin(qJ(5));
t255 = t205 * t217;
t241 = -rSges(6,2) * t255 - t206 * rSges(6,3);
t219 = cos(qJ(5));
t273 = rSges(6,1) * t219;
t131 = -t206 * pkin(8) + (pkin(4) + t273) * t205 + t241;
t125 = t131 + t278;
t251 = t206 * t219;
t252 = t206 * t217;
t236 = rSges(6,1) * t251 - rSges(6,2) * t252 + t205 * rSges(6,3);
t132 = t206 * pkin(4) + t205 * pkin(8) + t236;
t126 = t204 + t132;
t274 = -t125 * t132 + t131 * t126;
t307 = m(6) / 0.2e1;
t308 = m(5) / 0.2e1;
t233 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t216);
t115 = t125 + t233;
t238 = pkin(2) * cos(t216) + cos(qJ(1)) * pkin(1);
t116 = t126 + t238;
t50 = t115 * t132 - t131 * t116;
t145 = t158 + t233;
t146 = t159 + t238;
t84 = t145 * t169 - t168 * t146;
t275 = (t274 + t50) * t307 + (t246 + t84) * t308;
t276 = (-t274 + t50) * t307 + (-t246 + t84) * t308;
t5 = t275 - t276;
t327 = t5 * qJD(1);
t214 = Icges(6,4) * t219;
t190 = -Icges(6,2) * t217 + t214;
t320 = Icges(6,1) * t217 + t214;
t326 = t190 + t320;
t281 = m(4) * (-t238 * (rSges(4,1) * t208 + rSges(4,2) * t209) + t233 * (t209 * rSges(4,1) - rSges(4,2) * t208));
t195 = -rSges(6,2) * t217 + t273;
t324 = m(6) * t195;
t271 = Icges(6,4) * t217;
t189 = Icges(6,2) * t219 + t271;
t192 = Icges(6,1) * t219 - t271;
t234 = t326 * t219 / 0.2e1 + (-t189 / 0.2e1 + t192 / 0.2e1) * t217;
t193 = rSges(6,1) * t217 + rSges(6,2) * t219;
t166 = t193 * t205;
t167 = t193 * t206;
t62 = -t115 * t166 - t116 * t167;
t265 = t132 * t167;
t67 = -t131 * t166 - t265;
t297 = m(6) * (t67 + t62);
t322 = t234 + t297 / 0.2e1;
t64 = -t125 * t166 - t126 * t167;
t295 = m(6) * (t67 + t64);
t321 = t234 + t295 / 0.2e1;
t77 = t145 * t159 - t158 * t146;
t318 = t326 * t217 + (t189 - t192) * t219;
t184 = Icges(6,4) * t252;
t152 = Icges(6,1) * t251 + Icges(6,5) * t205 - t184;
t242 = -Icges(6,2) * t251 + t152 - t184;
t150 = Icges(6,4) * t251 - Icges(6,2) * t252 + Icges(6,6) * t205;
t244 = t206 * t320 + t150;
t317 = -t242 * t217 - t244 * t219;
t151 = -Icges(6,5) * t206 + t192 * t205;
t243 = -t189 * t205 + t151;
t149 = -Icges(6,6) * t206 + t190 * t205;
t245 = t205 * t320 + t149;
t316 = -t243 * t217 - t245 * t219;
t315 = t205 ^ 2;
t314 = t206 ^ 2;
t313 = 0.4e1 * qJD(1);
t310 = 2 * qJD(4);
t304 = m(5) * t77;
t303 = m(5) * t84;
t301 = m(5) * t246;
t298 = m(6) * (t64 + t62);
t296 = m(6) * ((-t116 + t126) * t206 + (-t115 + t125) * t205) * t193;
t294 = m(6) * (t265 + (-t116 * t206 + (-t115 + t131) * t205) * t193);
t293 = m(6) * (t265 + (-t126 * t206 + (-t125 + t131) * t205) * t193);
t48 = t115 * t126 - t125 * t116;
t292 = m(6) * t48;
t291 = m(6) * t50;
t289 = m(6) * t274;
t287 = m(6) * t62;
t286 = m(6) * t64;
t285 = m(6) * t67;
t284 = -t205 / 0.2e1;
t283 = -t206 / 0.2e1;
t282 = t206 / 0.2e1;
t254 = t205 * t219;
t137 = t151 * t254;
t138 = t152 * t254;
t139 = t149 * t252;
t188 = Icges(6,5) * t219 - Icges(6,6) * t217;
t256 = t188 * t205;
t147 = -Icges(6,3) * t206 + t256;
t263 = t150 * t217;
t225 = t152 * t219 - t263;
t262 = t151 * t219;
t264 = t149 * t217;
t75 = -t205 * t147 - t151 * t251 + t139;
t148 = Icges(6,5) * t251 - Icges(6,6) * t252 + Icges(6,3) * t205;
t76 = t205 * t148 + t206 * t225;
t15 = (t75 + t138 - t139 + (t147 - t263) * t205) * t205 + (-t137 - t76 + (t147 + t225) * t206 + (t262 + t264) * t205) * t206;
t73 = -t206 * t147 - t149 * t255 + t137;
t74 = t206 * t148 + t150 * t255 - t138;
t16 = (t139 - t74 + (t148 - t262) * t206) * t206 + (-t137 + t73 + (t148 + t264) * t205) * t205;
t46 = -t205 * t74 - t206 * t73;
t47 = -t205 * t76 - t206 * t75;
t2 = (-t16 / 0.2e1 - t47 / 0.2e1) * t206 + (t46 / 0.2e1 - t15 / 0.2e1) * t205;
t277 = t2 * qJD(5);
t231 = t298 / 0.2e1 + t234;
t227 = Icges(6,5) * t217 + Icges(6,6) * t219;
t223 = t205 * t15 / 0.2e1 + (-t188 * t206 - t318 * t205 - t245 * t217 + t243 * t219) * t283 + (t16 + t47) * t282 + (t318 * t206 + t244 * t217 - t242 * t219 - t256 + t46) * t284;
t222 = -t234 + (t282 + t283) * (t219 * t150 + t217 * t152);
t161 = t227 * t206;
t160 = t205 * t227;
t121 = -t166 * t205 - t167 * t206;
t54 = t234 + t285;
t52 = t234 + t286;
t49 = t234 + t287;
t40 = t293 / 0.2e1;
t38 = t294 / 0.2e1;
t36 = -t289 - t301;
t33 = t296 / 0.2e1;
t28 = t291 + t303;
t19 = t281 + t292 + t304;
t18 = -t293 / 0.2e1 + t321;
t17 = t40 + t321;
t12 = -t294 / 0.2e1 + t322;
t11 = t38 + t322;
t10 = -t296 / 0.2e1 + t231;
t9 = t33 + t231;
t8 = t40 - t295 / 0.2e1 + t222;
t7 = t38 - t297 / 0.2e1 + t222;
t6 = t33 - t298 / 0.2e1 + t222;
t3 = t275 + t276;
t1 = [qJD(3) * t19 + qJD(4) * t28 + qJD(5) * t49, 0, t19 * qJD(1) + t3 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t281 / 0.2e1 + t48 * t307 + t77 * t308) * qJD(3), t28 * qJD(1) + t3 * qJD(3) + t11 * qJD(5) + (t50 * t307 + t84 * t308) * t310, t49 * qJD(1) + t9 * qJD(3) + t11 * qJD(4) + ((t115 * t206 - t116 * t205) * t324 + t223) * qJD(5); 0, 0, 0, 0, m(6) * t121 * qJD(5); -t5 * qJD(4) + t10 * qJD(5) + (-t281 / 0.4e1 - t304 / 0.4e1 - t292 / 0.4e1) * t313, 0, qJD(4) * t36 + qJD(5) * t52, -t327 + t36 * qJD(3) + t17 * qJD(5) + (-t246 * t308 - t274 * t307) * t310, t10 * qJD(1) + t52 * qJD(3) + t17 * qJD(4) + ((t125 * t206 - t126 * t205) * t324 + t223) * qJD(5); t5 * qJD(3) + t12 * qJD(5) + (-t303 / 0.4e1 - t291 / 0.4e1) * t313, 0, t327 + t18 * qJD(5) + 0.4e1 * (t301 / 0.4e1 + t289 / 0.4e1) * qJD(3), qJD(5) * t54, t12 * qJD(1) + t18 * qJD(3) + t54 * qJD(4) + ((t131 * t206 - t132 * t205) * t324 + t223) * qJD(5); (t222 - t287) * qJD(1) + t6 * qJD(3) + t7 * qJD(4) + t277, 0, t6 * qJD(1) + (t222 - t286) * qJD(3) + t8 * qJD(4) + t277, t7 * qJD(1) + t8 * qJD(3) + (t222 - t285) * qJD(4) + t277, (m(6) * ((t206 * t236 + t205 * (rSges(6,1) * t254 + t241)) * t121 + (t314 + t315) * t195 * t193) + (-t314 * t160 + (t317 * t205 + (t161 - t316) * t206) * t205) * t283 + (t315 * t161 + (t316 * t206 + (-t160 - t317) * t205) * t206) * t284) * qJD(5) + (qJD(1) + qJD(3) + qJD(4)) * t2;];
Cq = t1;
