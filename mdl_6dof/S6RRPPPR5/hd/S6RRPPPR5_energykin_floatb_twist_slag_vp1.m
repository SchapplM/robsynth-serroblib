% Calculate kinetic energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:16
% EndTime: 2019-03-09 08:21:19
% DurationCPUTime: 2.88s
% Computational Cost: add. (1339->320), mult. (2764->437), div. (0->0), fcn. (2900->8), ass. (0->150)
t321 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t320 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t319 = -Icges(4,4) + Icges(6,5) - Icges(5,6);
t318 = Icges(4,6) + Icges(6,4) - Icges(5,5);
t317 = Icges(6,6) - Icges(4,5) + Icges(5,4);
t316 = -Icges(4,3) - Icges(6,2) - Icges(5,1);
t256 = sin(pkin(9));
t260 = sin(qJ(1));
t262 = cos(qJ(2));
t295 = t260 * t262;
t257 = cos(pkin(9));
t263 = cos(qJ(1));
t298 = t257 * t263;
t212 = t256 * t295 + t298;
t294 = t263 * t256;
t213 = t257 * t295 - t294;
t259 = sin(qJ(2));
t297 = t259 * t260;
t315 = t212 * t320 + t213 * t319 - t297 * t318;
t214 = -t257 * t260 + t262 * t294;
t215 = t256 * t260 + t262 * t298;
t296 = t259 * t263;
t314 = t214 * t320 + t215 * t319 - t296 * t318;
t313 = t212 * t319 + t213 * t321 - t297 * t317;
t312 = t214 * t319 + t215 * t321 - t296 * t317;
t311 = -t212 * t318 - t213 * t317 - t297 * t316;
t310 = -t214 * t318 - t215 * t317 - t296 * t316;
t309 = t317 * t262 + (t256 * t319 + t257 * t321) * t259;
t308 = t316 * t262 + (-t256 * t318 - t257 * t317) * t259;
t307 = t318 * t262 + (t256 * t320 + t257 * t319) * t259;
t303 = Icges(2,4) * t260;
t302 = Icges(3,4) * t259;
t301 = Icges(3,4) * t262;
t300 = t256 * t259;
t299 = t257 * t259;
t171 = pkin(3) * t215 + qJ(4) * t214;
t279 = pkin(2) * t262 + qJ(3) * t259;
t220 = t279 * t263;
t293 = -t171 - t220;
t218 = (pkin(3) * t257 + qJ(4) * t256) * t259;
t238 = pkin(2) * t259 - qJ(3) * t262;
t292 = -t218 - t238;
t219 = t279 * t260;
t242 = pkin(1) * t260 - pkin(7) * t263;
t291 = -t219 - t242;
t290 = qJD(3) * t259;
t289 = qJD(6) * t259;
t288 = V_base(5) * pkin(6) + V_base(1);
t170 = pkin(3) * t213 + qJ(4) * t212;
t285 = -t170 + t291;
t178 = pkin(4) * t296 + qJ(5) * t215;
t284 = -t178 + t293;
t222 = -pkin(4) * t262 + qJ(5) * t299;
t283 = -t222 + t292;
t245 = qJD(2) * t260 + V_base(4);
t252 = V_base(6) + qJD(1);
t177 = pkin(4) * t297 + qJ(5) * t213;
t282 = -t177 + t285;
t244 = -qJD(2) * t263 + V_base(5);
t281 = t238 * t244 + t263 * t290 + t288;
t280 = rSges(3,1) * t262 - rSges(3,2) * t259;
t278 = Icges(3,1) * t262 - t302;
t277 = -Icges(3,2) * t259 + t301;
t276 = Icges(3,5) * t262 - Icges(3,6) * t259;
t243 = pkin(1) * t263 + pkin(7) * t260;
t275 = -V_base(4) * pkin(6) + t243 * t252 + V_base(2);
t274 = t242 * V_base(4) - t243 * V_base(5) + V_base(3);
t273 = qJD(4) * t214 + t218 * t244 + t281;
t272 = (-Icges(3,3) * t263 + t260 * t276) * t244 + (Icges(3,3) * t260 + t263 * t276) * t245 + (Icges(3,5) * t259 + Icges(3,6) * t262) * t252;
t271 = t220 * t252 + t260 * t290 + t275;
t270 = qJD(5) * t215 + t222 * t244 + t273;
t269 = -qJD(3) * t262 + t219 * t245 + t274;
t268 = qJD(4) * t212 + t171 * t252 + t271;
t267 = qJD(4) * t300 + t170 * t245 + t269;
t266 = qJD(5) * t213 + t178 * t252 + t268;
t265 = qJD(5) * t299 + t177 * t245 + t267;
t199 = -Icges(3,6) * t263 + t260 * t277;
t200 = Icges(3,6) * t260 + t263 * t277;
t201 = -Icges(3,5) * t263 + t260 * t278;
t202 = Icges(3,5) * t260 + t263 * t278;
t231 = Icges(3,2) * t262 + t302;
t234 = Icges(3,1) * t259 + t301;
t264 = (-t200 * t259 + t202 * t262) * t245 + (-t199 * t259 + t201 * t262) * t244 + (-t231 * t259 + t234 * t262) * t252;
t261 = cos(qJ(6));
t258 = sin(qJ(6));
t254 = Icges(2,4) * t263;
t241 = rSges(2,1) * t263 - rSges(2,2) * t260;
t240 = rSges(2,1) * t260 + rSges(2,2) * t263;
t239 = rSges(3,1) * t259 + rSges(3,2) * t262;
t237 = -qJD(6) * t262 + t252;
t236 = Icges(2,1) * t263 - t303;
t235 = Icges(2,1) * t260 + t254;
t233 = -Icges(2,2) * t260 + t254;
t232 = Icges(2,2) * t263 + t303;
t227 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t226 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t225 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t224 = pkin(5) * t300 - pkin(8) * t262;
t217 = t263 * t289 + t245;
t216 = t260 * t289 + t244;
t206 = (t256 * t261 + t257 * t258) * t259;
t205 = (-t256 * t258 + t257 * t261) * t259;
t204 = rSges(3,3) * t260 + t263 * t280;
t203 = -rSges(3,3) * t263 + t260 * t280;
t196 = -rSges(5,1) * t262 + (-rSges(5,2) * t257 + rSges(5,3) * t256) * t259;
t195 = -rSges(4,3) * t262 + (rSges(4,1) * t257 - rSges(4,2) * t256) * t259;
t194 = rSges(6,2) * t262 + (rSges(6,1) * t256 + rSges(6,3) * t257) * t259;
t180 = pkin(5) * t214 + pkin(8) * t296;
t179 = pkin(5) * t212 + pkin(8) * t297;
t176 = V_base(5) * rSges(2,3) - t240 * t252 + t288;
t175 = t241 * t252 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t174 = t240 * V_base(4) - t241 * V_base(5) + V_base(3);
t169 = t214 * t261 + t215 * t258;
t168 = -t214 * t258 + t215 * t261;
t167 = t212 * t261 + t213 * t258;
t166 = -t212 * t258 + t213 * t261;
t164 = rSges(4,1) * t215 - rSges(4,2) * t214 + rSges(4,3) * t296;
t163 = rSges(6,1) * t214 - rSges(6,2) * t296 + rSges(6,3) * t215;
t162 = rSges(4,1) * t213 - rSges(4,2) * t212 + rSges(4,3) * t297;
t161 = rSges(6,1) * t212 - rSges(6,2) * t297 + rSges(6,3) * t213;
t160 = rSges(5,1) * t296 - rSges(5,2) * t215 + rSges(5,3) * t214;
t159 = rSges(5,1) * t297 - rSges(5,2) * t213 + rSges(5,3) * t212;
t139 = rSges(7,1) * t206 + rSges(7,2) * t205 - rSges(7,3) * t262;
t138 = Icges(7,1) * t206 + Icges(7,4) * t205 - Icges(7,5) * t262;
t137 = Icges(7,4) * t206 + Icges(7,2) * t205 - Icges(7,6) * t262;
t136 = Icges(7,5) * t206 + Icges(7,6) * t205 - Icges(7,3) * t262;
t135 = t239 * t244 + (-t203 - t242) * t252 + t288;
t134 = t204 * t252 - t239 * t245 + t275;
t133 = t203 * t245 - t204 * t244 + t274;
t132 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t296;
t131 = rSges(7,1) * t167 + rSges(7,2) * t166 + rSges(7,3) * t297;
t130 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t296;
t129 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t297;
t128 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t296;
t127 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t297;
t126 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t296;
t125 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t297;
t124 = t195 * t244 + (-t162 + t291) * t252 + t281;
t123 = t164 * t252 + (-t195 - t238) * t245 + t271;
t122 = t162 * t245 + (-t164 - t220) * t244 + t269;
t121 = t196 * t244 + (-t159 + t285) * t252 + t273;
t120 = t160 * t252 + (-t196 + t292) * t245 + t268;
t119 = t159 * t245 + (-t160 + t293) * t244 + t267;
t118 = t194 * t244 + (-t161 + t282) * t252 + t270;
t117 = t163 * t252 + (-t194 + t283) * t245 + t266;
t116 = t161 * t245 + (-t163 + t284) * t244 + t265;
t115 = -t131 * t237 + t139 * t216 + t224 * t244 + (-t179 + t282) * t252 + t270;
t114 = t132 * t237 - t139 * t217 + t180 * t252 + (-t224 + t283) * t245 + t266;
t113 = t131 * t217 - t132 * t216 + t179 * t245 + (-t180 + t284) * t244 + t265;
t1 = t216 * ((t126 * t297 + t128 * t166 + t130 * t167) * t217 + (t125 * t297 + t166 * t127 + t167 * t129) * t216 + (t136 * t297 + t137 * t166 + t138 * t167) * t237) / 0.2e1 + t217 * ((t126 * t296 + t168 * t128 + t169 * t130) * t217 + (t125 * t296 + t127 * t168 + t129 * t169) * t216 + (t136 * t296 + t137 * t168 + t138 * t169) * t237) / 0.2e1 + m(1) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + m(2) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + m(3) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + t237 * ((-t126 * t262 + t128 * t205 + t130 * t206) * t217 + (-t125 * t262 + t127 * t205 + t129 * t206) * t216 + (-t262 * t136 + t205 * t137 + t206 * t138) * t237) / 0.2e1 + ((-t232 * t260 + t235 * t263 + Icges(1,4)) * V_base(5) + (-t260 * t233 + t263 * t236 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t263 * t232 + t260 * t235 + Icges(1,2)) * V_base(5) + (t233 * t263 + t236 * t260 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t260 * t264 - t263 * t272 + (t212 * t307 + t213 * t309 + t297 * t308) * t252 + (t212 * t314 + t213 * t312 + t297 * t310) * t245 + (t315 * t212 + t313 * t213 + t311 * t297) * t244) * t244 / 0.2e1 + (t260 * t272 + t263 * t264 + (t214 * t307 + t215 * t309 + t296 * t308) * t252 + (t314 * t214 + t312 * t215 + t310 * t296) * t245 + (t214 * t315 + t215 * t313 + t296 * t311) * t244) * t245 / 0.2e1 + (Icges(2,3) * t252 + ((t231 - t308) * t252 + (t200 - t310) * t245 + (t199 - t311) * t244) * t262 + ((t256 * t307 + t257 * t309 + t234) * t252 + (t256 * t314 + t257 * t312 + t202) * t245 + (t256 * t315 + t257 * t313 + t201) * t244) * t259) * t252 / 0.2e1 + t252 * V_base(4) * (Icges(2,5) * t263 - Icges(2,6) * t260) + t252 * V_base(5) * (Icges(2,5) * t260 + Icges(2,6) * t263) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
