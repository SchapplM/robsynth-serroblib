% Calculate kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:29
% EndTime: 2019-03-09 08:33:32
% DurationCPUTime: 3.32s
% Computational Cost: add. (1043->263), mult. (1939->358), div. (0->0), fcn. (1777->6), ass. (0->138)
t340 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t339 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t338 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t235 = sin(qJ(2));
t337 = t340 * t235;
t238 = cos(qJ(2));
t336 = t340 * t238;
t335 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t334 = Icges(5,5) + Icges(3,6) - Icges(4,6);
t333 = t338 * t235 - t336;
t332 = t339 * t238 - t337;
t331 = Icges(6,1) + Icges(7,1);
t330 = Icges(6,4) + Icges(7,4);
t329 = Icges(7,5) + Icges(6,5);
t328 = Icges(6,2) + Icges(7,2);
t327 = Icges(7,6) + Icges(6,6);
t326 = Icges(7,3) + Icges(6,3);
t325 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t324 = t333 * t236 + t334 * t239;
t323 = -t334 * t236 + t333 * t239;
t322 = t332 * t236 - t335 * t239;
t321 = t335 * t236 + t332 * t239;
t320 = -t338 * t238 - t337;
t319 = t339 * t235 + t336;
t318 = -t334 * t235 + t335 * t238;
t237 = cos(qJ(5));
t286 = t237 * t239;
t234 = sin(qJ(5));
t289 = t236 * t234;
t182 = -t235 * t289 + t286;
t288 = t236 * t237;
t290 = t234 * t239;
t183 = t235 * t288 + t290;
t287 = t236 * t238;
t317 = t327 * t182 + t329 * t183 + t326 * t287;
t184 = -t235 * t290 - t288;
t185 = t235 * t286 - t289;
t285 = t238 * t239;
t316 = t327 * t184 + t329 * t185 + t326 * t285;
t315 = t328 * t182 + t330 * t183 + t327 * t287;
t314 = t328 * t184 + t330 * t185 + t327 * t285;
t313 = t330 * t182 + t331 * t183 + t329 * t287;
t312 = t330 * t184 + t331 * t185 + t329 * t285;
t311 = (t327 * t234 - t329 * t237) * t238 + t326 * t235;
t310 = (t328 * t234 - t330 * t237) * t238 + t327 * t235;
t309 = (t330 * t234 - t331 * t237) * t238 + t329 * t235;
t223 = -qJD(2) * t239 + V_base(5);
t224 = qJD(2) * t236 + V_base(4);
t228 = V_base(6) + qJD(1);
t308 = (t320 * t235 + t319 * t238) * t228 + (t323 * t235 + t321 * t238) * t224 + (t324 * t235 + t322 * t238) * t223;
t307 = (t335 * t235 + t334 * t238) * t228 + (t325 * t236 + t318 * t239) * t224 + (t318 * t236 - t325 * t239) * t223;
t300 = pkin(3) * t235;
t299 = pkin(5) * t237;
t297 = Icges(2,4) * t236;
t254 = qJ(6) * t238 + t235 * t299;
t284 = rSges(7,1) * t183 + rSges(7,2) * t182 + rSges(7,3) * t287 + pkin(5) * t290 + t236 * t254;
t283 = t185 * rSges(7,1) + t184 * rSges(7,2) + rSges(7,3) * t285 - pkin(5) * t289 + t239 * t254;
t282 = (-rSges(7,1) * t237 + rSges(7,2) * t234 - t299) * t238 + (qJ(6) + rSges(7,3)) * t235;
t266 = pkin(2) * t238 + qJ(3) * t235;
t186 = t266 * t236;
t220 = t236 * pkin(1) - pkin(7) * t239;
t281 = -t186 - t220;
t187 = t266 * t239;
t194 = pkin(3) * t285 - t236 * qJ(4);
t280 = -t187 - t194;
t279 = qJD(3) * t235;
t278 = qJD(5) * t238;
t277 = qJD(6) * t238;
t276 = V_base(5) * pkin(6) + V_base(1);
t193 = pkin(3) * t287 + qJ(4) * t239;
t273 = -t193 + t281;
t214 = pkin(2) * t235 - qJ(3) * t238;
t272 = -t214 - t300;
t271 = t223 * t214 + t239 * t279 + t276;
t270 = pkin(4) * t235 + pkin(8) * t238;
t269 = rSges(3,1) * t238 - rSges(3,2) * t235;
t268 = rSges(4,1) * t238 + rSges(4,3) * t235;
t267 = rSges(5,1) * t235 - rSges(5,2) * t238;
t222 = pkin(1) * t239 + t236 * pkin(7);
t256 = -V_base(4) * pkin(6) + t228 * t222 + V_base(2);
t255 = V_base(4) * t220 - t222 * V_base(5) + V_base(3);
t250 = t228 * t187 + t236 * t279 + t256;
t249 = -qJD(4) * t236 + t223 * t300 + t271;
t248 = -qJD(3) * t238 + t224 * t186 + t255;
t247 = qJD(4) * t239 + t228 * t194 + t250;
t246 = t224 * t193 + t248;
t189 = t270 * t236;
t221 = -pkin(4) * t238 + pkin(8) * t235;
t245 = t223 * t221 + (-t189 + t273) * t228 + t249;
t190 = t270 * t239;
t244 = t224 * t189 + (-t190 + t280) * t223 + t246;
t243 = t228 * t190 + (-t221 + t272) * t224 + t247;
t231 = Icges(2,4) * t239;
t219 = rSges(2,1) * t239 - t236 * rSges(2,2);
t218 = -rSges(5,1) * t238 - rSges(5,2) * t235;
t217 = t236 * rSges(2,1) + rSges(2,2) * t239;
t216 = rSges(3,1) * t235 + rSges(3,2) * t238;
t215 = rSges(4,1) * t235 - rSges(4,3) * t238;
t213 = qJD(5) * t235 + t228;
t212 = Icges(2,1) * t239 - t297;
t211 = Icges(2,1) * t236 + t231;
t207 = -Icges(2,2) * t236 + t231;
t206 = Icges(2,2) * t239 + t297;
t197 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t196 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t195 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t181 = t239 * t278 + t224;
t180 = t236 * t278 + t223;
t176 = t236 * rSges(3,3) + t239 * t269;
t175 = t236 * rSges(4,2) + t239 * t268;
t174 = -t236 * rSges(5,3) + t239 * t267;
t173 = rSges(6,3) * t235 + (-rSges(6,1) * t237 + rSges(6,2) * t234) * t238;
t171 = -rSges(3,3) * t239 + t236 * t269;
t170 = -rSges(4,2) * t239 + t236 * t268;
t169 = rSges(5,3) * t239 + t236 * t267;
t138 = V_base(5) * rSges(2,3) - t217 * t228 + t276;
t137 = t219 * t228 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t136 = t217 * V_base(4) - t219 * V_base(5) + V_base(3);
t135 = t185 * rSges(6,1) + t184 * rSges(6,2) + rSges(6,3) * t285;
t133 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t287;
t117 = t216 * t223 + (-t171 - t220) * t228 + t276;
t116 = t176 * t228 - t216 * t224 + t256;
t115 = t171 * t224 - t176 * t223 + t255;
t114 = t215 * t223 + (-t170 + t281) * t228 + t271;
t113 = t175 * t228 + (-t214 - t215) * t224 + t250;
t112 = t170 * t224 + (-t175 - t187) * t223 + t248;
t111 = t218 * t223 + (-t169 + t273) * t228 + t249;
t110 = t174 * t228 + (-t218 + t272) * t224 + t247;
t109 = t169 * t224 + (-t174 + t280) * t223 + t246;
t108 = -t133 * t213 + t173 * t180 + t245;
t107 = t135 * t213 - t173 * t181 + t243;
t106 = t133 * t181 - t135 * t180 + t244;
t105 = t180 * t282 - t213 * t284 + t239 * t277 + t245;
t104 = -t181 * t282 + t213 * t283 + t236 * t277 + t243;
t103 = qJD(6) * t235 - t180 * t283 + t181 * t284 + t244;
t1 = m(2) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(1) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + ((t182 * t310 + t183 * t309 + t287 * t311) * t213 + (t182 * t314 + t183 * t312 + t287 * t316) * t181 + (t315 * t182 + t313 * t183 + t317 * t287) * t180) * t180 / 0.2e1 + ((t184 * t310 + t185 * t309 + t285 * t311) * t213 + (t314 * t184 + t312 * t185 + t316 * t285) * t181 + (t315 * t184 + t313 * t185 + t285 * t317) * t180) * t181 / 0.2e1 + (((t234 * t310 - t237 * t309) * t213 + (t234 * t314 - t237 * t312) * t181 + (t234 * t315 - t237 * t313) * t180) * t238 + (t180 * t317 + t316 * t181 + t311 * t213) * t235) * t213 / 0.2e1 + ((-t236 * t206 + t211 * t239 + Icges(1,4)) * V_base(5) + (-t236 * t207 + t239 * t212 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t239 * t206 + t236 * t211 + Icges(1,2)) * V_base(5) + (t207 * t239 + t236 * t212 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t308 * t236 - t307 * t239) * t223 / 0.2e1 + (t307 * t236 + t308 * t239) * t224 / 0.2e1 + ((t321 * t235 - t323 * t238) * t224 + (t322 * t235 - t324 * t238) * t223 + (t319 * t235 - t320 * t238 + Icges(2,3)) * t228) * t228 / 0.2e1 + t228 * V_base(4) * (Icges(2,5) * t239 - Icges(2,6) * t236) + t228 * V_base(5) * (Icges(2,5) * t236 + Icges(2,6) * t239) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
