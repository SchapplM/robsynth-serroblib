% Calculate kinetic energy for
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:29
% EndTime: 2019-03-09 03:27:32
% DurationCPUTime: 3.18s
% Computational Cost: add. (1371->308), mult. (2001->424), div. (0->0), fcn. (1905->8), ass. (0->151)
t325 = Icges(2,4) + Icges(3,6);
t324 = Icges(2,1) + Icges(3,2);
t323 = Icges(6,1) + Icges(7,1);
t322 = -Icges(3,4) + Icges(2,5);
t321 = Icges(6,4) - Icges(7,5);
t320 = Icges(7,4) + Icges(6,5);
t319 = Icges(3,5) - Icges(2,6);
t318 = Icges(2,2) + Icges(3,3);
t317 = Icges(6,2) + Icges(7,3);
t316 = Icges(7,6) - Icges(6,6);
t315 = Icges(6,3) + Icges(7,2);
t314 = rSges(7,1) + pkin(5);
t313 = rSges(7,3) + qJ(6);
t238 = cos(qJ(1));
t312 = t325 * t238;
t236 = sin(qJ(1));
t311 = t325 * t236;
t235 = sin(qJ(3));
t231 = pkin(9) + qJ(5);
t222 = cos(t231);
t271 = t238 * t222;
t221 = sin(t231);
t277 = t236 * t221;
t163 = t235 * t277 - t271;
t276 = t236 * t222;
t164 = t221 * t238 + t235 * t276;
t237 = cos(qJ(3));
t273 = t236 * t237;
t310 = t317 * t163 - t321 * t164 - t316 * t273;
t278 = t235 * t238;
t165 = t221 * t278 + t276;
t166 = -t235 * t271 + t277;
t272 = t237 * t238;
t309 = -t317 * t165 - t321 * t166 + t316 * t272;
t308 = t316 * t163 + t320 * t164 - t315 * t273;
t307 = -t316 * t165 + t320 * t166 + t315 * t272;
t306 = -t321 * t163 + t323 * t164 - t320 * t273;
t305 = t321 * t165 + t323 * t166 + t320 * t272;
t304 = (t317 * t221 - t321 * t222) * t237 + t316 * t235;
t303 = (t316 * t221 + t320 * t222) * t237 + t315 * t235;
t302 = (-t321 * t221 + t323 * t222) * t237 + t320 * t235;
t301 = -t318 * t238 - t311;
t300 = t318 * t236 - t312;
t299 = t324 * t236 + t312;
t298 = t324 * t238 - t311;
t233 = cos(pkin(9));
t285 = pkin(4) * t233;
t295 = -pkin(8) * t237 + t235 * t285;
t283 = Icges(4,4) * t235;
t253 = Icges(4,2) * t237 + t283;
t169 = Icges(4,6) * t238 + t236 * t253;
t170 = Icges(4,6) * t236 - t238 * t253;
t282 = Icges(4,4) * t237;
t254 = Icges(4,1) * t235 + t282;
t171 = Icges(4,5) * t238 + t236 * t254;
t172 = Icges(4,5) * t236 - t238 * t254;
t197 = -Icges(4,2) * t235 + t282;
t202 = Icges(4,1) * t237 - t283;
t216 = qJD(3) * t236 + V_base(5);
t217 = qJD(3) * t238 + V_base(4);
t223 = V_base(6) + qJD(1);
t294 = (t169 * t237 + t171 * t235) * t217 + (t170 * t237 + t172 * t235) * t216 + (t197 * t237 + t202 * t235) * t223;
t287 = pkin(7) * t236;
t286 = pkin(7) * t238;
t232 = sin(pkin(9));
t279 = t232 * t238;
t275 = t236 * t232;
t274 = t236 * t233;
t269 = -rSges(7,2) * t273 + t313 * t163 + t164 * t314;
t268 = rSges(7,2) * t272 - t313 * t165 + t166 * t314;
t267 = rSges(7,2) * t235 + (t313 * t221 + t222 * t314) * t237;
t266 = qJD(4) * t237;
t265 = qJD(5) * t237;
t206 = t236 * pkin(1) - qJ(2) * t238;
t264 = V_base(4) * t206 + V_base(3);
t263 = V_base(5) * pkin(6) + V_base(1);
t260 = -t206 - t287;
t259 = qJD(2) * t236 + t263;
t255 = pkin(3) * t235 - qJ(4) * t237;
t185 = t255 * t238;
t258 = t185 + t260;
t257 = V_base(5) * pkin(2) + t259;
t256 = rSges(4,1) * t235 + rSges(4,2) * t237;
t252 = Icges(4,5) * t235 + Icges(4,6) * t237;
t211 = pkin(1) * t238 + t236 * qJ(2);
t248 = -qJD(2) * t238 + t223 * t211 + V_base(2);
t247 = (Icges(4,3) * t238 + t236 * t252) * t217 + (Icges(4,3) * t236 - t238 * t252) * t216 + (Icges(4,5) * t237 - Icges(4,6) * t235) * t223;
t209 = pkin(3) * t237 + qJ(4) * t235;
t246 = t216 * t209 - t236 * t266 + t257;
t245 = V_base(4) * t287 + (-t211 - t286) * V_base(5) + t264;
t244 = t223 * t286 + (-pkin(2) - pkin(6)) * V_base(4) + t248;
t243 = qJD(4) * t235 - t217 * t185 + t245;
t184 = t255 * t236;
t242 = t223 * t184 + t238 * t266 + t244;
t142 = pkin(4) * t279 + t236 * t295;
t143 = pkin(4) * t275 - t238 * t295;
t241 = t217 * t143 + (-t142 - t184) * t216 + t243;
t148 = pkin(8) * t235 + t237 * t285;
t240 = t216 * t148 + (-t143 + t258) * t223 + t246;
t239 = t223 * t142 + (-t148 - t209) * t217 + t242;
t213 = rSges(2,1) * t238 - t236 * rSges(2,2);
t212 = -rSges(3,2) * t238 + t236 * rSges(3,3);
t210 = rSges(4,1) * t237 - rSges(4,2) * t235;
t208 = t236 * rSges(2,1) + rSges(2,2) * t238;
t207 = -t236 * rSges(3,2) - rSges(3,3) * t238;
t205 = qJD(5) * t235 + t223;
t189 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t188 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t187 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t182 = -t236 * t265 + t217;
t181 = t238 * t265 + t216;
t180 = -t233 * t278 + t275;
t179 = t232 * t278 + t274;
t178 = t235 * t274 + t279;
t177 = t233 * t238 - t235 * t275;
t175 = t236 * rSges(4,3) - t238 * t256;
t174 = rSges(4,3) * t238 + t236 * t256;
t162 = rSges(5,3) * t235 + (rSges(5,1) * t233 - rSges(5,2) * t232) * t237;
t160 = Icges(5,5) * t235 + (Icges(5,1) * t233 - Icges(5,4) * t232) * t237;
t159 = Icges(5,6) * t235 + (Icges(5,4) * t233 - Icges(5,2) * t232) * t237;
t158 = Icges(5,3) * t235 + (Icges(5,5) * t233 - Icges(5,6) * t232) * t237;
t156 = rSges(6,3) * t235 + (rSges(6,1) * t222 - rSges(6,2) * t221) * t237;
t147 = V_base(5) * rSges(2,3) - t208 * t223 + t263;
t146 = t213 * t223 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t208 * V_base(4) - t213 * V_base(5) + V_base(3);
t141 = t180 * rSges(5,1) + t179 * rSges(5,2) + rSges(5,3) * t272;
t140 = rSges(5,1) * t178 + rSges(5,2) * t177 - rSges(5,3) * t273;
t139 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t272;
t138 = Icges(5,1) * t178 + Icges(5,4) * t177 - Icges(5,5) * t273;
t137 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t272;
t136 = Icges(5,4) * t178 + Icges(5,2) * t177 - Icges(5,6) * t273;
t135 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t272;
t134 = Icges(5,5) * t178 + Icges(5,6) * t177 - Icges(5,3) * t273;
t131 = V_base(5) * rSges(3,1) + (-t206 - t207) * t223 + t259;
t130 = t223 * t212 + (-rSges(3,1) - pkin(6)) * V_base(4) + t248;
t128 = t166 * rSges(6,1) + t165 * rSges(6,2) + rSges(6,3) * t272;
t126 = rSges(6,1) * t164 - rSges(6,2) * t163 - rSges(6,3) * t273;
t111 = t207 * V_base(4) + (-t211 - t212) * V_base(5) + t264;
t110 = t210 * t216 + (-t175 + t260) * t223 + t257;
t109 = t223 * t174 - t217 * t210 + t244;
t108 = -t216 * t174 + t217 * t175 + t245;
t107 = t162 * t216 + (-t141 + t258) * t223 + t246;
t106 = t223 * t140 + (-t162 - t209) * t217 + t242;
t105 = t217 * t141 + (-t140 - t184) * t216 + t243;
t104 = -t128 * t205 + t156 * t181 + t240;
t103 = t205 * t126 - t182 * t156 + t239;
t102 = -t181 * t126 + t182 * t128 + t241;
t101 = qJD(6) * t163 + t181 * t267 - t205 * t268 + t240;
t100 = -qJD(6) * t165 - t182 * t267 + t205 * t269 + t239;
t99 = qJD(6) * t237 * t221 - t181 * t269 + t182 * t268 + t241;
t1 = m(1) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(2) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t111 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + ((-t165 * t304 + t166 * t302 + t272 * t303) * t205 + (-t165 * t310 + t306 * t166 + t308 * t272) * t182 + (-t309 * t165 + t305 * t166 + t307 * t272) * t181) * t181 / 0.2e1 + ((t163 * t304 + t164 * t302 - t273 * t303) * t205 + (t310 * t163 + t306 * t164 - t308 * t273) * t182 + (t163 * t309 + t164 * t305 - t273 * t307) * t181) * t182 / 0.2e1 + (((t221 * t304 + t222 * t302) * t205 + (t221 * t310 + t306 * t222) * t182 + (t221 * t309 + t222 * t305) * t181) * t237 + (t181 * t307 + t182 * t308 + t205 * t303) * t235) * t205 / 0.2e1 + ((t134 * t272 + t179 * t136 + t180 * t138) * t217 + (t135 * t272 + t179 * t137 + t180 * t139) * t216 + (t158 * t272 + t179 * t159 + t180 * t160) * t223 + t247 * t236 - t294 * t238) * t216 / 0.2e1 + ((-t134 * t273 + t177 * t136 + t178 * t138) * t217 + (-t135 * t273 + t137 * t177 + t139 * t178) * t216 + (-t158 * t273 + t159 * t177 + t160 * t178) * t223 + t247 * t238 + t294 * t236) * t217 / 0.2e1 + ((t236 * t301 + t238 * t299 + Icges(1,4)) * V_base(5) + (t300 * t236 + t298 * t238 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t299 * t236 - t301 * t238 + Icges(1,2)) * V_base(5) + (t236 * t298 - t238 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t134 * t217 + t135 * t216) * t235 + ((-t136 * t232 + t138 * t233) * t217 + (-t137 * t232 + t139 * t233) * t216) * t237 + (-t169 * t235 + t171 * t237) * t217 + (-t170 * t235 + t172 * t237) * t216 + (Icges(2,3) + Icges(3,1) + (-t159 * t232 + t160 * t233 + t202) * t237 + (t158 - t197) * t235) * t223) * t223 / 0.2e1 + t223 * V_base(5) * (t322 * t236 - t319 * t238) + t223 * V_base(4) * (t319 * t236 + t322 * t238) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
