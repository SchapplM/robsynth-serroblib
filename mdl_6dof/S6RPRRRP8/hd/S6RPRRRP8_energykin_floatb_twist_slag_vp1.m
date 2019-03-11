% Calculate kinetic energy for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:46
% EndTime: 2019-03-09 06:22:49
% DurationCPUTime: 2.80s
% Computational Cost: add. (1403->281), mult. (1766->394), div. (0->0), fcn. (1614->8), ass. (0->149)
t326 = Icges(2,4) + Icges(3,6);
t325 = Icges(2,1) + Icges(3,2);
t324 = Icges(6,1) + Icges(7,1);
t323 = -Icges(3,4) + Icges(2,5);
t322 = Icges(6,4) - Icges(7,5);
t321 = Icges(7,4) + Icges(6,5);
t320 = Icges(3,5) - Icges(2,6);
t319 = Icges(2,2) + Icges(3,3);
t318 = Icges(6,2) + Icges(7,3);
t317 = Icges(7,6) - Icges(6,6);
t316 = Icges(6,3) + Icges(7,2);
t315 = rSges(7,1) + pkin(5);
t314 = rSges(7,3) + qJ(6);
t233 = cos(qJ(1));
t313 = t326 * t233;
t230 = sin(qJ(1));
t312 = t326 * t230;
t227 = qJ(3) + qJ(4);
t222 = sin(t227);
t231 = cos(qJ(5));
t273 = t233 * t231;
t228 = sin(qJ(5));
t276 = t228 * t230;
t171 = t222 * t276 - t273;
t274 = t230 * t231;
t275 = t228 * t233;
t172 = t222 * t274 + t275;
t223 = cos(t227);
t278 = t223 * t230;
t311 = t171 * t318 - t172 * t322 - t278 * t317;
t173 = t222 * t275 + t274;
t174 = -t222 * t273 + t276;
t277 = t223 * t233;
t310 = -t173 * t318 - t174 * t322 + t277 * t317;
t309 = t171 * t317 + t172 * t321 - t278 * t316;
t308 = -t173 * t317 + t174 * t321 + t277 * t316;
t307 = -t171 * t322 + t172 * t324 - t278 * t321;
t306 = t173 * t322 + t174 * t324 + t277 * t321;
t305 = (t228 * t318 - t231 * t322) * t223 + t317 * t222;
t304 = (t228 * t317 + t231 * t321) * t223 + t316 * t222;
t303 = (-t228 * t322 + t231 * t324) * t223 + t321 * t222;
t302 = -t233 * t319 - t312;
t301 = t230 * t319 - t313;
t300 = t230 * t325 + t313;
t299 = t233 * t325 - t312;
t282 = Icges(5,4) * t222;
t254 = Icges(5,2) * t223 + t282;
t149 = Icges(5,6) * t233 + t230 * t254;
t150 = Icges(5,6) * t230 - t233 * t254;
t281 = Icges(5,4) * t223;
t256 = Icges(5,1) * t222 + t281;
t151 = Icges(5,5) * t233 + t230 * t256;
t152 = Icges(5,5) * t230 - t233 * t256;
t177 = -Icges(5,2) * t222 + t281;
t178 = Icges(5,1) * t223 - t282;
t212 = qJD(3) * t230 + V_base(5);
t181 = qJD(4) * t230 + t212;
t213 = qJD(3) * t233 + V_base(4);
t182 = qJD(4) * t233 + t213;
t216 = V_base(6) + qJD(1);
t296 = (t149 * t223 + t151 * t222) * t182 + (t150 * t223 + t152 * t222) * t181 + (t177 * t223 + t178 * t222) * t216;
t232 = cos(qJ(3));
t229 = sin(qJ(3));
t284 = Icges(4,4) * t229;
t255 = Icges(4,2) * t232 + t284;
t160 = Icges(4,6) * t233 + t230 * t255;
t161 = Icges(4,6) * t230 - t233 * t255;
t283 = Icges(4,4) * t232;
t257 = Icges(4,1) * t229 + t283;
t162 = Icges(4,5) * t233 + t230 * t257;
t163 = Icges(4,5) * t230 - t233 * t257;
t196 = -Icges(4,2) * t229 + t283;
t201 = Icges(4,1) * t232 - t284;
t295 = (t160 * t232 + t162 * t229) * t213 + (t161 * t232 + t163 * t229) * t212 + (t196 * t232 + t201 * t229) * t216;
t290 = pkin(3) * t229;
t289 = pkin(3) * t232;
t288 = t230 * pkin(7);
t287 = t233 * pkin(7);
t272 = -rSges(7,2) * t278 + t314 * t171 + t315 * t172;
t271 = rSges(7,2) * t277 - t314 * t173 + t315 * t174;
t270 = rSges(7,2) * t222 + (t314 * t228 + t315 * t231) * t223;
t269 = qJD(5) * t223;
t204 = pkin(1) * t230 - qJ(2) * t233;
t268 = V_base(4) * t204 + V_base(3);
t267 = V_base(5) * pkin(6) + V_base(1);
t264 = -t204 - t288;
t263 = qJD(2) * t230 + t267;
t169 = pkin(8) * t230 - t233 * t290;
t262 = -t169 + t264;
t261 = V_base(5) * pkin(2) + t263;
t260 = pkin(4) * t222 - pkin(9) * t223;
t259 = rSges(4,1) * t229 + rSges(4,2) * t232;
t258 = rSges(5,1) * t222 + rSges(5,2) * t223;
t253 = Icges(4,5) * t229 + Icges(4,6) * t232;
t252 = Icges(5,5) * t222 + Icges(5,6) * t223;
t208 = pkin(1) * t233 + qJ(2) * t230;
t245 = -qJD(2) * t233 + t216 * t208 + V_base(2);
t244 = t212 * t289 + t261;
t243 = (Icges(5,3) * t233 + t230 * t252) * t182 + (Icges(5,3) * t230 - t233 * t252) * t181 + (Icges(5,5) * t223 - Icges(5,6) * t222) * t216;
t242 = (Icges(4,3) * t233 + t230 * t253) * t213 + (Icges(4,3) * t230 - t233 * t253) * t212 + (Icges(4,5) * t232 - Icges(4,6) * t229) * t216;
t241 = V_base(4) * t288 + (-t208 - t287) * V_base(5) + t268;
t240 = t216 * t287 + (-pkin(2) - pkin(6)) * V_base(4) + t245;
t170 = pkin(8) * t233 + t230 * t290;
t239 = t213 * t169 - t170 * t212 + t241;
t168 = t260 * t233;
t180 = pkin(4) * t223 + pkin(9) * t222;
t238 = t181 * t180 + (t168 + t262) * t216 + t244;
t237 = t216 * t170 - t213 * t289 + t240;
t167 = t260 * t230;
t236 = -t167 * t181 - t182 * t168 + t239;
t235 = t216 * t167 - t180 * t182 + t237;
t210 = rSges(2,1) * t233 - rSges(2,2) * t230;
t209 = -rSges(3,2) * t233 + rSges(3,3) * t230;
t207 = rSges(4,1) * t232 - rSges(4,2) * t229;
t206 = rSges(2,1) * t230 + rSges(2,2) * t233;
t205 = -rSges(3,2) * t230 - rSges(3,3) * t233;
t188 = qJD(5) * t222 + t216;
t187 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t186 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t185 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t179 = rSges(5,1) * t223 - rSges(5,2) * t222;
t165 = rSges(4,3) * t230 - t233 * t259;
t164 = rSges(4,3) * t233 + t230 * t259;
t157 = -t230 * t269 + t182;
t156 = t233 * t269 + t181;
t154 = rSges(5,3) * t230 - t233 * t258;
t153 = rSges(5,3) * t233 + t230 * t258;
t145 = rSges(6,3) * t222 + (rSges(6,1) * t231 - rSges(6,2) * t228) * t223;
t137 = V_base(5) * rSges(2,3) - t206 * t216 + t267;
t136 = t210 * t216 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t133 = t206 * V_base(4) - t210 * V_base(5) + V_base(3);
t129 = V_base(5) * rSges(3,1) + (-t204 - t205) * t216 + t263;
t128 = t209 * t216 + (-rSges(3,1) - pkin(6)) * V_base(4) + t245;
t127 = rSges(6,1) * t174 + rSges(6,2) * t173 + rSges(6,3) * t277;
t125 = rSges(6,1) * t172 - rSges(6,2) * t171 - rSges(6,3) * t278;
t111 = t205 * V_base(4) + (-t208 - t209) * V_base(5) + t268;
t110 = t207 * t212 + (-t165 + t264) * t216 + t261;
t109 = t164 * t216 - t207 * t213 + t240;
t108 = -t164 * t212 + t165 * t213 + t241;
t107 = t179 * t181 + (-t154 + t262) * t216 + t244;
t106 = t153 * t216 - t179 * t182 + t237;
t105 = -t153 * t181 + t154 * t182 + t239;
t104 = -t127 * t188 + t145 * t156 + t238;
t103 = t125 * t188 - t145 * t157 + t235;
t102 = -t125 * t156 + t127 * t157 + t236;
t101 = qJD(6) * t171 + t156 * t270 - t188 * t271 + t238;
t100 = -qJD(6) * t173 - t157 * t270 + t188 * t272 + t235;
t99 = qJD(6) * t223 * t228 - t156 * t272 + t157 * t271 + t236;
t1 = t182 * (t296 * t230 + t243 * t233) / 0.2e1 + t181 * (t243 * t230 - t296 * t233) / 0.2e1 + t213 * (t295 * t230 + t242 * t233) / 0.2e1 + t212 * (t242 * t230 - t295 * t233) / 0.2e1 + m(1) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(2) * (t133 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(3) * (t111 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + ((-t173 * t305 + t174 * t303 + t277 * t304) * t188 + (-t173 * t311 + t307 * t174 + t309 * t277) * t157 + (-t310 * t173 + t306 * t174 + t308 * t277) * t156) * t156 / 0.2e1 + ((t171 * t305 + t172 * t303 - t278 * t304) * t188 + (t311 * t171 + t307 * t172 - t309 * t278) * t157 + (t171 * t310 + t172 * t306 - t278 * t308) * t156) * t157 / 0.2e1 + (((t228 * t305 + t231 * t303) * t188 + (t228 * t311 + t307 * t231) * t157 + (t228 * t310 + t231 * t306) * t156) * t223 + (t156 * t308 + t157 * t309 + t188 * t304) * t222) * t188 / 0.2e1 + ((t230 * t302 + t233 * t300 + Icges(1,4)) * V_base(5) + (t230 * t301 + t233 * t299 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t230 * t300 - t233 * t302 + Icges(1,2)) * V_base(5) + (t230 * t299 - t233 * t301 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t160 * t229 + t162 * t232) * t213 + (-t161 * t229 + t163 * t232) * t212 + (-t149 * t222 + t151 * t223) * t182 + (-t150 * t222 + t152 * t223) * t181 + (-t177 * t222 + t178 * t223 - t196 * t229 + t201 * t232 + Icges(3,1) + Icges(2,3)) * t216) * t216 / 0.2e1 + t216 * V_base(5) * (t230 * t323 - t233 * t320) + t216 * V_base(4) * (t230 * t320 + t233 * t323) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
