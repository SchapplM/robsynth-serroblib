% Calculate kinetic energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:49
% EndTime: 2019-03-09 06:02:52
% DurationCPUTime: 3.26s
% Computational Cost: add. (2229->317), mult. (2094->454), div. (0->0), fcn. (1996->10), ass. (0->156)
t324 = Icges(6,1) + Icges(7,1);
t323 = -Icges(6,4) + Icges(7,5);
t322 = Icges(7,4) + Icges(6,5);
t321 = Icges(6,2) + Icges(7,3);
t320 = -Icges(7,6) + Icges(6,6);
t319 = -Icges(6,3) - Icges(7,2);
t318 = rSges(7,1) + pkin(5);
t317 = rSges(7,3) + qJ(6);
t248 = qJ(1) + pkin(10);
t240 = sin(t248);
t241 = cos(t248);
t249 = qJ(4) + qJ(5);
t244 = cos(t249);
t243 = sin(t249);
t254 = cos(qJ(3));
t288 = t243 * t254;
t174 = t240 * t288 + t241 * t244;
t287 = t244 * t254;
t175 = t240 * t287 - t241 * t243;
t251 = sin(qJ(3));
t291 = t240 * t251;
t316 = t321 * t174 + t323 * t175 - t320 * t291;
t176 = -t240 * t244 + t241 * t288;
t177 = t240 * t243 + t241 * t287;
t289 = t241 * t251;
t315 = t321 * t176 + t323 * t177 - t320 * t289;
t314 = -t320 * t174 + t322 * t175 - t319 * t291;
t313 = -t320 * t176 + t322 * t177 - t319 * t289;
t312 = t323 * t174 + t324 * t175 + t322 * t291;
t311 = t323 * t176 + t324 * t177 + t322 * t289;
t310 = t320 * t254 + (t321 * t243 + t323 * t244) * t251;
t309 = t319 * t254 + (-t320 * t243 + t322 * t244) * t251;
t308 = -t322 * t254 + (t323 * t243 + t324 * t244) * t251;
t252 = sin(qJ(1));
t301 = pkin(1) * t252;
t255 = cos(qJ(1));
t300 = pkin(1) * t255;
t253 = cos(qJ(4));
t299 = pkin(4) * t253;
t298 = -pkin(6) - qJ(2);
t296 = Icges(2,4) * t252;
t295 = Icges(3,4) * t240;
t294 = Icges(4,4) * t251;
t293 = Icges(4,4) * t254;
t250 = sin(qJ(4));
t292 = t240 * t250;
t290 = t241 * t250;
t286 = t250 * t254;
t285 = t253 * t254;
t284 = rSges(7,2) * t291 + t317 * t174 + t318 * t175;
t283 = rSges(7,2) * t289 + t317 * t176 + t318 * t177;
t282 = -rSges(7,2) * t254 + (t317 * t243 + t318 * t244) * t251;
t281 = qJD(4) * t251;
t280 = qJD(5) * t251;
t242 = V_base(6) + qJD(1);
t279 = t242 * t300 + V_base(2);
t278 = V_base(5) * pkin(6) + V_base(1);
t217 = qJD(3) * t240 + V_base(4);
t211 = pkin(2) * t240 - pkin(7) * t241;
t275 = -t211 - t301;
t274 = V_base(5) * qJ(2) + t278;
t273 = V_base(4) * t301 + qJD(2) + V_base(3);
t188 = t241 * t281 + t217;
t272 = pkin(3) * t254 + pkin(8) * t251;
t216 = -qJD(3) * t241 + V_base(5);
t271 = rSges(4,1) * t254 - rSges(4,2) * t251;
t270 = Icges(4,1) * t254 - t294;
t269 = -Icges(4,2) * t251 + t293;
t268 = Icges(4,5) * t254 - Icges(4,6) * t251;
t187 = t240 * t281 + t216;
t267 = pkin(9) * t251 + t254 * t299;
t266 = (-Icges(4,3) * t241 + t240 * t268) * t216 + (Icges(4,3) * t240 + t241 * t268) * t217 + (Icges(4,5) * t251 + Icges(4,6) * t254) * t242;
t212 = pkin(2) * t241 + pkin(7) * t240;
t265 = t242 * t212 + t298 * V_base(4) + t279;
t198 = t272 * t240;
t234 = t251 * pkin(3) - t254 * pkin(8);
t264 = t216 * t234 + (-t198 + t275) * t242 + t274;
t263 = V_base(4) * t211 + (-t212 - t300) * V_base(5) + t273;
t199 = t272 * t241;
t262 = t242 * t199 - t217 * t234 + t265;
t150 = -pkin(4) * t290 + t240 * t267;
t173 = -pkin(9) * t254 + t251 * t299;
t230 = -qJD(4) * t254 + t242;
t261 = -t150 * t230 + t187 * t173 + t264;
t260 = t217 * t198 - t199 * t216 + t263;
t151 = pkin(4) * t292 + t241 * t267;
t259 = t230 * t151 - t173 * t188 + t262;
t258 = t188 * t150 - t151 * t187 + t260;
t164 = -Icges(4,6) * t241 + t240 * t269;
t165 = Icges(4,6) * t240 + t241 * t269;
t166 = -Icges(4,5) * t241 + t240 * t270;
t167 = Icges(4,5) * t240 + t241 * t270;
t222 = Icges(4,2) * t254 + t294;
t225 = Icges(4,1) * t251 + t293;
t257 = (-t165 * t251 + t167 * t254) * t217 + (-t164 * t251 + t166 * t254) * t216 + (-t222 * t251 + t225 * t254) * t242;
t246 = Icges(2,4) * t255;
t238 = Icges(3,4) * t241;
t233 = rSges(2,1) * t255 - rSges(2,2) * t252;
t232 = rSges(2,1) * t252 + rSges(2,2) * t255;
t231 = rSges(4,1) * t251 + rSges(4,2) * t254;
t227 = Icges(2,1) * t255 - t296;
t226 = Icges(2,1) * t252 + t246;
t224 = -Icges(2,2) * t252 + t246;
t223 = Icges(2,2) * t255 + t296;
t215 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t214 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t213 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t210 = rSges(3,1) * t241 - rSges(3,2) * t240;
t209 = rSges(3,1) * t240 + rSges(3,2) * t241;
t208 = (-qJD(4) - qJD(5)) * t254 + t242;
t207 = Icges(3,1) * t241 - t295;
t206 = Icges(3,1) * t240 + t238;
t205 = -Icges(3,2) * t240 + t238;
t204 = Icges(3,2) * t241 + t295;
t196 = -rSges(5,3) * t254 + (rSges(5,1) * t253 - rSges(5,2) * t250) * t251;
t195 = -Icges(5,5) * t254 + (Icges(5,1) * t253 - Icges(5,4) * t250) * t251;
t194 = -Icges(5,6) * t254 + (Icges(5,4) * t253 - Icges(5,2) * t250) * t251;
t193 = -Icges(5,3) * t254 + (Icges(5,5) * t253 - Icges(5,6) * t250) * t251;
t192 = t241 * t285 + t292;
t191 = t240 * t253 - t241 * t286;
t190 = t240 * t285 - t290;
t189 = -t240 * t286 - t241 * t253;
t185 = -rSges(6,3) * t254 + (rSges(6,1) * t244 - rSges(6,2) * t243) * t251;
t171 = rSges(4,3) * t240 + t241 * t271;
t170 = -rSges(4,3) * t241 + t240 * t271;
t169 = V_base(5) * rSges(2,3) - t232 * t242 + t278;
t168 = t233 * t242 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t161 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t160 = t241 * t280 + t188;
t159 = t240 * t280 + t187;
t157 = V_base(5) * rSges(3,3) + (-t209 - t301) * t242 + t274;
t156 = t210 * t242 + (-rSges(3,3) + t298) * V_base(4) + t279;
t152 = t209 * V_base(4) + (-t210 - t300) * V_base(5) + t273;
t149 = rSges(5,1) * t192 + rSges(5,2) * t191 + rSges(5,3) * t289;
t148 = rSges(5,1) * t190 + rSges(5,2) * t189 + rSges(5,3) * t291;
t147 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t289;
t146 = Icges(5,1) * t190 + Icges(5,4) * t189 + Icges(5,5) * t291;
t145 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t289;
t144 = Icges(5,4) * t190 + Icges(5,2) * t189 + Icges(5,6) * t291;
t143 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t289;
t142 = Icges(5,5) * t190 + Icges(5,6) * t189 + Icges(5,3) * t291;
t141 = rSges(6,1) * t177 - rSges(6,2) * t176 + rSges(6,3) * t289;
t139 = rSges(6,1) * t175 - rSges(6,2) * t174 + rSges(6,3) * t291;
t123 = t216 * t231 + (-t170 + t275) * t242 + t274;
t122 = t171 * t242 - t217 * t231 + t265;
t121 = t170 * t217 - t171 * t216 + t263;
t120 = -t148 * t230 + t187 * t196 + t264;
t119 = t149 * t230 - t188 * t196 + t262;
t118 = t148 * t188 - t149 * t187 + t260;
t117 = -t139 * t208 + t159 * t185 + t261;
t116 = t141 * t208 - t160 * t185 + t259;
t115 = t139 * t160 - t141 * t159 + t258;
t114 = qJD(6) * t176 + t159 * t282 - t208 * t284 + t261;
t113 = qJD(6) * t174 - t160 * t282 + t208 * t283 + t259;
t112 = qJD(6) * t243 * t251 - t159 * t283 + t160 * t284 + t258;
t1 = m(3) * (t152 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(7) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(2) * (t161 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + t217 * (t266 * t240 + t257 * t241) / 0.2e1 + t216 * (t257 * t240 - t266 * t241) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t187 * ((t143 * t291 + t145 * t189 + t147 * t190) * t188 + (t142 * t291 + t189 * t144 + t190 * t146) * t187 + (t189 * t194 + t190 * t195 + t193 * t291) * t230) / 0.2e1 + t188 * ((t143 * t289 + t191 * t145 + t192 * t147) * t188 + (t142 * t289 + t144 * t191 + t146 * t192) * t187 + (t191 * t194 + t192 * t195 + t193 * t289) * t230) / 0.2e1 + t230 * ((-t142 * t187 - t143 * t188 - t193 * t230) * t254 + ((-t145 * t250 + t147 * t253) * t188 + (-t144 * t250 + t146 * t253) * t187 + (-t194 * t250 + t195 * t253) * t230) * t251) / 0.2e1 + m(6) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + ((t174 * t310 + t175 * t308 + t291 * t309) * t208 + (t174 * t315 + t175 * t311 + t291 * t313) * t160 + (t316 * t174 + t312 * t175 + t314 * t291) * t159) * t159 / 0.2e1 + ((t176 * t310 + t177 * t308 + t289 * t309) * t208 + (t315 * t176 + t311 * t177 + t313 * t289) * t160 + (t176 * t316 + t312 * t177 + t314 * t289) * t159) * t160 / 0.2e1 + ((-t159 * t314 - t160 * t313 - t208 * t309) * t254 + ((t243 * t310 + t244 * t308) * t208 + (t243 * t315 + t244 * t311) * t160 + (t243 * t316 + t312 * t244) * t159) * t251) * t208 / 0.2e1 + ((t165 * t254 + t167 * t251) * t217 + (t164 * t254 + t166 * t251) * t216 + (t222 * t254 + t225 * t251 + Icges(2,3) + Icges(3,3)) * t242) * t242 / 0.2e1 + ((-t204 * t240 + t206 * t241 - t223 * t252 + t226 * t255 + Icges(1,4)) * V_base(5) + (-t205 * t240 + t207 * t241 - t224 * t252 + t227 * t255 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t204 * t241 + t206 * t240 + t223 * t255 + t226 * t252 + Icges(1,2)) * V_base(5) + (t205 * t241 + t207 * t240 + t224 * t255 + t227 * t252 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t242 * (Icges(2,5) * t255 + Icges(3,5) * t241 - Icges(2,6) * t252 - Icges(3,6) * t240) + V_base(5) * t242 * (Icges(2,5) * t252 + Icges(3,5) * t240 + Icges(2,6) * t255 + Icges(3,6) * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
