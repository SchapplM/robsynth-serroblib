% Calculate kinetic energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [8x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S7RRRRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(6,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [7x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: m has to be [8x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [8,3]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [8x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [8 6]), ...
  'S7RRRRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:07
% EndTime: 2019-03-10 06:29:12
% DurationCPUTime: 5.21s
% Computational Cost: add. (4030->452), mult. (9598->694), div. (0->0), fcn. (12138->14), ass. (0->189)
t329 = cos(qJ(4));
t328 = cos(qJ(5));
t297 = cos(qJ(1));
t279 = -qJD(2) * t297 + V_base(5);
t327 = pkin(2) * t279;
t281 = V_base(6) + qJD(1);
t296 = cos(qJ(2));
t273 = qJD(3) * t296 + t281;
t290 = sin(qJ(3));
t291 = sin(qJ(2));
t321 = t290 * t291;
t249 = -qJD(4) * t321 + t273;
t326 = pkin(3) * t249;
t289 = sin(qJ(4));
t295 = cos(qJ(3));
t252 = t291 * t295 * t289 + t296 * t329;
t325 = pkin(3) * t252;
t292 = sin(qJ(1));
t324 = Icges(2,4) * t292;
t323 = Icges(3,4) * t291;
t322 = Icges(3,4) * t296;
t320 = t291 * t292;
t319 = t291 * t297;
t318 = t292 * t296;
t317 = t296 * t297;
t316 = qJD(3) * t291;
t315 = V_base(5) * pkin(1) + V_base(1);
t314 = pkin(2) * t320;
t311 = t291 * t329;
t280 = qJD(2) * t292 + V_base(4);
t310 = t281 * t314 + t296 * t327 + t315;
t309 = rSges(3,1) * t296 - rSges(3,2) * t291;
t308 = -V_base(4) * pkin(1) + V_base(2);
t307 = Icges(3,1) * t296 - t323;
t306 = -Icges(3,2) * t291 + t322;
t305 = Icges(3,5) * t296 - Icges(3,6) * t291;
t251 = -t297 * t316 + t280;
t304 = -t280 * t314 + t319 * t327 + V_base(3);
t303 = (-Icges(3,3) * t297 + t292 * t305) * t279 + (Icges(3,3) * t292 + t297 * t305) * t280 + (Icges(3,5) * t291 + Icges(3,6) * t296) * t281;
t256 = -t290 * t317 - t292 * t295;
t219 = qJD(4) * t256 + t251;
t250 = -t292 * t316 + t279;
t257 = -t292 * t290 + t295 * t317;
t228 = t257 * t289 - t297 * t311;
t190 = qJD(5) * t228 + t219;
t254 = -t290 * t318 + t295 * t297;
t218 = qJD(4) * t254 + t250;
t255 = t290 * t297 + t295 * t318;
t226 = t255 * t289 - t292 * t311;
t302 = t218 * t325 - t226 * t326 + t310;
t217 = qJD(5) * t252 + t249;
t229 = t257 * t329 + t289 * t319;
t288 = sin(qJ(5));
t199 = t229 * t288 - t256 * t328;
t169 = qJD(6) * t199 + t190;
t253 = -t296 * t289 + t295 * t311;
t224 = t253 * t288 + t321 * t328;
t188 = qJD(6) * t224 + t217;
t189 = qJD(5) * t226 + t218;
t301 = t304 + (-t218 * t228 + t219 * t226) * pkin(3);
t227 = t255 * t329 + t289 * t320;
t197 = t227 * t288 - t254 * t328;
t168 = qJD(6) * t197 + t189;
t300 = (-t280 * t296 - t281 * t319) * pkin(2) + t308;
t299 = -t219 * t325 + t228 * t326 + t300;
t236 = -Icges(3,6) * t297 + t292 * t306;
t237 = Icges(3,6) * t292 + t297 * t306;
t239 = -Icges(3,5) * t297 + t292 * t307;
t240 = Icges(3,5) * t292 + t297 * t307;
t267 = Icges(3,2) * t296 + t323;
t270 = Icges(3,1) * t291 + t322;
t298 = (-t237 * t291 + t240 * t296) * t280 + (-t236 * t291 + t239 * t296) * t279 + (-t267 * t291 + t270 * t296) * t281;
t294 = cos(qJ(6));
t293 = cos(qJ(7));
t287 = sin(qJ(6));
t286 = sin(qJ(7));
t284 = Icges(2,4) * t297;
t276 = rSges(2,1) * t297 - t292 * rSges(2,2);
t275 = t292 * rSges(2,1) + rSges(2,2) * t297;
t274 = rSges(3,1) * t291 + rSges(3,2) * t296;
t272 = Icges(2,1) * t297 - t324;
t271 = Icges(2,1) * t292 + t284;
t269 = -Icges(2,2) * t292 + t284;
t268 = Icges(2,2) * t297 + t324;
t263 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t262 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t261 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t243 = t292 * rSges(3,3) + t297 * t309;
t242 = -rSges(3,3) * t297 + t292 * t309;
t241 = rSges(4,3) * t296 + (rSges(4,1) * t295 - rSges(4,2) * t290) * t291;
t238 = Icges(4,5) * t296 + (Icges(4,1) * t295 - Icges(4,4) * t290) * t291;
t235 = Icges(4,6) * t296 + (Icges(4,4) * t295 - Icges(4,2) * t290) * t291;
t232 = Icges(4,3) * t296 + (Icges(4,5) * t295 - Icges(4,6) * t290) * t291;
t231 = V_base(5) * rSges(2,3) - t275 * t281 + t315;
t230 = t276 * t281 + V_base(2) + (-rSges(2,3) - pkin(1)) * V_base(4);
t225 = t253 * t328 - t288 * t321;
t223 = t275 * V_base(4) - t276 * V_base(5) + V_base(3);
t216 = t257 * rSges(4,1) + t256 * rSges(4,2) - rSges(4,3) * t319;
t215 = rSges(4,1) * t255 + rSges(4,2) * t254 - rSges(4,3) * t320;
t214 = rSges(5,1) * t253 - rSges(5,2) * t252 - rSges(5,3) * t321;
t213 = Icges(4,1) * t257 + Icges(4,4) * t256 - Icges(4,5) * t319;
t212 = Icges(4,1) * t255 + Icges(4,4) * t254 - Icges(4,5) * t320;
t211 = Icges(5,1) * t253 - Icges(5,4) * t252 - Icges(5,5) * t321;
t210 = Icges(4,4) * t257 + Icges(4,2) * t256 - Icges(4,6) * t319;
t209 = Icges(4,4) * t255 + Icges(4,2) * t254 - Icges(4,6) * t320;
t208 = Icges(5,4) * t253 - Icges(5,2) * t252 - Icges(5,6) * t321;
t207 = Icges(4,5) * t257 + Icges(4,6) * t256 - Icges(4,3) * t319;
t206 = Icges(4,5) * t255 + Icges(4,6) * t254 - Icges(4,3) * t320;
t205 = Icges(5,5) * t253 - Icges(5,6) * t252 - Icges(5,3) * t321;
t202 = -t242 * t281 + t274 * t279 + t315;
t201 = t243 * t281 - t274 * t280 + t308;
t200 = t229 * t328 + t256 * t288;
t198 = t227 * t328 + t254 * t288;
t196 = t225 * t294 + t252 * t287;
t195 = -t225 * t287 + t252 * t294;
t192 = t242 * t280 - t243 * t279 + V_base(3);
t187 = rSges(5,1) * t229 - rSges(5,2) * t228 + rSges(5,3) * t256;
t186 = rSges(5,1) * t227 - rSges(5,2) * t226 + rSges(5,3) * t254;
t185 = rSges(6,1) * t225 - rSges(6,2) * t224 + rSges(6,3) * t252;
t184 = Icges(5,1) * t229 - Icges(5,4) * t228 + Icges(5,5) * t256;
t183 = Icges(5,1) * t227 - Icges(5,4) * t226 + Icges(5,5) * t254;
t182 = Icges(6,1) * t225 - Icges(6,4) * t224 + Icges(6,5) * t252;
t181 = Icges(5,4) * t229 - Icges(5,2) * t228 + Icges(5,6) * t256;
t180 = Icges(5,4) * t227 - Icges(5,2) * t226 + Icges(5,6) * t254;
t179 = Icges(6,4) * t225 - Icges(6,2) * t224 + Icges(6,6) * t252;
t178 = Icges(5,5) * t229 - Icges(5,6) * t228 + Icges(5,3) * t256;
t177 = Icges(5,5) * t227 - Icges(5,6) * t226 + Icges(5,3) * t254;
t176 = Icges(6,5) * t225 - Icges(6,6) * t224 + Icges(6,3) * t252;
t175 = t200 * t294 + t228 * t287;
t174 = -t200 * t287 + t228 * t294;
t173 = t198 * t294 + t226 * t287;
t172 = -t198 * t287 + t226 * t294;
t171 = t196 * t293 - t224 * t286;
t170 = -t196 * t286 - t224 * t293;
t167 = qJD(7) * t195 + t188;
t166 = -t215 * t273 + t241 * t250 + t310;
t165 = t273 * t216 - t251 * t241 + t300;
t164 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t228;
t163 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t226;
t162 = rSges(7,1) * t196 + rSges(7,2) * t195 + rSges(7,3) * t224;
t161 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t228;
t160 = Icges(6,1) * t198 - Icges(6,4) * t197 + Icges(6,5) * t226;
t159 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t224;
t158 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t228;
t157 = Icges(6,4) * t198 - Icges(6,2) * t197 + Icges(6,6) * t226;
t156 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t224;
t155 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t228;
t154 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t226;
t153 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t224;
t152 = t175 * t293 - t199 * t286;
t151 = -t175 * t286 - t199 * t293;
t150 = t173 * t293 - t197 * t286;
t149 = -t173 * t286 - t197 * t293;
t148 = t215 * t251 - t216 * t250 + t304;
t147 = -t186 * t249 + t214 * t218 + t310;
t146 = t249 * t187 - t219 * t214 + t300;
t145 = qJD(7) * t174 + t169;
t144 = qJD(7) * t172 + t168;
t143 = rSges(7,1) * t175 + rSges(7,2) * t174 + rSges(7,3) * t199;
t142 = rSges(7,1) * t173 + rSges(7,2) * t172 + rSges(7,3) * t197;
t141 = rSges(8,1) * t171 + rSges(8,2) * t170 + rSges(8,3) * t195;
t140 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t199;
t139 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t197;
t138 = Icges(8,1) * t171 + Icges(8,4) * t170 + Icges(8,5) * t195;
t137 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t199;
t136 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t197;
t135 = Icges(8,4) * t171 + Icges(8,2) * t170 + Icges(8,6) * t195;
t134 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t199;
t133 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t197;
t132 = Icges(8,5) * t171 + Icges(8,6) * t170 + Icges(8,3) * t195;
t131 = t186 * t219 - t187 * t218 + t304;
t130 = rSges(8,1) * t152 + rSges(8,2) * t151 + rSges(8,3) * t174;
t129 = rSges(8,1) * t150 + rSges(8,2) * t149 + rSges(8,3) * t172;
t128 = Icges(8,1) * t152 + Icges(8,4) * t151 + Icges(8,5) * t174;
t127 = Icges(8,1) * t150 + Icges(8,4) * t149 + Icges(8,5) * t172;
t126 = Icges(8,4) * t152 + Icges(8,2) * t151 + Icges(8,6) * t174;
t125 = Icges(8,4) * t150 + Icges(8,2) * t149 + Icges(8,6) * t172;
t124 = Icges(8,5) * t152 + Icges(8,6) * t151 + Icges(8,3) * t174;
t123 = Icges(8,5) * t150 + Icges(8,6) * t149 + Icges(8,3) * t172;
t122 = -t163 * t217 + t185 * t189 + t302;
t121 = t217 * t164 - t190 * t185 + t299;
t120 = t163 * t190 - t164 * t189 + t301;
t119 = -t142 * t188 + t162 * t168 + t302;
t118 = t188 * t143 - t169 * t162 + t299;
t117 = t142 * t169 - t143 * t168 + t301;
t116 = -t129 * t167 + t141 * t144 + (t168 * t195 - t172 * t188) * pkin(4) + t302;
t115 = t167 * t130 - t145 * t141 + (-t169 * t195 + t174 * t188) * pkin(4) + t299;
t114 = t129 * t145 - t130 * t144 + (-t168 * t174 + t169 * t172) * pkin(4) + t301;
t1 = ((-t292 * t268 + t271 * t297 + Icges(1,4)) * V_base(5) + (-t292 * t269 + t272 * t297 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t219 * ((t178 * t256 - t181 * t228 + t184 * t229) * t219 + (t177 * t256 - t180 * t228 + t183 * t229) * t218 + (t205 * t256 - t208 * t228 + t211 * t229) * t249) / 0.2e1 + m(1) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t217 * ((t155 * t252 - t158 * t224 + t161 * t225) * t190 + (t154 * t252 - t157 * t224 + t160 * t225) * t189 + (t176 * t252 - t179 * t224 + t182 * t225) * t217) / 0.2e1 + t218 * ((t178 * t254 - t181 * t226 + t184 * t227) * t219 + (t177 * t254 - t180 * t226 + t183 * t227) * t218 + (t205 * t254 - t208 * t226 + t211 * t227) * t249) / 0.2e1 + t190 * ((t155 * t228 - t158 * t199 + t161 * t200) * t190 + (t154 * t228 - t157 * t199 + t160 * t200) * t189 + (t176 * t228 - t179 * t199 + t182 * t200) * t217) / 0.2e1 + m(2) * (t223 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t188 * ((t134 * t224 + t137 * t195 + t140 * t196) * t169 + (t133 * t224 + t136 * t195 + t139 * t196) * t168 + (t153 * t224 + t156 * t195 + t159 * t196) * t188) / 0.2e1 + t189 * ((t155 * t226 - t158 * t197 + t161 * t198) * t190 + (t154 * t226 - t157 * t197 + t160 * t198) * t189 + (t176 * t226 - t179 * t197 + t182 * t198) * t217) / 0.2e1 + t169 * ((t134 * t199 + t137 * t174 + t140 * t175) * t169 + (t133 * t199 + t136 * t174 + t139 * t175) * t168 + (t153 * t199 + t156 * t174 + t159 * t175) * t188) / 0.2e1 + m(3) * (t192 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + t167 * ((t124 * t195 + t126 * t170 + t128 * t171) * t145 + (t123 * t195 + t125 * t170 + t127 * t171) * t144 + (t132 * t195 + t135 * t170 + t138 * t171) * t167) / 0.2e1 + t168 * ((t134 * t197 + t137 * t172 + t140 * t173) * t169 + (t133 * t197 + t136 * t172 + t139 * t173) * t168 + (t153 * t197 + t156 * t172 + t159 * t173) * t188) / 0.2e1 + t144 * ((t124 * t172 + t126 * t149 + t128 * t150) * t145 + (t172 * t123 + t149 * t125 + t150 * t127) * t144 + (t132 * t172 + t135 * t149 + t138 * t150) * t167) / 0.2e1 + t145 * ((t174 * t124 + t151 * t126 + t152 * t128) * t145 + (t123 * t174 + t125 * t151 + t127 * t152) * t144 + (t132 * t174 + t135 * t151 + t138 * t152) * t167) / 0.2e1 + m(4) * (t148 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(5) * (t131 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(6) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(8) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + V_base(5) * t281 * (Icges(2,5) * t292 + Icges(2,6) * t297) + t273 * ((t206 * t250 + t207 * t251 + t232 * t273) * t296 + ((-t210 * t290 + t213 * t295) * t251 + (-t209 * t290 + t212 * t295) * t250 + (-t235 * t290 + t238 * t295) * t273) * t291) / 0.2e1 + ((t237 * t296 + t240 * t291) * t280 + (t236 * t296 + t239 * t291) * t279 + (t267 * t296 + t270 * t291 + Icges(2,3)) * t281) * t281 / 0.2e1 + ((t268 * t297 + t292 * t271 + Icges(1,2)) * V_base(5) + (t269 * t297 + t292 * t272 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t281 * V_base(4) * (Icges(2,5) * t297 - Icges(2,6) * t292) + t249 * ((-t178 * t321 - t181 * t252 + t184 * t253) * t219 + (-t177 * t321 - t180 * t252 + t183 * t253) * t218 + (-t205 * t321 - t208 * t252 + t211 * t253) * t249) / 0.2e1 + t251 * ((-t207 * t319 + t256 * t210 + t257 * t213) * t251 + (-t206 * t319 + t256 * t209 + t257 * t212) * t250 + (-t232 * t319 + t256 * t235 + t257 * t238) * t273) / 0.2e1 + t250 * ((-t207 * t320 + t210 * t254 + t213 * t255) * t251 + (-t206 * t320 + t209 * t254 + t212 * t255) * t250 + (-t232 * t320 + t235 * t254 + t238 * t255) * t273) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t280 * (t292 * t303 + t297 * t298) / 0.2e1 + t279 * (t292 * t298 - t303 * t297) / 0.2e1;
T  = t1;
