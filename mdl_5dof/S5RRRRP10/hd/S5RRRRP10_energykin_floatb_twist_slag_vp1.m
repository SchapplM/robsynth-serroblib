% Calculate kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:04
% EndTime: 2019-12-31 22:08:07
% DurationCPUTime: 3.18s
% Computational Cost: add. (2045->309), mult. (4686->449), div. (0->0), fcn. (5615->10), ass. (0->141)
t305 = Icges(5,1) + Icges(6,1);
t304 = Icges(5,4) + Icges(6,4);
t303 = Icges(5,5) + Icges(6,5);
t302 = Icges(5,2) + Icges(6,2);
t301 = Icges(5,6) + Icges(6,6);
t300 = Icges(5,3) + Icges(6,3);
t299 = rSges(6,3) + qJ(5);
t246 = cos(pkin(5));
t251 = sin(qJ(1));
t253 = cos(qJ(2));
t274 = t251 * t253;
t250 = sin(qJ(2));
t254 = cos(qJ(1));
t275 = t250 * t254;
t214 = t246 * t275 + t274;
t249 = sin(qJ(3));
t245 = sin(pkin(5));
t277 = t245 * t254;
t286 = cos(qJ(3));
t196 = t214 * t286 - t249 * t277;
t273 = t253 * t254;
t276 = t250 * t251;
t213 = -t246 * t273 + t276;
t248 = sin(qJ(4));
t252 = cos(qJ(4));
t169 = -t196 * t248 + t213 * t252;
t281 = t213 * t248;
t170 = t196 * t252 + t281;
t264 = t245 * t286;
t195 = t214 * t249 + t254 * t264;
t298 = t301 * t169 + t303 * t170 + t300 * t195;
t216 = -t246 * t276 + t273;
t279 = t245 * t251;
t198 = t216 * t286 + t249 * t279;
t215 = t246 * t274 + t275;
t171 = -t198 * t248 + t215 * t252;
t280 = t215 * t248;
t172 = t198 * t252 + t280;
t197 = t216 * t249 - t251 * t264;
t297 = t301 * t171 + t303 * t172 + t300 * t197;
t296 = t302 * t169 + t304 * t170 + t301 * t195;
t295 = t302 * t171 + t304 * t172 + t301 * t197;
t294 = t304 * t169 + t305 * t170 + t303 * t195;
t293 = t304 * t171 + t305 * t172 + t303 * t197;
t212 = t246 * t249 + t250 * t264;
t278 = t245 * t253;
t191 = -t212 * t248 - t252 * t278;
t265 = t248 * t278;
t192 = t212 * t252 - t265;
t211 = t245 * t249 * t250 - t246 * t286;
t292 = t301 * t191 + t303 * t192 + t300 * t211;
t291 = t302 * t191 + t304 * t192 + t301 * t211;
t290 = t304 * t191 + t305 * t192 + t303 * t211;
t285 = pkin(7) * t246;
t284 = pkin(4) * t252;
t282 = Icges(2,4) * t251;
t272 = rSges(6,1) * t170 + rSges(6,2) * t169 + pkin(4) * t281 + t299 * t195 + t196 * t284;
t271 = rSges(6,1) * t172 + rSges(6,2) * t171 + pkin(4) * t280 + t299 * t197 + t198 * t284;
t270 = rSges(6,1) * t192 + rSges(6,2) * t191 - pkin(4) * t265 + t299 * t211 + t212 * t284;
t269 = qJD(2) * t245;
t268 = V_base(5) * pkin(6) + V_base(1);
t225 = t251 * t269 + V_base(4);
t242 = V_base(6) + qJD(1);
t194 = qJD(3) * t215 + t225;
t226 = qJD(2) * t246 + t242;
t224 = -t254 * t269 + V_base(5);
t219 = t251 * pkin(1) - pkin(7) * t277;
t263 = -t219 * t242 + V_base(5) * t285 + t268;
t220 = pkin(1) * t254 + pkin(7) * t279;
t262 = V_base(4) * t219 - t220 * V_base(5) + V_base(3);
t193 = qJD(3) * t213 + t224;
t209 = -qJD(3) * t278 + t226;
t261 = t242 * t220 + V_base(2) + (-pkin(6) - t285) * V_base(4);
t189 = pkin(2) * t214 + pkin(8) * t213;
t218 = (pkin(2) * t250 - pkin(8) * t253) * t245;
t260 = -t189 * t226 + t224 * t218 + t263;
t190 = pkin(2) * t216 + pkin(8) * t215;
t259 = t225 * t189 - t190 * t224 + t262;
t258 = t226 * t190 - t218 * t225 + t261;
t165 = pkin(3) * t196 + pkin(9) * t195;
t187 = pkin(3) * t212 + pkin(9) * t211;
t257 = -t165 * t209 + t193 * t187 + t260;
t166 = pkin(3) * t198 + pkin(9) * t197;
t256 = t194 * t165 - t166 * t193 + t259;
t255 = t209 * t166 - t187 * t194 + t258;
t243 = Icges(2,4) * t254;
t234 = rSges(2,1) * t254 - t251 * rSges(2,2);
t233 = t251 * rSges(2,1) + rSges(2,2) * t254;
t232 = Icges(2,1) * t254 - t282;
t231 = Icges(2,1) * t251 + t243;
t230 = -Icges(2,2) * t251 + t243;
t229 = Icges(2,2) * t254 + t282;
t223 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t222 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t221 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t206 = rSges(3,3) * t246 + (rSges(3,1) * t250 + rSges(3,2) * t253) * t245;
t205 = Icges(3,5) * t246 + (Icges(3,1) * t250 + Icges(3,4) * t253) * t245;
t204 = Icges(3,6) * t246 + (Icges(3,4) * t250 + Icges(3,2) * t253) * t245;
t203 = Icges(3,3) * t246 + (Icges(3,5) * t250 + Icges(3,6) * t253) * t245;
t202 = V_base(5) * rSges(2,3) - t233 * t242 + t268;
t201 = t234 * t242 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t199 = t233 * V_base(4) - t234 * V_base(5) + V_base(3);
t188 = qJD(4) * t211 + t209;
t186 = rSges(3,1) * t216 - rSges(3,2) * t215 + rSges(3,3) * t279;
t185 = t214 * rSges(3,1) - t213 * rSges(3,2) - rSges(3,3) * t277;
t184 = Icges(3,1) * t216 - Icges(3,4) * t215 + Icges(3,5) * t279;
t183 = Icges(3,1) * t214 - Icges(3,4) * t213 - Icges(3,5) * t277;
t182 = Icges(3,4) * t216 - Icges(3,2) * t215 + Icges(3,6) * t279;
t181 = Icges(3,4) * t214 - Icges(3,2) * t213 - Icges(3,6) * t277;
t180 = Icges(3,5) * t216 - Icges(3,6) * t215 + Icges(3,3) * t279;
t179 = Icges(3,5) * t214 - Icges(3,6) * t213 - Icges(3,3) * t277;
t178 = rSges(4,1) * t212 - rSges(4,2) * t211 - rSges(4,3) * t278;
t177 = Icges(4,1) * t212 - Icges(4,4) * t211 - Icges(4,5) * t278;
t176 = Icges(4,4) * t212 - Icges(4,2) * t211 - Icges(4,6) * t278;
t175 = Icges(4,5) * t212 - Icges(4,6) * t211 - Icges(4,3) * t278;
t168 = qJD(4) * t197 + t194;
t167 = qJD(4) * t195 + t193;
t163 = rSges(4,1) * t198 - rSges(4,2) * t197 + rSges(4,3) * t215;
t162 = rSges(4,1) * t196 - rSges(4,2) * t195 + rSges(4,3) * t213;
t160 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t215;
t159 = Icges(4,1) * t196 - Icges(4,4) * t195 + Icges(4,5) * t213;
t158 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t215;
t157 = Icges(4,4) * t196 - Icges(4,2) * t195 + Icges(4,6) * t213;
t156 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t215;
t155 = Icges(4,5) * t196 - Icges(4,6) * t195 + Icges(4,3) * t213;
t154 = rSges(5,1) * t192 + rSges(5,2) * t191 + rSges(5,3) * t211;
t144 = -t185 * t226 + t206 * t224 + t263;
t143 = t186 * t226 - t206 * t225 + t261;
t142 = rSges(5,1) * t172 + rSges(5,2) * t171 + rSges(5,3) * t197;
t140 = rSges(5,1) * t170 + rSges(5,2) * t169 + rSges(5,3) * t195;
t126 = t185 * t225 - t186 * t224 + t262;
t123 = -t162 * t209 + t178 * t193 + t260;
t122 = t163 * t209 - t178 * t194 + t258;
t121 = t162 * t194 - t163 * t193 + t259;
t120 = -t140 * t188 + t154 * t167 + t257;
t119 = t142 * t188 - t154 * t168 + t255;
t118 = t140 * t168 - t142 * t167 + t256;
t117 = qJD(5) * t197 + t167 * t270 - t188 * t272 + t257;
t116 = qJD(5) * t195 - t168 * t270 + t188 * t271 + t255;
t115 = qJD(5) * t211 - t167 * t271 + t168 * t272 + t256;
t1 = m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t199 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(3) * (t126 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + t225 * ((t180 * t279 - t182 * t215 + t184 * t216) * t225 + (t179 * t279 - t181 * t215 + t183 * t216) * t224 + (t203 * t279 - t204 * t215 + t205 * t216) * t226) / 0.2e1 + t224 * ((-t180 * t277 - t213 * t182 + t214 * t184) * t225 + (-t179 * t277 - t213 * t181 + t214 * t183) * t224 + (-t203 * t277 - t213 * t204 + t214 * t205) * t226) / 0.2e1 + t226 * ((t179 * t224 + t180 * t225 + t203 * t226) * t246 + ((t182 * t253 + t184 * t250) * t225 + (t181 * t253 + t183 * t250) * t224 + (t204 * t253 + t205 * t250) * t226) * t245) / 0.2e1 + m(4) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + t194 * ((t156 * t215 - t158 * t197 + t160 * t198) * t194 + (t155 * t215 - t157 * t197 + t159 * t198) * t193 + (t175 * t215 - t176 * t197 + t177 * t198) * t209) / 0.2e1 + t193 * ((t156 * t213 - t158 * t195 + t160 * t196) * t194 + (t155 * t213 - t157 * t195 + t159 * t196) * t193 + (t175 * t213 - t176 * t195 + t177 * t196) * t209) / 0.2e1 + t209 * ((-t156 * t278 - t158 * t211 + t160 * t212) * t194 + (-t155 * t278 - t157 * t211 + t159 * t212) * t193 + (-t175 * t278 - t176 * t211 + t177 * t212) * t209) / 0.2e1 + m(5) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(6) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + ((t169 * t291 + t170 * t290 + t195 * t292) * t188 + (t169 * t295 + t170 * t293 + t195 * t297) * t168 + (t296 * t169 + t294 * t170 + t298 * t195) * t167) * t167 / 0.2e1 + ((t171 * t291 + t172 * t290 + t197 * t292) * t188 + (t295 * t171 + t293 * t172 + t297 * t197) * t168 + (t296 * t171 + t294 * t172 + t197 * t298) * t167) * t168 / 0.2e1 + ((t291 * t191 + t290 * t192 + t292 * t211) * t188 + (t191 * t295 + t192 * t293 + t211 * t297) * t168 + (t296 * t191 + t294 * t192 + t211 * t298) * t167) * t188 / 0.2e1 + ((-t251 * t229 + t231 * t254 + Icges(1,4)) * V_base(5) + (-t251 * t230 + t232 * t254 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t229 * t254 + t251 * t231 + Icges(1,2)) * V_base(5) + (t230 * t254 + t251 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t251 + Icges(2,6) * t254) * V_base(5) + (Icges(2,5) * t254 - Icges(2,6) * t251) * V_base(4) + Icges(2,3) * t242 / 0.2e1) * t242;
T = t1;
