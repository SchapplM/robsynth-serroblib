% Calculate kinetic energy for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:02
% EndTime: 2019-12-05 16:54:05
% DurationCPUTime: 3.34s
% Computational Cost: add. (2000->309), mult. (4686->447), div. (0->0), fcn. (5615->10), ass. (0->140)
t303 = Icges(5,1) + Icges(6,1);
t302 = Icges(5,4) + Icges(6,4);
t301 = Icges(5,5) + Icges(6,5);
t300 = Icges(5,2) + Icges(6,2);
t299 = Icges(5,6) + Icges(6,6);
t298 = Icges(5,3) + Icges(6,3);
t297 = rSges(6,3) + qJ(5);
t240 = sin(pkin(9));
t242 = cos(pkin(9));
t249 = cos(qJ(2));
t243 = cos(pkin(5));
t247 = sin(qJ(2));
t271 = t243 * t247;
t206 = t240 * t249 + t242 * t271;
t241 = sin(pkin(5));
t246 = sin(qJ(3));
t273 = t241 * t246;
t282 = cos(qJ(3));
t191 = t206 * t282 - t242 * t273;
t270 = t243 * t249;
t205 = t240 * t247 - t242 * t270;
t245 = sin(qJ(4));
t248 = cos(qJ(4));
t165 = -t191 * t245 + t205 * t248;
t277 = t205 * t245;
t166 = t191 * t248 + t277;
t259 = t241 * t282;
t190 = t206 * t246 + t242 * t259;
t294 = t299 * t165 + t301 * t166 + t298 * t190;
t208 = -t240 * t271 + t242 * t249;
t193 = t208 * t282 + t240 * t273;
t207 = t240 * t270 + t242 * t247;
t167 = -t193 * t245 + t207 * t248;
t276 = t207 * t245;
t168 = t193 * t248 + t276;
t192 = t208 * t246 - t240 * t259;
t293 = t299 * t167 + t301 * t168 + t298 * t192;
t292 = t300 * t165 + t302 * t166 + t299 * t190;
t291 = t300 * t167 + t302 * t168 + t299 * t192;
t290 = t302 * t165 + t303 * t166 + t301 * t190;
t289 = t302 * t167 + t303 * t168 + t301 * t192;
t213 = t243 * t246 + t247 * t259;
t272 = t241 * t249;
t194 = -t213 * t245 - t248 * t272;
t260 = t245 * t272;
t195 = t213 * t248 - t260;
t212 = -t243 * t282 + t247 * t273;
t288 = t299 * t194 + t301 * t195 + t298 * t212;
t287 = t300 * t194 + t302 * t195 + t299 * t212;
t286 = t302 * t194 + t303 * t195 + t301 * t212;
t281 = pkin(6) * t243;
t280 = pkin(4) * t248;
t278 = Icges(2,4) * t240;
t275 = t240 * t241;
t274 = t241 * t242;
t269 = rSges(6,1) * t166 + rSges(6,2) * t165 + pkin(4) * t277 + t297 * t190 + t191 * t280;
t268 = rSges(6,1) * t168 + rSges(6,2) * t167 + pkin(4) * t276 + t297 * t192 + t193 * t280;
t267 = rSges(6,1) * t195 + rSges(6,2) * t194 - pkin(4) * t260 + t297 * t212 + t213 * t280;
t266 = qJD(2) * t241;
t265 = V_base(5) * qJ(1) + V_base(1);
t261 = qJD(1) + V_base(3);
t221 = t240 * t266 + V_base(4);
t232 = qJD(2) * t243 + V_base(6);
t189 = qJD(3) * t207 + t221;
t220 = -t242 * t266 + V_base(5);
t188 = qJD(3) * t205 + t220;
t209 = -qJD(3) * t272 + t232;
t215 = pkin(1) * t240 - pkin(6) * t274;
t258 = -t215 * V_base(6) + V_base(5) * t281 + t265;
t216 = pkin(1) * t242 + pkin(6) * t275;
t257 = V_base(4) * t215 - t216 * V_base(5) + t261;
t256 = V_base(6) * t216 + V_base(2) + (-qJ(1) - t281) * V_base(4);
t183 = pkin(2) * t206 + pkin(7) * t205;
t214 = (pkin(2) * t247 - pkin(7) * t249) * t241;
t255 = -t183 * t232 + t220 * t214 + t258;
t184 = pkin(2) * t208 + pkin(7) * t207;
t254 = t221 * t183 - t184 * t220 + t257;
t253 = t232 * t184 - t214 * t221 + t256;
t161 = pkin(3) * t191 + pkin(8) * t190;
t185 = t213 * pkin(3) + t212 * pkin(8);
t252 = -t161 * t209 + t188 * t185 + t255;
t162 = pkin(3) * t193 + pkin(8) * t192;
t251 = t189 * t161 - t162 * t188 + t254;
t250 = t209 * t162 - t185 * t189 + t253;
t238 = Icges(2,4) * t242;
t229 = rSges(2,1) * t242 - rSges(2,2) * t240;
t228 = rSges(2,1) * t240 + rSges(2,2) * t242;
t227 = Icges(2,1) * t242 - t278;
t226 = Icges(2,1) * t240 + t238;
t225 = -Icges(2,2) * t240 + t238;
t224 = Icges(2,2) * t242 + t278;
t219 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t218 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t217 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t202 = t243 * rSges(3,3) + (rSges(3,1) * t247 + rSges(3,2) * t249) * t241;
t201 = Icges(3,5) * t243 + (Icges(3,1) * t247 + Icges(3,4) * t249) * t241;
t200 = Icges(3,6) * t243 + (Icges(3,4) * t247 + Icges(3,2) * t249) * t241;
t199 = Icges(3,3) * t243 + (Icges(3,5) * t247 + Icges(3,6) * t249) * t241;
t198 = V_base(5) * rSges(2,3) - t228 * V_base(6) + t265;
t197 = t229 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t187 = t228 * V_base(4) - t229 * V_base(5) + t261;
t186 = qJD(4) * t212 + t209;
t182 = t213 * rSges(4,1) - t212 * rSges(4,2) - rSges(4,3) * t272;
t181 = Icges(4,1) * t213 - Icges(4,4) * t212 - Icges(4,5) * t272;
t180 = Icges(4,4) * t213 - Icges(4,2) * t212 - Icges(4,6) * t272;
t179 = Icges(4,5) * t213 - Icges(4,6) * t212 - Icges(4,3) * t272;
t178 = rSges(3,1) * t208 - rSges(3,2) * t207 + rSges(3,3) * t275;
t177 = rSges(3,1) * t206 - rSges(3,2) * t205 - rSges(3,3) * t274;
t176 = Icges(3,1) * t208 - Icges(3,4) * t207 + Icges(3,5) * t275;
t175 = Icges(3,1) * t206 - Icges(3,4) * t205 - Icges(3,5) * t274;
t174 = Icges(3,4) * t208 - Icges(3,2) * t207 + Icges(3,6) * t275;
t173 = Icges(3,4) * t206 - Icges(3,2) * t205 - Icges(3,6) * t274;
t172 = Icges(3,5) * t208 - Icges(3,6) * t207 + Icges(3,3) * t275;
t171 = Icges(3,5) * t206 - Icges(3,6) * t205 - Icges(3,3) * t274;
t164 = qJD(4) * t192 + t189;
t163 = qJD(4) * t190 + t188;
t159 = rSges(5,1) * t195 + rSges(5,2) * t194 + rSges(5,3) * t212;
t150 = rSges(4,1) * t193 - rSges(4,2) * t192 + rSges(4,3) * t207;
t149 = rSges(4,1) * t191 - rSges(4,2) * t190 + rSges(4,3) * t205;
t148 = Icges(4,1) * t193 - Icges(4,4) * t192 + Icges(4,5) * t207;
t147 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t205;
t146 = Icges(4,4) * t193 - Icges(4,2) * t192 + Icges(4,6) * t207;
t145 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t205;
t144 = Icges(4,5) * t193 - Icges(4,6) * t192 + Icges(4,3) * t207;
t143 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t205;
t140 = -t177 * t232 + t202 * t220 + t258;
t139 = t178 * t232 - t202 * t221 + t256;
t138 = rSges(5,1) * t168 + rSges(5,2) * t167 + rSges(5,3) * t192;
t136 = rSges(5,1) * t166 + rSges(5,2) * t165 + rSges(5,3) * t190;
t122 = t177 * t221 - t178 * t220 + t257;
t119 = -t149 * t209 + t182 * t188 + t255;
t118 = t150 * t209 - t182 * t189 + t253;
t117 = t149 * t189 - t150 * t188 + t254;
t116 = -t136 * t186 + t159 * t163 + t252;
t115 = t138 * t186 - t159 * t164 + t250;
t114 = t136 * t164 - t138 * t163 + t251;
t113 = qJD(5) * t192 + t163 * t267 - t186 * t269 + t252;
t112 = qJD(5) * t190 - t164 * t267 + t186 * t268 + t250;
t111 = qJD(5) * t212 - t163 * t268 + t164 * t269 + t251;
t1 = m(1) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + m(2) * (t187 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + t221 * ((t172 * t275 - t174 * t207 + t176 * t208) * t221 + (t171 * t275 - t173 * t207 + t175 * t208) * t220 + (t199 * t275 - t200 * t207 + t201 * t208) * t232) / 0.2e1 + t220 * ((-t172 * t274 - t174 * t205 + t176 * t206) * t221 + (-t171 * t274 - t173 * t205 + t175 * t206) * t220 + (-t199 * t274 - t200 * t205 + t201 * t206) * t232) / 0.2e1 + t232 * ((t171 * t220 + t172 * t221 + t199 * t232) * t243 + ((t174 * t249 + t176 * t247) * t221 + (t173 * t249 + t175 * t247) * t220 + (t200 * t249 + t201 * t247) * t232) * t241) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t189 * ((t144 * t207 - t146 * t192 + t148 * t193) * t189 + (t143 * t207 - t145 * t192 + t147 * t193) * t188 + (t179 * t207 - t180 * t192 + t181 * t193) * t209) / 0.2e1 + t188 * ((t144 * t205 - t146 * t190 + t148 * t191) * t189 + (t143 * t205 - t145 * t190 + t147 * t191) * t188 + (t179 * t205 - t180 * t190 + t181 * t191) * t209) / 0.2e1 + t209 * ((-t144 * t272 - t212 * t146 + t213 * t148) * t189 + (-t143 * t272 - t212 * t145 + t213 * t147) * t188 + (-t179 * t272 - t212 * t180 + t213 * t181) * t209) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + ((t165 * t287 + t166 * t286 + t190 * t288) * t186 + (t165 * t291 + t166 * t289 + t190 * t293) * t164 + (t292 * t165 + t290 * t166 + t294 * t190) * t163) * t163 / 0.2e1 + ((t167 * t287 + t168 * t286 + t192 * t288) * t186 + (t291 * t167 + t289 * t168 + t293 * t192) * t164 + (t167 * t292 + t168 * t290 + t192 * t294) * t163) * t164 / 0.2e1 + ((t287 * t194 + t286 * t195 + t288 * t212) * t186 + (t194 * t291 + t195 * t289 + t212 * t293) * t164 + (t194 * t292 + t195 * t290 + t212 * t294) * t163) * t186 / 0.2e1 + ((-t224 * t240 + t226 * t242 + Icges(1,4)) * V_base(5) + (-t225 * t240 + t227 * t242 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t224 * t242 + t226 * t240 + Icges(1,2)) * V_base(5) + (t225 * t242 + t227 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t240 + Icges(2,6) * t242 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t242 - Icges(2,6) * t240 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
