% Calculate kinetic energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:54
% EndTime: 2019-12-05 16:34:58
% DurationCPUTime: 3.97s
% Computational Cost: add. (2234->351), mult. (5356->503), div. (0->0), fcn. (6570->12), ass. (0->147)
t308 = Icges(4,2) + Icges(5,3);
t265 = sin(pkin(9));
t267 = cos(pkin(9));
t273 = cos(qJ(2));
t268 = cos(pkin(5));
t271 = sin(qJ(2));
t291 = t268 * t271;
t231 = t265 * t273 + t267 * t291;
t266 = sin(pkin(5));
t270 = sin(qJ(3));
t293 = t266 * t270;
t299 = cos(qJ(3));
t215 = t231 * t299 - t267 * t293;
t290 = t268 * t273;
t230 = t265 * t271 - t267 * t290;
t264 = sin(pkin(10));
t297 = cos(pkin(10));
t183 = t215 * t264 - t230 * t297;
t184 = t215 * t297 + t230 * t264;
t283 = t266 * t299;
t214 = t231 * t270 + t267 * t283;
t305 = -Icges(4,4) * t215 + Icges(5,5) * t184 - Icges(4,6) * t230 - Icges(5,6) * t183 + t308 * t214;
t233 = -t265 * t291 + t267 * t273;
t217 = t233 * t299 + t265 * t293;
t232 = t265 * t290 + t267 * t271;
t185 = t217 * t264 - t232 * t297;
t186 = t217 * t297 + t232 * t264;
t216 = t233 * t270 - t265 * t283;
t304 = -Icges(4,4) * t217 + Icges(5,5) * t186 - Icges(4,6) * t232 - Icges(5,6) * t185 + t308 * t216;
t238 = t268 * t270 + t271 * t283;
t292 = t266 * t273;
t212 = t238 * t264 + t292 * t297;
t213 = t238 * t297 - t264 * t292;
t237 = -t268 * t299 + t271 * t293;
t303 = -Icges(4,4) * t238 + Icges(5,5) * t213 + Icges(4,6) * t292 - Icges(5,6) * t212 + t308 * t237;
t298 = pkin(6) * t268;
t296 = Icges(2,4) * t265;
t295 = t265 * t266;
t294 = t266 * t267;
t289 = qJD(2) * t266;
t288 = V_base(5) * qJ(1) + V_base(1);
t284 = qJD(1) + V_base(3);
t246 = t265 * t289 + V_base(4);
t257 = qJD(2) * t268 + V_base(6);
t211 = qJD(3) * t232 + t246;
t245 = -t267 * t289 + V_base(5);
t210 = qJD(3) * t230 + t245;
t234 = -qJD(3) * t292 + t257;
t240 = pkin(1) * t265 - pkin(6) * t294;
t282 = -t240 * V_base(6) + V_base(5) * t298 + t288;
t241 = pkin(1) * t267 + pkin(6) * t295;
t281 = V_base(4) * t240 - t241 * V_base(5) + t284;
t280 = V_base(6) * t241 + V_base(2) + (-qJ(1) - t298) * V_base(4);
t204 = pkin(2) * t231 + pkin(7) * t230;
t239 = (pkin(2) * t271 - pkin(7) * t273) * t266;
t279 = -t204 * t257 + t245 * t239 + t282;
t205 = pkin(2) * t233 + pkin(7) * t232;
t278 = t246 * t204 - t205 * t245 + t281;
t277 = t257 * t205 - t239 * t246 + t280;
t206 = pkin(3) * t238 + qJ(4) * t237;
t276 = qJD(4) * t216 + t210 * t206 + t279;
t180 = pkin(3) * t215 + qJ(4) * t214;
t275 = qJD(4) * t237 + t211 * t180 + t278;
t181 = pkin(3) * t217 + qJ(4) * t216;
t274 = qJD(4) * t214 + t234 * t181 + t277;
t272 = cos(qJ(5));
t269 = sin(qJ(5));
t262 = Icges(2,4) * t267;
t254 = rSges(2,1) * t267 - rSges(2,2) * t265;
t253 = rSges(2,1) * t265 + rSges(2,2) * t267;
t252 = Icges(2,1) * t267 - t296;
t251 = Icges(2,1) * t265 + t262;
t250 = -Icges(2,2) * t265 + t262;
t249 = Icges(2,2) * t267 + t296;
t244 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t243 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t242 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t224 = t268 * rSges(3,3) + (rSges(3,1) * t271 + rSges(3,2) * t273) * t266;
t223 = Icges(3,5) * t268 + (Icges(3,1) * t271 + Icges(3,4) * t273) * t266;
t222 = Icges(3,6) * t268 + (Icges(3,4) * t271 + Icges(3,2) * t273) * t266;
t221 = Icges(3,3) * t268 + (Icges(3,5) * t271 + Icges(3,6) * t273) * t266;
t220 = V_base(5) * rSges(2,3) - t253 * V_base(6) + t288;
t219 = t254 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t209 = t253 * V_base(4) - t254 * V_base(5) + t284;
t203 = t238 * rSges(4,1) - t237 * rSges(4,2) - rSges(4,3) * t292;
t202 = Icges(4,1) * t238 - Icges(4,4) * t237 - Icges(4,5) * t292;
t200 = Icges(4,5) * t238 - Icges(4,6) * t237 - Icges(4,3) * t292;
t199 = rSges(3,1) * t233 - rSges(3,2) * t232 + rSges(3,3) * t295;
t198 = rSges(3,1) * t231 - rSges(3,2) * t230 - rSges(3,3) * t294;
t197 = Icges(3,1) * t233 - Icges(3,4) * t232 + Icges(3,5) * t295;
t196 = Icges(3,1) * t231 - Icges(3,4) * t230 - Icges(3,5) * t294;
t195 = Icges(3,4) * t233 - Icges(3,2) * t232 + Icges(3,6) * t295;
t194 = Icges(3,4) * t231 - Icges(3,2) * t230 - Icges(3,6) * t294;
t193 = Icges(3,5) * t233 - Icges(3,6) * t232 + Icges(3,3) * t295;
t192 = Icges(3,5) * t231 - Icges(3,6) * t230 - Icges(3,3) * t294;
t190 = qJD(5) * t212 + t234;
t188 = t213 * t272 + t237 * t269;
t187 = -t213 * t269 + t237 * t272;
t182 = pkin(4) * t213 + pkin(8) * t212;
t177 = rSges(5,1) * t213 - rSges(5,2) * t212 + rSges(5,3) * t237;
t176 = rSges(4,1) * t217 - rSges(4,2) * t216 + rSges(4,3) * t232;
t175 = rSges(4,1) * t215 - rSges(4,2) * t214 + rSges(4,3) * t230;
t174 = Icges(5,1) * t213 - Icges(5,4) * t212 + Icges(5,5) * t237;
t173 = Icges(5,4) * t213 - Icges(5,2) * t212 + Icges(5,6) * t237;
t171 = Icges(4,1) * t217 - Icges(4,4) * t216 + Icges(4,5) * t232;
t170 = Icges(4,1) * t215 - Icges(4,4) * t214 + Icges(4,5) * t230;
t167 = Icges(4,5) * t217 - Icges(4,6) * t216 + Icges(4,3) * t232;
t166 = Icges(4,5) * t215 - Icges(4,6) * t214 + Icges(4,3) * t230;
t165 = qJD(5) * t185 + t211;
t164 = qJD(5) * t183 + t210;
t163 = t186 * t272 + t216 * t269;
t162 = -t186 * t269 + t216 * t272;
t161 = t184 * t272 + t214 * t269;
t160 = -t184 * t269 + t214 * t272;
t158 = pkin(4) * t186 + pkin(8) * t185;
t157 = pkin(4) * t184 + pkin(8) * t183;
t156 = -t198 * t257 + t224 * t245 + t282;
t155 = t199 * t257 - t224 * t246 + t280;
t154 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t212;
t153 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t212;
t152 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t212;
t151 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t212;
t150 = rSges(5,1) * t186 - rSges(5,2) * t185 + rSges(5,3) * t216;
t149 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t214;
t148 = Icges(5,1) * t186 - Icges(5,4) * t185 + Icges(5,5) * t216;
t147 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t214;
t146 = Icges(5,4) * t186 - Icges(5,2) * t185 + Icges(5,6) * t216;
t145 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t214;
t142 = t198 * t246 - t199 * t245 + t281;
t141 = rSges(6,1) * t163 + rSges(6,2) * t162 + rSges(6,3) * t185;
t140 = rSges(6,1) * t161 + rSges(6,2) * t160 + rSges(6,3) * t183;
t139 = Icges(6,1) * t163 + Icges(6,4) * t162 + Icges(6,5) * t185;
t138 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t183;
t137 = Icges(6,4) * t163 + Icges(6,2) * t162 + Icges(6,6) * t185;
t136 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t183;
t135 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t185;
t134 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t183;
t133 = -t175 * t234 + t203 * t210 + t279;
t132 = t176 * t234 - t203 * t211 + t277;
t131 = t175 * t211 - t176 * t210 + t278;
t130 = t177 * t210 + (-t149 - t180) * t234 + t276;
t129 = t150 * t234 + (-t177 - t206) * t211 + t274;
t128 = t149 * t211 + (-t150 - t181) * t210 + t275;
t127 = -t140 * t190 + t154 * t164 + t182 * t210 + (-t157 - t180) * t234 + t276;
t126 = t141 * t190 - t154 * t165 + t158 * t234 + (-t182 - t206) * t211 + t274;
t125 = t140 * t165 - t141 * t164 + t157 * t211 + (-t158 - t181) * t210 + t275;
t1 = m(1) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(2) * (t209 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(3) * (t142 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + t246 * ((t193 * t295 - t195 * t232 + t197 * t233) * t246 + (t192 * t295 - t194 * t232 + t196 * t233) * t245 + (t221 * t295 - t222 * t232 + t223 * t233) * t257) / 0.2e1 + t245 * ((-t193 * t294 - t195 * t230 + t197 * t231) * t246 + (-t192 * t294 - t194 * t230 + t196 * t231) * t245 + (-t221 * t294 - t222 * t230 + t223 * t231) * t257) / 0.2e1 + t257 * ((t192 * t245 + t193 * t246 + t221 * t257) * t268 + ((t195 * t273 + t197 * t271) * t246 + (t194 * t273 + t196 * t271) * t245 + (t222 * t273 + t223 * t271) * t257) * t266) / 0.2e1 + m(4) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(5) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(6) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + t165 * ((t185 * t135 + t162 * t137 + t163 * t139) * t165 + (t134 * t185 + t136 * t162 + t138 * t163) * t164 + (t151 * t185 + t152 * t162 + t153 * t163) * t190) / 0.2e1 + t164 * ((t135 * t183 + t137 * t160 + t139 * t161) * t165 + (t183 * t134 + t160 * t136 + t161 * t138) * t164 + (t151 * t183 + t152 * t160 + t153 * t161) * t190) / 0.2e1 + t190 * ((t135 * t212 + t137 * t187 + t139 * t188) * t165 + (t134 * t212 + t136 * t187 + t138 * t188) * t164 + (t151 * t212 + t152 * t187 + t153 * t188) * t190) / 0.2e1 + ((-t173 * t183 + t174 * t184 + t200 * t230 + t202 * t215 + t214 * t303) * t234 + (-t146 * t183 + t148 * t184 + t167 * t230 + t171 * t215 + t214 * t304) * t211 + (-t145 * t183 + t147 * t184 + t166 * t230 + t170 * t215 + t305 * t214) * t210) * t210 / 0.2e1 + ((-t173 * t185 + t174 * t186 + t200 * t232 + t202 * t217 + t216 * t303) * t234 + (-t146 * t185 + t148 * t186 + t167 * t232 + t171 * t217 + t304 * t216) * t211 + (-t145 * t185 + t147 * t186 + t166 * t232 + t170 * t217 + t216 * t305) * t210) * t211 / 0.2e1 + ((-t212 * t173 + t213 * t174 - t200 * t292 + t238 * t202 + t237 * t303) * t234 + (-t146 * t212 + t148 * t213 - t167 * t292 + t238 * t171 + t237 * t304) * t211 + (-t145 * t212 + t147 * t213 - t166 * t292 + t238 * t170 + t237 * t305) * t210) * t234 / 0.2e1 + ((-t249 * t265 + t251 * t267 + Icges(1,4)) * V_base(5) + (-t250 * t265 + t252 * t267 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t249 * t267 + t251 * t265 + Icges(1,2)) * V_base(5) + (t250 * t267 + t252 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t265 + Icges(2,6) * t267 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t267 - Icges(2,6) * t265 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
