% Calculate kinetic energy for
% S5RRPRR10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:42
% EndTime: 2019-12-31 20:23:45
% DurationCPUTime: 3.35s
% Computational Cost: add. (2456->357), mult. (5755->520), div. (0->0), fcn. (7113->12), ass. (0->160)
t275 = sin(qJ(1));
t278 = cos(qJ(1));
t274 = sin(qJ(2));
t277 = cos(qJ(2));
t306 = sin(pkin(10));
t307 = cos(pkin(10));
t240 = -t274 * t306 + t277 * t307;
t271 = cos(pkin(5));
t283 = t271 * t240;
t287 = t274 * t307 + t277 * t306;
t212 = -t275 * t287 + t278 * t283;
t231 = t287 * t271;
t213 = t231 * t278 + t275 * t240;
t270 = sin(pkin(5));
t303 = t270 * t278;
t168 = Icges(4,5) * t213 + Icges(4,6) * t212 - Icges(4,3) * t303;
t299 = t277 * t278;
t301 = t275 * t274;
t234 = t271 * t299 - t301;
t300 = t275 * t277;
t302 = t274 * t278;
t235 = t271 * t302 + t300;
t199 = Icges(3,5) * t235 + Icges(3,6) * t234 - Icges(3,3) * t303;
t317 = t168 + t199;
t214 = -t275 * t283 - t278 * t287;
t215 = -t275 * t231 + t240 * t278;
t304 = t270 * t275;
t169 = Icges(4,5) * t215 + Icges(4,6) * t214 + Icges(4,3) * t304;
t236 = -t271 * t300 - t302;
t237 = -t271 * t301 + t299;
t200 = Icges(3,5) * t237 + Icges(3,6) * t236 + Icges(3,3) * t304;
t316 = t169 + t200;
t229 = t240 * t270;
t230 = t287 * t270;
t195 = Icges(4,5) * t230 + Icges(4,6) * t229 + Icges(4,3) * t271;
t225 = Icges(3,3) * t271 + (Icges(3,5) * t274 + Icges(3,6) * t277) * t270;
t315 = t195 + t225;
t311 = cos(qJ(4));
t310 = pkin(2) * t274;
t309 = pkin(7) * t271;
t308 = pkin(2) * t277;
t305 = Icges(2,4) * t275;
t298 = qJD(2) * t270;
t297 = qJD(3) * t270;
t296 = V_base(5) * pkin(6) + V_base(1);
t293 = t270 * t311;
t248 = t275 * t298 + V_base(4);
t267 = V_base(6) + qJD(1);
t292 = -qJ(3) * t270 + t271 * t310;
t189 = -qJD(4) * t214 + t248;
t249 = qJD(2) * t271 + t267;
t216 = -qJD(4) * t229 + t249;
t247 = -t278 * t298 + V_base(5);
t242 = t275 * pkin(1) - pkin(7) * t303;
t289 = -t242 * t267 + V_base(5) * t309 + t296;
t243 = pkin(1) * t278 + pkin(7) * t304;
t288 = V_base(4) * t242 - t243 * V_base(5) + V_base(3);
t188 = -qJD(4) * t212 + t247;
t241 = qJ(3) * t271 + t270 * t310;
t286 = t247 * t241 + t275 * t297 + t289;
t210 = t275 * t308 + t278 * t292;
t285 = qJD(3) * t271 + t248 * t210 + t288;
t284 = t267 * t243 + V_base(2) + (-pkin(6) - t309) * V_base(4);
t180 = pkin(3) * t213 - pkin(8) * t212;
t207 = pkin(3) * t230 - pkin(8) * t229;
t282 = t247 * t207 + (-t180 - t210) * t249 + t286;
t181 = pkin(3) * t215 - pkin(8) * t214;
t211 = -t275 * t292 + t278 * t308;
t281 = t248 * t180 + (-t181 - t211) * t247 + t285;
t280 = t249 * t211 - t278 * t297 + t284;
t279 = t249 * t181 + (-t207 - t241) * t248 + t280;
t276 = cos(qJ(5));
t273 = sin(qJ(4));
t272 = sin(qJ(5));
t268 = Icges(2,4) * t278;
t257 = rSges(2,1) * t278 - t275 * rSges(2,2);
t256 = t275 * rSges(2,1) + rSges(2,2) * t278;
t255 = Icges(2,1) * t278 - t305;
t254 = Icges(2,1) * t275 + t268;
t253 = -Icges(2,2) * t275 + t268;
t252 = Icges(2,2) * t278 + t305;
t246 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t245 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t244 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t228 = rSges(3,3) * t271 + (rSges(3,1) * t274 + rSges(3,2) * t277) * t270;
t227 = Icges(3,5) * t271 + (Icges(3,1) * t274 + Icges(3,4) * t277) * t270;
t226 = Icges(3,6) * t271 + (Icges(3,4) * t274 + Icges(3,2) * t277) * t270;
t221 = V_base(5) * rSges(2,3) - t256 * t267 + t296;
t220 = t257 * t267 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t219 = t256 * V_base(4) - t257 * V_base(5) + V_base(3);
t218 = t230 * t311 + t271 * t273;
t217 = t230 * t273 - t271 * t311;
t206 = rSges(3,1) * t237 + rSges(3,2) * t236 + rSges(3,3) * t304;
t205 = t235 * rSges(3,1) + t234 * rSges(3,2) - rSges(3,3) * t303;
t204 = Icges(3,1) * t237 + Icges(3,4) * t236 + Icges(3,5) * t304;
t203 = Icges(3,1) * t235 + Icges(3,4) * t234 - Icges(3,5) * t303;
t202 = Icges(3,4) * t237 + Icges(3,2) * t236 + Icges(3,6) * t304;
t201 = Icges(3,4) * t235 + Icges(3,2) * t234 - Icges(3,6) * t303;
t198 = rSges(4,1) * t230 + rSges(4,2) * t229 + rSges(4,3) * t271;
t197 = Icges(4,1) * t230 + Icges(4,4) * t229 + Icges(4,5) * t271;
t196 = Icges(4,4) * t230 + Icges(4,2) * t229 + Icges(4,6) * t271;
t193 = t215 * t311 + t273 * t304;
t192 = t215 * t273 - t275 * t293;
t191 = t213 * t311 - t273 * t303;
t190 = t213 * t273 + t278 * t293;
t185 = t218 * t276 - t229 * t272;
t184 = -t218 * t272 - t229 * t276;
t183 = qJD(5) * t217 + t216;
t182 = pkin(4) * t218 + pkin(9) * t217;
t179 = rSges(5,1) * t218 - rSges(5,2) * t217 - rSges(5,3) * t229;
t178 = Icges(5,1) * t218 - Icges(5,4) * t217 - Icges(5,5) * t229;
t177 = Icges(5,4) * t218 - Icges(5,2) * t217 - Icges(5,6) * t229;
t176 = Icges(5,5) * t218 - Icges(5,6) * t217 - Icges(5,3) * t229;
t175 = rSges(4,1) * t215 + rSges(4,2) * t214 + rSges(4,3) * t304;
t174 = t213 * rSges(4,1) + t212 * rSges(4,2) - rSges(4,3) * t303;
t173 = Icges(4,1) * t215 + Icges(4,4) * t214 + Icges(4,5) * t304;
t172 = Icges(4,1) * t213 + Icges(4,4) * t212 - Icges(4,5) * t303;
t171 = Icges(4,4) * t215 + Icges(4,2) * t214 + Icges(4,6) * t304;
t170 = Icges(4,4) * t213 + Icges(4,2) * t212 - Icges(4,6) * t303;
t165 = t193 * t276 - t214 * t272;
t164 = -t193 * t272 - t214 * t276;
t163 = t191 * t276 - t212 * t272;
t162 = -t191 * t272 - t212 * t276;
t161 = qJD(5) * t192 + t189;
t160 = qJD(5) * t190 + t188;
t159 = pkin(4) * t193 + pkin(9) * t192;
t158 = pkin(4) * t191 + pkin(9) * t190;
t157 = -t205 * t249 + t228 * t247 + t289;
t156 = t206 * t249 - t228 * t248 + t284;
t155 = t205 * t248 - t206 * t247 + t288;
t154 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t217;
t153 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t217;
t152 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t217;
t151 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t217;
t150 = rSges(5,1) * t193 - rSges(5,2) * t192 - rSges(5,3) * t214;
t149 = rSges(5,1) * t191 - rSges(5,2) * t190 - rSges(5,3) * t212;
t148 = Icges(5,1) * t193 - Icges(5,4) * t192 - Icges(5,5) * t214;
t147 = Icges(5,1) * t191 - Icges(5,4) * t190 - Icges(5,5) * t212;
t146 = Icges(5,4) * t193 - Icges(5,2) * t192 - Icges(5,6) * t214;
t145 = Icges(5,4) * t191 - Icges(5,2) * t190 - Icges(5,6) * t212;
t144 = Icges(5,5) * t193 - Icges(5,6) * t192 - Icges(5,3) * t214;
t143 = Icges(5,5) * t191 - Icges(5,6) * t190 - Icges(5,3) * t212;
t142 = rSges(6,1) * t165 + rSges(6,2) * t164 + rSges(6,3) * t192;
t141 = rSges(6,1) * t163 + rSges(6,2) * t162 + rSges(6,3) * t190;
t140 = Icges(6,1) * t165 + Icges(6,4) * t164 + Icges(6,5) * t192;
t139 = Icges(6,1) * t163 + Icges(6,4) * t162 + Icges(6,5) * t190;
t138 = Icges(6,4) * t165 + Icges(6,2) * t164 + Icges(6,6) * t192;
t137 = Icges(6,4) * t163 + Icges(6,2) * t162 + Icges(6,6) * t190;
t136 = Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t192;
t135 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t190;
t134 = t198 * t247 + (-t174 - t210) * t249 + t286;
t133 = t249 * t175 + (-t198 - t241) * t248 + t280;
t132 = t174 * t248 + (-t175 - t211) * t247 + t285;
t131 = -t149 * t216 + t179 * t188 + t282;
t130 = t216 * t150 - t189 * t179 + t279;
t129 = t149 * t189 - t150 * t188 + t281;
t128 = -t141 * t183 + t154 * t160 - t158 * t216 + t182 * t188 + t282;
t127 = t183 * t142 - t161 * t154 + t216 * t159 - t189 * t182 + t279;
t126 = t141 * t161 - t142 * t160 + t158 * t189 - t159 * t188 + t281;
t1 = m(1) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + m(2) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(3) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(4) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(5) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + t189 * ((-t214 * t144 - t192 * t146 + t193 * t148) * t189 + (-t143 * t214 - t145 * t192 + t147 * t193) * t188 + (-t176 * t214 - t177 * t192 + t178 * t193) * t216) / 0.2e1 + t188 * ((-t144 * t212 - t146 * t190 + t148 * t191) * t189 + (-t212 * t143 - t190 * t145 + t191 * t147) * t188 + (-t176 * t212 - t177 * t190 + t178 * t191) * t216) / 0.2e1 + t216 * ((-t144 * t229 - t146 * t217 + t148 * t218) * t189 + (-t143 * t229 - t145 * t217 + t147 * t218) * t188 + (-t229 * t176 - t217 * t177 + t218 * t178) * t216) / 0.2e1 + m(6) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + t161 * ((t192 * t136 + t164 * t138 + t165 * t140) * t161 + (t135 * t192 + t137 * t164 + t139 * t165) * t160 + (t151 * t192 + t152 * t164 + t153 * t165) * t183) / 0.2e1 + t160 * ((t136 * t190 + t138 * t162 + t140 * t163) * t161 + (t190 * t135 + t162 * t137 + t163 * t139) * t160 + (t151 * t190 + t152 * t162 + t153 * t163) * t183) / 0.2e1 + t183 * ((t136 * t217 + t138 * t184 + t140 * t185) * t161 + (t135 * t217 + t137 * t184 + t139 * t185) * t160 + (t217 * t151 + t184 * t152 + t185 * t153) * t183) / 0.2e1 + ((t212 * t196 + t213 * t197 + t234 * t226 + t235 * t227 - t303 * t315) * t249 + (t212 * t171 + t213 * t173 + t234 * t202 + t235 * t204 - t303 * t316) * t248 + (t212 * t170 + t213 * t172 + t234 * t201 + t235 * t203 - t303 * t317) * t247) * t247 / 0.2e1 + ((t196 * t214 + t197 * t215 + t226 * t236 + t227 * t237 + t304 * t315) * t249 + (t171 * t214 + t173 * t215 + t202 * t236 + t204 * t237 + t304 * t316) * t248 + (t170 * t214 + t172 * t215 + t201 * t236 + t203 * t237 + t304 * t317) * t247) * t248 / 0.2e1 + ((t199 * t247 + t200 * t248 + t225 * t249) * t271 + ((t202 * t277 + t204 * t274) * t248 + (t201 * t277 + t203 * t274) * t247 + (t226 * t277 + t227 * t274) * t249) * t270 + (t169 * t271 + t171 * t229 + t173 * t230) * t248 + (t168 * t271 + t170 * t229 + t172 * t230) * t247 + (t195 * t271 + t196 * t229 + t197 * t230) * t249) * t249 / 0.2e1 + ((-t275 * t252 + t254 * t278 + Icges(1,4)) * V_base(5) + (-t275 * t253 + t255 * t278 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t252 * t278 + t275 * t254 + Icges(1,2)) * V_base(5) + (t253 * t278 + t275 * t255 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t275 + Icges(2,6) * t278) * V_base(5) + (Icges(2,5) * t278 - Icges(2,6) * t275) * V_base(4) + Icges(2,3) * t267 / 0.2e1) * t267;
T = t1;
