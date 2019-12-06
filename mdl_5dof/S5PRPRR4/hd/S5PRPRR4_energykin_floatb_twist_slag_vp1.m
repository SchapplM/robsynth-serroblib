% Calculate kinetic energy for
% S5PRPRR4
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:32
% EndTime: 2019-12-05 15:49:36
% DurationCPUTime: 3.55s
% Computational Cost: add. (2411->357), mult. (5755->519), div. (0->0), fcn. (7113->12), ass. (0->158)
t266 = sin(pkin(9));
t268 = cos(pkin(9));
t272 = sin(qJ(2));
t274 = cos(qJ(2));
t303 = sin(pkin(10));
t304 = cos(pkin(10));
t237 = -t272 * t303 + t274 * t304;
t269 = cos(pkin(5));
t279 = t269 * t237;
t283 = t272 * t304 + t274 * t303;
t209 = -t266 * t283 + t268 * t279;
t228 = t283 * t269;
t210 = t228 * t268 + t237 * t266;
t267 = sin(pkin(5));
t300 = t267 * t268;
t165 = Icges(4,5) * t210 + Icges(4,6) * t209 - Icges(4,3) * t300;
t297 = t269 * t274;
t230 = -t266 * t272 + t268 * t297;
t298 = t269 * t272;
t231 = t266 * t274 + t268 * t298;
t196 = Icges(3,5) * t231 + Icges(3,6) * t230 - Icges(3,3) * t300;
t313 = t165 + t196;
t211 = -t266 * t279 - t268 * t283;
t212 = -t228 * t266 + t237 * t268;
t301 = t266 * t267;
t166 = Icges(4,5) * t212 + Icges(4,6) * t211 + Icges(4,3) * t301;
t232 = -t266 * t297 - t268 * t272;
t233 = -t266 * t298 + t268 * t274;
t197 = Icges(3,5) * t233 + Icges(3,6) * t232 + Icges(3,3) * t301;
t312 = t166 + t197;
t226 = t237 * t267;
t227 = t283 * t267;
t192 = Icges(4,5) * t227 + Icges(4,6) * t226 + Icges(4,3) * t269;
t222 = Icges(3,3) * t269 + (Icges(3,5) * t272 + Icges(3,6) * t274) * t267;
t311 = t192 + t222;
t307 = cos(qJ(4));
t306 = pkin(6) * t269;
t305 = pkin(2) * t274;
t302 = Icges(2,4) * t266;
t271 = sin(qJ(4));
t299 = t267 * t271;
t296 = qJD(2) * t267;
t295 = qJD(3) * t267;
t294 = V_base(5) * qJ(1) + V_base(1);
t290 = qJD(1) + V_base(3);
t289 = t267 * t307;
t245 = t266 * t296 + V_base(4);
t256 = qJD(2) * t269 + V_base(6);
t288 = pkin(2) * t298 - qJ(3) * t267;
t186 = -qJD(4) * t211 + t245;
t214 = -qJD(4) * t226 + t256;
t244 = -t268 * t296 + V_base(5);
t185 = -qJD(4) * t209 + t244;
t239 = pkin(1) * t266 - pkin(6) * t300;
t285 = -t239 * V_base(6) + V_base(5) * t306 + t294;
t240 = pkin(1) * t268 + pkin(6) * t301;
t284 = V_base(4) * t239 - t240 * V_base(5) + t290;
t282 = V_base(6) * t240 + V_base(2) + (-qJ(1) - t306) * V_base(4);
t238 = pkin(2) * t267 * t272 + qJ(3) * t269;
t281 = t244 * t238 + t266 * t295 + t285;
t204 = t266 * t305 + t268 * t288;
t280 = qJD(3) * t269 + t245 * t204 + t284;
t205 = -t266 * t288 + t268 * t305;
t278 = t256 * t205 - t268 * t295 + t282;
t173 = pkin(3) * t210 - pkin(7) * t209;
t208 = pkin(3) * t227 - pkin(7) * t226;
t277 = t244 * t208 + (-t173 - t204) * t256 + t281;
t174 = pkin(3) * t212 - pkin(7) * t211;
t276 = t245 * t173 + (-t174 - t205) * t244 + t280;
t275 = t256 * t174 + (-t208 - t238) * t245 + t278;
t273 = cos(qJ(5));
t270 = sin(qJ(5));
t264 = Icges(2,4) * t268;
t253 = rSges(2,1) * t268 - rSges(2,2) * t266;
t252 = rSges(2,1) * t266 + rSges(2,2) * t268;
t251 = Icges(2,1) * t268 - t302;
t250 = Icges(2,1) * t266 + t264;
t249 = -Icges(2,2) * t266 + t264;
t248 = Icges(2,2) * t268 + t302;
t243 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t242 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t241 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t225 = t269 * rSges(3,3) + (rSges(3,1) * t272 + rSges(3,2) * t274) * t267;
t224 = Icges(3,5) * t269 + (Icges(3,1) * t272 + Icges(3,4) * t274) * t267;
t223 = Icges(3,6) * t269 + (Icges(3,4) * t272 + Icges(3,2) * t274) * t267;
t218 = V_base(5) * rSges(2,3) - t252 * V_base(6) + t294;
t217 = t253 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t216 = t227 * t307 + t269 * t271;
t215 = t227 * t271 - t269 * t307;
t213 = t252 * V_base(4) - t253 * V_base(5) + t290;
t203 = rSges(3,1) * t233 + rSges(3,2) * t232 + rSges(3,3) * t301;
t202 = rSges(3,1) * t231 + rSges(3,2) * t230 - rSges(3,3) * t300;
t201 = Icges(3,1) * t233 + Icges(3,4) * t232 + Icges(3,5) * t301;
t200 = Icges(3,1) * t231 + Icges(3,4) * t230 - Icges(3,5) * t300;
t199 = Icges(3,4) * t233 + Icges(3,2) * t232 + Icges(3,6) * t301;
t198 = Icges(3,4) * t231 + Icges(3,2) * t230 - Icges(3,6) * t300;
t195 = rSges(4,1) * t227 + rSges(4,2) * t226 + rSges(4,3) * t269;
t194 = Icges(4,1) * t227 + Icges(4,4) * t226 + Icges(4,5) * t269;
t193 = Icges(4,4) * t227 + Icges(4,2) * t226 + Icges(4,6) * t269;
t190 = t212 * t307 + t266 * t299;
t189 = t212 * t271 - t266 * t289;
t188 = t210 * t307 - t268 * t299;
t187 = t210 * t271 + t268 * t289;
t182 = t216 * t273 - t226 * t270;
t181 = -t216 * t270 - t226 * t273;
t180 = qJD(5) * t215 + t214;
t179 = pkin(4) * t216 + pkin(8) * t215;
t178 = rSges(5,1) * t216 - rSges(5,2) * t215 - rSges(5,3) * t226;
t177 = Icges(5,1) * t216 - Icges(5,4) * t215 - Icges(5,5) * t226;
t176 = Icges(5,4) * t216 - Icges(5,2) * t215 - Icges(5,6) * t226;
t175 = Icges(5,5) * t216 - Icges(5,6) * t215 - Icges(5,3) * t226;
t172 = rSges(4,1) * t212 + rSges(4,2) * t211 + rSges(4,3) * t301;
t171 = rSges(4,1) * t210 + rSges(4,2) * t209 - rSges(4,3) * t300;
t170 = Icges(4,1) * t212 + Icges(4,4) * t211 + Icges(4,5) * t301;
t169 = Icges(4,1) * t210 + Icges(4,4) * t209 - Icges(4,5) * t300;
t168 = Icges(4,4) * t212 + Icges(4,2) * t211 + Icges(4,6) * t301;
t167 = Icges(4,4) * t210 + Icges(4,2) * t209 - Icges(4,6) * t300;
t162 = t190 * t273 - t211 * t270;
t161 = -t190 * t270 - t211 * t273;
t160 = t188 * t273 - t209 * t270;
t159 = -t188 * t270 - t209 * t273;
t158 = qJD(5) * t189 + t186;
t157 = qJD(5) * t187 + t185;
t156 = pkin(4) * t190 + pkin(8) * t189;
t155 = pkin(4) * t188 + pkin(8) * t187;
t154 = -t202 * t256 + t225 * t244 + t285;
t153 = t203 * t256 - t225 * t245 + t282;
t152 = rSges(6,1) * t182 + rSges(6,2) * t181 + rSges(6,3) * t215;
t151 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t215;
t150 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t215;
t149 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t215;
t148 = t202 * t245 - t203 * t244 + t284;
t147 = rSges(5,1) * t190 - rSges(5,2) * t189 - rSges(5,3) * t211;
t146 = rSges(5,1) * t188 - rSges(5,2) * t187 - rSges(5,3) * t209;
t145 = Icges(5,1) * t190 - Icges(5,4) * t189 - Icges(5,5) * t211;
t144 = Icges(5,1) * t188 - Icges(5,4) * t187 - Icges(5,5) * t209;
t143 = Icges(5,4) * t190 - Icges(5,2) * t189 - Icges(5,6) * t211;
t142 = Icges(5,4) * t188 - Icges(5,2) * t187 - Icges(5,6) * t209;
t141 = Icges(5,5) * t190 - Icges(5,6) * t189 - Icges(5,3) * t211;
t140 = Icges(5,5) * t188 - Icges(5,6) * t187 - Icges(5,3) * t209;
t139 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t189;
t138 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t187;
t137 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t189;
t136 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t187;
t135 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t189;
t134 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t187;
t133 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t189;
t132 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t187;
t131 = t195 * t244 + (-t171 - t204) * t256 + t281;
t130 = t172 * t256 + (-t195 - t238) * t245 + t278;
t129 = t171 * t245 + (-t172 - t205) * t244 + t280;
t128 = -t146 * t214 + t178 * t185 + t277;
t127 = t147 * t214 - t178 * t186 + t275;
t126 = t146 * t186 - t147 * t185 + t276;
t125 = -t138 * t180 + t152 * t157 - t155 * t214 + t179 * t185 + t277;
t124 = t139 * t180 - t152 * t158 + t156 * t214 - t179 * t186 + t275;
t123 = t138 * t158 - t139 * t157 + t155 * t186 - t156 * t185 + t276;
t1 = m(1) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + m(2) * (t213 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + m(3) * (t148 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(4) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(5) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + t186 * ((-t141 * t211 - t143 * t189 + t145 * t190) * t186 + (-t140 * t211 - t142 * t189 + t144 * t190) * t185 + (-t175 * t211 - t176 * t189 + t177 * t190) * t214) / 0.2e1 + t185 * ((-t141 * t209 - t143 * t187 + t145 * t188) * t186 + (-t140 * t209 - t142 * t187 + t144 * t188) * t185 + (-t175 * t209 - t176 * t187 + t177 * t188) * t214) / 0.2e1 + t214 * ((-t141 * t226 - t143 * t215 + t145 * t216) * t186 + (-t140 * t226 - t142 * t215 + t144 * t216) * t185 + (-t175 * t226 - t176 * t215 + t177 * t216) * t214) / 0.2e1 + m(6) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + t158 * ((t133 * t189 + t135 * t161 + t137 * t162) * t158 + (t132 * t189 + t134 * t161 + t136 * t162) * t157 + (t149 * t189 + t150 * t161 + t151 * t162) * t180) / 0.2e1 + t157 * ((t133 * t187 + t135 * t159 + t137 * t160) * t158 + (t132 * t187 + t134 * t159 + t136 * t160) * t157 + (t149 * t187 + t150 * t159 + t151 * t160) * t180) / 0.2e1 + t180 * ((t133 * t215 + t135 * t181 + t137 * t182) * t158 + (t132 * t215 + t134 * t181 + t136 * t182) * t157 + (t149 * t215 + t150 * t181 + t151 * t182) * t180) / 0.2e1 + ((t193 * t209 + t194 * t210 + t223 * t230 + t224 * t231 - t311 * t300) * t256 + (t168 * t209 + t170 * t210 + t199 * t230 + t201 * t231 - t312 * t300) * t245 + (t167 * t209 + t169 * t210 + t198 * t230 + t200 * t231 - t313 * t300) * t244) * t244 / 0.2e1 + ((t193 * t211 + t194 * t212 + t223 * t232 + t224 * t233 + t311 * t301) * t256 + (t168 * t211 + t170 * t212 + t199 * t232 + t201 * t233 + t312 * t301) * t245 + (t167 * t211 + t169 * t212 + t198 * t232 + t200 * t233 + t313 * t301) * t244) * t245 / 0.2e1 + ((t196 * t244 + t197 * t245 + t222 * t256) * t269 + ((t199 * t274 + t201 * t272) * t245 + (t198 * t274 + t200 * t272) * t244 + (t223 * t274 + t224 * t272) * t256) * t267 + (t166 * t269 + t168 * t226 + t170 * t227) * t245 + (t165 * t269 + t167 * t226 + t169 * t227) * t244 + (t192 * t269 + t193 * t226 + t194 * t227) * t256) * t256 / 0.2e1 + ((-t248 * t266 + t250 * t268 + Icges(1,4)) * V_base(5) + (-t249 * t266 + t251 * t268 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t248 * t268 + t250 * t266 + Icges(1,2)) * V_base(5) + (t249 * t268 + t251 * t266 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t266 + Icges(2,6) * t268 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t268 - Icges(2,6) * t266 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
