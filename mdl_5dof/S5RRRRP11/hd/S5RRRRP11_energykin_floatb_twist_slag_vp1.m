% Calculate kinetic energy for
% S5RRRRP11
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:21
% EndTime: 2019-12-31 22:14:24
% DurationCPUTime: 3.16s
% Computational Cost: add. (2016->304), mult. (4647->442), div. (0->0), fcn. (5577->10), ass. (0->138)
t315 = Icges(5,1) + Icges(6,1);
t314 = -Icges(5,4) + Icges(6,5);
t313 = Icges(6,4) + Icges(5,5);
t312 = Icges(5,2) + Icges(6,3);
t311 = Icges(6,2) + Icges(5,3);
t310 = -Icges(5,6) + Icges(6,6);
t309 = rSges(6,1) + pkin(4);
t308 = rSges(6,3) + qJ(5);
t261 = cos(pkin(5));
t265 = sin(qJ(1));
t266 = cos(qJ(2));
t286 = t265 * t266;
t264 = sin(qJ(2));
t267 = cos(qJ(1));
t287 = t264 * t267;
t230 = t261 * t287 + t286;
t263 = sin(qJ(3));
t260 = sin(pkin(5));
t289 = t260 * t267;
t295 = cos(qJ(3));
t210 = t230 * t295 - t263 * t289;
t285 = t266 * t267;
t288 = t264 * t265;
t229 = -t261 * t285 + t288;
t262 = sin(qJ(4));
t294 = cos(qJ(4));
t183 = t210 * t262 - t229 * t294;
t184 = t210 * t294 + t229 * t262;
t277 = t260 * t295;
t209 = t230 * t263 + t267 * t277;
t307 = t183 * t312 + t184 * t314 + t209 * t310;
t232 = -t261 * t288 + t285;
t291 = t260 * t265;
t212 = t232 * t295 + t263 * t291;
t231 = t261 * t286 + t287;
t185 = t212 * t262 - t231 * t294;
t186 = t212 * t294 + t231 * t262;
t211 = t232 * t263 - t265 * t277;
t306 = t185 * t312 + t186 * t314 + t211 * t310;
t305 = t183 * t310 + t184 * t313 + t209 * t311;
t304 = t185 * t310 + t186 * t313 + t211 * t311;
t303 = t183 * t314 + t184 * t315 + t209 * t313;
t302 = t185 * t314 + t186 * t315 + t211 * t313;
t228 = t261 * t263 + t264 * t277;
t290 = t260 * t266;
t205 = t228 * t262 + t290 * t294;
t206 = t228 * t294 - t262 * t290;
t227 = t260 * t263 * t264 - t261 * t295;
t301 = t205 * t312 + t206 * t314 + t227 * t310;
t300 = t205 * t310 + t206 * t313 + t227 * t311;
t299 = t205 * t314 + t206 * t315 + t227 * t313;
t293 = pkin(7) * t261;
t292 = Icges(2,4) * t265;
t284 = rSges(6,2) * t209 + t183 * t308 + t184 * t309;
t283 = rSges(6,2) * t211 + t185 * t308 + t186 * t309;
t282 = rSges(6,2) * t227 + t205 * t308 + t206 * t309;
t281 = qJD(2) * t260;
t280 = V_base(5) * pkin(6) + V_base(1);
t241 = t265 * t281 + V_base(4);
t257 = V_base(6) + qJD(1);
t208 = qJD(3) * t231 + t241;
t242 = qJD(2) * t261 + t257;
t240 = -t267 * t281 + V_base(5);
t235 = pkin(1) * t265 - pkin(7) * t289;
t276 = -t235 * t257 + V_base(5) * t293 + t280;
t236 = pkin(1) * t267 + pkin(7) * t291;
t275 = V_base(4) * t235 - t236 * V_base(5) + V_base(3);
t207 = qJD(3) * t229 + t240;
t225 = -qJD(3) * t290 + t242;
t274 = t257 * t236 + V_base(2) + (-pkin(6) - t293) * V_base(4);
t203 = pkin(2) * t230 + pkin(8) * t229;
t234 = (pkin(2) * t264 - pkin(8) * t266) * t260;
t273 = -t203 * t242 + t240 * t234 + t276;
t204 = pkin(2) * t232 + pkin(8) * t231;
t272 = t203 * t241 - t204 * t240 + t275;
t271 = t204 * t242 - t234 * t241 + t274;
t179 = pkin(3) * t210 + pkin(9) * t209;
t201 = pkin(3) * t228 + pkin(9) * t227;
t270 = -t179 * t225 + t201 * t207 + t273;
t180 = pkin(3) * t212 + pkin(9) * t211;
t269 = t179 * t208 - t180 * t207 + t272;
t268 = t180 * t225 - t201 * t208 + t271;
t258 = Icges(2,4) * t267;
t250 = rSges(2,1) * t267 - rSges(2,2) * t265;
t249 = rSges(2,1) * t265 + rSges(2,2) * t267;
t248 = Icges(2,1) * t267 - t292;
t247 = Icges(2,1) * t265 + t258;
t246 = -Icges(2,2) * t265 + t258;
t245 = Icges(2,2) * t267 + t292;
t239 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t238 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t237 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t220 = rSges(3,3) * t261 + (rSges(3,1) * t264 + rSges(3,2) * t266) * t260;
t219 = Icges(3,5) * t261 + (Icges(3,1) * t264 + Icges(3,4) * t266) * t260;
t218 = Icges(3,6) * t261 + (Icges(3,4) * t264 + Icges(3,2) * t266) * t260;
t217 = Icges(3,3) * t261 + (Icges(3,5) * t264 + Icges(3,6) * t266) * t260;
t216 = V_base(5) * rSges(2,3) - t249 * t257 + t280;
t215 = t250 * t257 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t213 = t249 * V_base(4) - t250 * V_base(5) + V_base(3);
t202 = qJD(4) * t227 + t225;
t200 = rSges(3,1) * t232 - rSges(3,2) * t231 + rSges(3,3) * t291;
t199 = rSges(3,1) * t230 - rSges(3,2) * t229 - rSges(3,3) * t289;
t198 = Icges(3,1) * t232 - Icges(3,4) * t231 + Icges(3,5) * t291;
t197 = Icges(3,1) * t230 - Icges(3,4) * t229 - Icges(3,5) * t289;
t196 = Icges(3,4) * t232 - Icges(3,2) * t231 + Icges(3,6) * t291;
t195 = Icges(3,4) * t230 - Icges(3,2) * t229 - Icges(3,6) * t289;
t194 = Icges(3,5) * t232 - Icges(3,6) * t231 + Icges(3,3) * t291;
t193 = Icges(3,5) * t230 - Icges(3,6) * t229 - Icges(3,3) * t289;
t192 = rSges(4,1) * t228 - rSges(4,2) * t227 - rSges(4,3) * t290;
t191 = Icges(4,1) * t228 - Icges(4,4) * t227 - Icges(4,5) * t290;
t190 = Icges(4,4) * t228 - Icges(4,2) * t227 - Icges(4,6) * t290;
t189 = Icges(4,5) * t228 - Icges(4,6) * t227 - Icges(4,3) * t290;
t182 = qJD(4) * t211 + t208;
t181 = qJD(4) * t209 + t207;
t176 = rSges(4,1) * t212 - rSges(4,2) * t211 + rSges(4,3) * t231;
t175 = rSges(4,1) * t210 - rSges(4,2) * t209 + rSges(4,3) * t229;
t173 = Icges(4,1) * t212 - Icges(4,4) * t211 + Icges(4,5) * t231;
t172 = Icges(4,1) * t210 - Icges(4,4) * t209 + Icges(4,5) * t229;
t171 = Icges(4,4) * t212 - Icges(4,2) * t211 + Icges(4,6) * t231;
t170 = Icges(4,4) * t210 - Icges(4,2) * t209 + Icges(4,6) * t229;
t169 = Icges(4,5) * t212 - Icges(4,6) * t211 + Icges(4,3) * t231;
t168 = Icges(4,5) * t210 - Icges(4,6) * t209 + Icges(4,3) * t229;
t167 = rSges(5,1) * t206 - rSges(5,2) * t205 + rSges(5,3) * t227;
t156 = -t199 * t242 + t220 * t240 + t276;
t155 = t200 * t242 - t220 * t241 + t274;
t154 = rSges(5,1) * t186 - rSges(5,2) * t185 + rSges(5,3) * t211;
t152 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t209;
t138 = t199 * t241 - t200 * t240 + t275;
t137 = -t175 * t225 + t192 * t207 + t273;
t136 = t176 * t225 - t192 * t208 + t271;
t135 = t175 * t208 - t176 * t207 + t272;
t134 = -t152 * t202 + t167 * t181 + t270;
t133 = t154 * t202 - t167 * t182 + t268;
t132 = t152 * t182 - t154 * t181 + t269;
t131 = qJD(5) * t185 + t181 * t282 - t202 * t284 + t270;
t130 = qJD(5) * t183 - t182 * t282 + t202 * t283 + t268;
t129 = qJD(5) * t205 - t181 * t283 + t182 * t284 + t269;
t1 = m(1) * (t237 ^ 2 + t238 ^ 2 + t239 ^ 2) / 0.2e1 + m(2) * (t213 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 + m(3) * (t138 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + t241 * ((t194 * t291 - t196 * t231 + t198 * t232) * t241 + (t193 * t291 - t195 * t231 + t197 * t232) * t240 + (t217 * t291 - t218 * t231 + t219 * t232) * t242) / 0.2e1 + t240 * ((-t194 * t289 - t196 * t229 + t198 * t230) * t241 + (-t193 * t289 - t229 * t195 + t230 * t197) * t240 + (-t217 * t289 - t218 * t229 + t219 * t230) * t242) / 0.2e1 + t242 * ((t193 * t240 + t194 * t241 + t217 * t242) * t261 + ((t196 * t266 + t198 * t264) * t241 + (t195 * t266 + t197 * t264) * t240 + (t218 * t266 + t219 * t264) * t242) * t260) / 0.2e1 + m(4) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + t208 * ((t169 * t231 - t171 * t211 + t173 * t212) * t208 + (t168 * t231 - t170 * t211 + t172 * t212) * t207 + (t189 * t231 - t190 * t211 + t191 * t212) * t225) / 0.2e1 + t207 * ((t169 * t229 - t171 * t209 + t173 * t210) * t208 + (t168 * t229 - t170 * t209 + t172 * t210) * t207 + (t189 * t229 - t190 * t209 + t191 * t210) * t225) / 0.2e1 + t225 * ((-t169 * t290 - t171 * t227 + t173 * t228) * t208 + (-t168 * t290 - t170 * t227 + t172 * t228) * t207 + (-t189 * t290 - t190 * t227 + t191 * t228) * t225) / 0.2e1 + m(5) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(6) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + ((t183 * t301 + t184 * t299 + t209 * t300) * t202 + (t183 * t306 + t184 * t302 + t209 * t304) * t182 + (t307 * t183 + t303 * t184 + t305 * t209) * t181) * t181 / 0.2e1 + ((t185 * t301 + t186 * t299 + t211 * t300) * t202 + (t306 * t185 + t302 * t186 + t304 * t211) * t182 + (t185 * t307 + t186 * t303 + t211 * t305) * t181) * t182 / 0.2e1 + ((t301 * t205 + t299 * t206 + t300 * t227) * t202 + (t205 * t306 + t206 * t302 + t227 * t304) * t182 + (t205 * t307 + t206 * t303 + t227 * t305) * t181) * t202 / 0.2e1 + ((-t245 * t265 + t247 * t267 + Icges(1,4)) * V_base(5) + (-t265 * t246 + t248 * t267 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t245 * t267 + t265 * t247 + Icges(1,2)) * V_base(5) + (t246 * t267 + t248 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t265 + Icges(2,6) * t267) * V_base(5) + (Icges(2,5) * t267 - Icges(2,6) * t265) * V_base(4) + Icges(2,3) * t257 / 0.2e1) * t257;
T = t1;
