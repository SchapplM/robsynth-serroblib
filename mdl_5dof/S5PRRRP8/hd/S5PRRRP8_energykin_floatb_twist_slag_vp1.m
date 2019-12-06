% Calculate kinetic energy for
% S5PRRRP8
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:08
% EndTime: 2019-12-05 16:58:12
% DurationCPUTime: 3.40s
% Computational Cost: add. (1971->304), mult. (4647->440), div. (0->0), fcn. (5577->10), ass. (0->137)
t315 = Icges(5,1) + Icges(6,1);
t314 = -Icges(5,4) + Icges(6,5);
t313 = Icges(6,4) + Icges(5,5);
t312 = Icges(5,2) + Icges(6,3);
t311 = Icges(6,2) + Icges(5,3);
t310 = -Icges(5,6) + Icges(6,6);
t309 = rSges(6,1) + pkin(4);
t308 = rSges(6,3) + qJ(5);
t257 = sin(pkin(9));
t259 = cos(pkin(9));
t264 = cos(qJ(2));
t260 = cos(pkin(5));
t263 = sin(qJ(2));
t285 = t260 * t263;
t224 = t257 * t264 + t259 * t285;
t258 = sin(pkin(5));
t262 = sin(qJ(3));
t287 = t258 * t262;
t293 = cos(qJ(3));
t207 = t224 * t293 - t259 * t287;
t284 = t260 * t264;
t223 = t257 * t263 - t259 * t284;
t261 = sin(qJ(4));
t292 = cos(qJ(4));
t181 = t207 * t261 - t223 * t292;
t182 = t207 * t292 + t223 * t261;
t274 = t258 * t293;
t206 = t224 * t262 + t259 * t274;
t305 = t312 * t181 + t314 * t182 + t310 * t206;
t226 = -t257 * t285 + t259 * t264;
t209 = t226 * t293 + t257 * t287;
t225 = t257 * t284 + t259 * t263;
t183 = t209 * t261 - t225 * t292;
t184 = t209 * t292 + t225 * t261;
t208 = t226 * t262 - t257 * t274;
t304 = t312 * t183 + t314 * t184 + t310 * t208;
t303 = t310 * t181 + t313 * t182 + t311 * t206;
t302 = t310 * t183 + t313 * t184 + t311 * t208;
t301 = t314 * t181 + t315 * t182 + t313 * t206;
t300 = t314 * t183 + t315 * t184 + t313 * t208;
t231 = t260 * t262 + t263 * t274;
t286 = t258 * t264;
t210 = t231 * t261 + t286 * t292;
t211 = t231 * t292 - t261 * t286;
t230 = -t260 * t293 + t263 * t287;
t299 = t312 * t210 + t314 * t211 + t310 * t230;
t298 = t310 * t210 + t313 * t211 + t311 * t230;
t297 = t314 * t210 + t315 * t211 + t313 * t230;
t291 = pkin(6) * t260;
t290 = Icges(2,4) * t257;
t289 = t257 * t258;
t288 = t258 * t259;
t283 = rSges(6,2) * t206 + t308 * t181 + t309 * t182;
t282 = rSges(6,2) * t208 + t308 * t183 + t309 * t184;
t281 = rSges(6,2) * t230 + t308 * t210 + t309 * t211;
t280 = qJD(2) * t258;
t279 = V_base(5) * qJ(1) + V_base(1);
t275 = qJD(1) + V_base(3);
t239 = t257 * t280 + V_base(4);
t250 = qJD(2) * t260 + V_base(6);
t205 = qJD(3) * t225 + t239;
t238 = -t259 * t280 + V_base(5);
t204 = qJD(3) * t223 + t238;
t227 = -qJD(3) * t286 + t250;
t233 = pkin(1) * t257 - pkin(6) * t288;
t273 = -t233 * V_base(6) + V_base(5) * t291 + t279;
t234 = pkin(1) * t259 + pkin(6) * t289;
t272 = V_base(4) * t233 - t234 * V_base(5) + t275;
t271 = V_base(6) * t234 + V_base(2) + (-qJ(1) - t291) * V_base(4);
t199 = pkin(2) * t224 + pkin(7) * t223;
t232 = (pkin(2) * t263 - pkin(7) * t264) * t258;
t270 = -t199 * t250 + t238 * t232 + t273;
t200 = pkin(2) * t226 + pkin(7) * t225;
t269 = t239 * t199 - t200 * t238 + t272;
t268 = t250 * t200 - t232 * t239 + t271;
t176 = pkin(3) * t207 + pkin(8) * t206;
t201 = pkin(3) * t231 + pkin(8) * t230;
t267 = -t176 * t227 + t204 * t201 + t270;
t177 = pkin(3) * t209 + pkin(8) * t208;
t266 = t205 * t176 - t177 * t204 + t269;
t265 = t227 * t177 - t201 * t205 + t268;
t255 = Icges(2,4) * t259;
t247 = rSges(2,1) * t259 - rSges(2,2) * t257;
t246 = rSges(2,1) * t257 + rSges(2,2) * t259;
t245 = Icges(2,1) * t259 - t290;
t244 = Icges(2,1) * t257 + t255;
t243 = -Icges(2,2) * t257 + t255;
t242 = Icges(2,2) * t259 + t290;
t237 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t236 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t235 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t218 = t260 * rSges(3,3) + (rSges(3,1) * t263 + rSges(3,2) * t264) * t258;
t217 = Icges(3,5) * t260 + (Icges(3,1) * t263 + Icges(3,4) * t264) * t258;
t216 = Icges(3,6) * t260 + (Icges(3,4) * t263 + Icges(3,2) * t264) * t258;
t215 = Icges(3,3) * t260 + (Icges(3,5) * t263 + Icges(3,6) * t264) * t258;
t214 = V_base(5) * rSges(2,3) - t246 * V_base(6) + t279;
t213 = t247 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t203 = t246 * V_base(4) - t247 * V_base(5) + t275;
t202 = qJD(4) * t230 + t227;
t198 = t231 * rSges(4,1) - t230 * rSges(4,2) - rSges(4,3) * t286;
t197 = Icges(4,1) * t231 - Icges(4,4) * t230 - Icges(4,5) * t286;
t196 = Icges(4,4) * t231 - Icges(4,2) * t230 - Icges(4,6) * t286;
t195 = Icges(4,5) * t231 - Icges(4,6) * t230 - Icges(4,3) * t286;
t194 = rSges(3,1) * t226 - rSges(3,2) * t225 + rSges(3,3) * t289;
t193 = rSges(3,1) * t224 - rSges(3,2) * t223 - rSges(3,3) * t288;
t192 = Icges(3,1) * t226 - Icges(3,4) * t225 + Icges(3,5) * t289;
t191 = Icges(3,1) * t224 - Icges(3,4) * t223 - Icges(3,5) * t288;
t190 = Icges(3,4) * t226 - Icges(3,2) * t225 + Icges(3,6) * t289;
t189 = Icges(3,4) * t224 - Icges(3,2) * t223 - Icges(3,6) * t288;
t188 = Icges(3,5) * t226 - Icges(3,6) * t225 + Icges(3,3) * t289;
t187 = Icges(3,5) * t224 - Icges(3,6) * t223 - Icges(3,3) * t288;
t180 = qJD(4) * t208 + t205;
t179 = qJD(4) * t206 + t204;
t174 = rSges(5,1) * t211 - rSges(5,2) * t210 + rSges(5,3) * t230;
t165 = rSges(4,1) * t209 - rSges(4,2) * t208 + rSges(4,3) * t225;
t164 = rSges(4,1) * t207 - rSges(4,2) * t206 + rSges(4,3) * t223;
t163 = Icges(4,1) * t209 - Icges(4,4) * t208 + Icges(4,5) * t225;
t162 = Icges(4,1) * t207 - Icges(4,4) * t206 + Icges(4,5) * t223;
t161 = Icges(4,4) * t209 - Icges(4,2) * t208 + Icges(4,6) * t225;
t160 = Icges(4,4) * t207 - Icges(4,2) * t206 + Icges(4,6) * t223;
t159 = Icges(4,5) * t209 - Icges(4,6) * t208 + Icges(4,3) * t225;
t158 = Icges(4,5) * t207 - Icges(4,6) * t206 + Icges(4,3) * t223;
t154 = -t193 * t250 + t218 * t238 + t273;
t153 = t194 * t250 - t218 * t239 + t271;
t152 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t208;
t150 = rSges(5,1) * t182 - rSges(5,2) * t181 + rSges(5,3) * t206;
t136 = t193 * t239 - t194 * t238 + t272;
t135 = -t164 * t227 + t198 * t204 + t270;
t134 = t165 * t227 - t198 * t205 + t268;
t133 = t164 * t205 - t165 * t204 + t269;
t132 = -t150 * t202 + t174 * t179 + t267;
t131 = t152 * t202 - t174 * t180 + t265;
t130 = t150 * t180 - t152 * t179 + t266;
t129 = qJD(5) * t183 + t179 * t281 - t202 * t283 + t267;
t128 = qJD(5) * t181 - t180 * t281 + t202 * t282 + t265;
t127 = qJD(5) * t210 - t179 * t282 + t180 * t283 + t266;
t1 = m(1) * (t235 ^ 2 + t236 ^ 2 + t237 ^ 2) / 0.2e1 + m(2) * (t203 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + t239 * ((t188 * t289 - t190 * t225 + t192 * t226) * t239 + (t187 * t289 - t189 * t225 + t191 * t226) * t238 + (t215 * t289 - t216 * t225 + t217 * t226) * t250) / 0.2e1 + t238 * ((-t188 * t288 - t190 * t223 + t192 * t224) * t239 + (-t187 * t288 - t189 * t223 + t191 * t224) * t238 + (-t215 * t288 - t216 * t223 + t217 * t224) * t250) / 0.2e1 + t250 * ((t187 * t238 + t188 * t239 + t215 * t250) * t260 + ((t190 * t264 + t192 * t263) * t239 + (t189 * t264 + t191 * t263) * t238 + (t216 * t264 + t217 * t263) * t250) * t258) / 0.2e1 + m(4) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + t205 * ((t159 * t225 - t161 * t208 + t163 * t209) * t205 + (t158 * t225 - t160 * t208 + t162 * t209) * t204 + (t195 * t225 - t196 * t208 + t197 * t209) * t227) / 0.2e1 + t204 * ((t159 * t223 - t161 * t206 + t163 * t207) * t205 + (t158 * t223 - t160 * t206 + t162 * t207) * t204 + (t195 * t223 - t196 * t206 + t197 * t207) * t227) / 0.2e1 + t227 * ((-t159 * t286 - t230 * t161 + t231 * t163) * t205 + (-t158 * t286 - t230 * t160 + t231 * t162) * t204 + (-t195 * t286 - t230 * t196 + t231 * t197) * t227) / 0.2e1 + m(5) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(6) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + ((t181 * t299 + t182 * t297 + t206 * t298) * t202 + (t181 * t304 + t182 * t300 + t206 * t302) * t180 + (t305 * t181 + t301 * t182 + t303 * t206) * t179) * t179 / 0.2e1 + ((t183 * t299 + t184 * t297 + t208 * t298) * t202 + (t304 * t183 + t300 * t184 + t302 * t208) * t180 + (t183 * t305 + t184 * t301 + t208 * t303) * t179) * t180 / 0.2e1 + ((t299 * t210 + t297 * t211 + t298 * t230) * t202 + (t210 * t304 + t211 * t300 + t230 * t302) * t180 + (t210 * t305 + t211 * t301 + t230 * t303) * t179) * t202 / 0.2e1 + ((-t242 * t257 + t244 * t259 + Icges(1,4)) * V_base(5) + (-t243 * t257 + t245 * t259 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t242 * t259 + t244 * t257 + Icges(1,2)) * V_base(5) + (t243 * t259 + t245 * t257 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t257 + Icges(2,6) * t259 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t259 - Icges(2,6) * t257 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
