% Calculate kinetic energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:17
% EndTime: 2019-12-31 19:23:19
% DurationCPUTime: 1.95s
% Computational Cost: add. (1289->240), mult. (3030->320), div. (0->0), fcn. (3391->8), ass. (0->120)
t301 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t300 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t299 = -Icges(4,4) - Icges(5,6) + Icges(6,6);
t298 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t297 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t296 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t295 = rSges(6,1) + pkin(4);
t294 = rSges(6,3) + qJ(5);
t237 = sin(pkin(5));
t239 = sin(qJ(1));
t238 = sin(qJ(2));
t240 = cos(qJ(2));
t279 = cos(pkin(8));
t280 = cos(pkin(5));
t258 = t280 * t279;
t278 = sin(pkin(8));
t245 = t238 * t258 + t240 * t278;
t241 = cos(qJ(1));
t260 = t241 * t279;
t180 = t237 * t260 + t239 * t245;
t257 = t280 * t278;
t252 = t238 * t257;
t261 = t239 * t279;
t262 = t237 * t278;
t181 = -t239 * t252 + t240 * t261 - t241 * t262;
t274 = t237 * t238;
t205 = t239 * t274 - t241 * t280;
t293 = t299 * t180 + t301 * t181 - t298 * t205;
t182 = -t237 * t261 + t241 * t245;
t183 = t239 * t262 + t240 * t260 - t241 * t252;
t206 = t239 * t280 + t241 * t274;
t292 = t299 * t182 + t301 * t183 - t298 * t206;
t291 = t296 * t180 + t299 * t181 + t297 * t205;
t290 = t296 * t182 + t299 * t183 + t297 * t206;
t289 = t297 * t180 - t298 * t181 + t300 * t205;
t288 = t297 * t182 - t298 * t183 + t300 * t206;
t203 = t238 * t278 - t240 * t258;
t204 = t238 * t279 + t240 * t257;
t273 = t237 * t240;
t287 = t299 * t203 + t301 * t204 + t298 * t273;
t286 = t296 * t203 + t299 * t204 - t297 * t273;
t285 = t297 * t203 - t298 * t204 - t300 * t273;
t281 = pkin(2) * t240;
t277 = Icges(2,4) * t239;
t276 = Icges(3,4) * t238;
t275 = Icges(3,4) * t240;
t272 = rSges(6,2) * t180 + t294 * t181 + t295 * t205;
t271 = rSges(6,2) * t182 + t294 * t183 + t295 * t206;
t158 = pkin(3) * t183 + qJ(4) * t182;
t190 = qJ(3) * t206 + t241 * t281;
t270 = -t158 - t190;
t269 = rSges(6,2) * t203 + t294 * t204 - t295 * t273;
t176 = pkin(3) * t204 + qJ(4) * t203;
t208 = pkin(2) * t238 - qJ(3) * t273;
t268 = -t176 - t208;
t189 = qJ(3) * t205 + t239 * t281;
t230 = t239 * pkin(1) - pkin(7) * t241;
t267 = -t189 - t230;
t266 = V_base(5) * pkin(6) + V_base(1);
t157 = pkin(3) * t181 + qJ(4) * t180;
t263 = -t157 + t267;
t232 = -qJD(2) * t241 + V_base(5);
t259 = qJD(3) * t206 + t232 * t208 + t266;
t256 = rSges(3,1) * t240 - rSges(3,2) * t238;
t255 = Icges(3,1) * t240 - t276;
t254 = -Icges(3,2) * t238 + t275;
t253 = Icges(3,5) * t240 - Icges(3,6) * t238;
t231 = pkin(1) * t241 + t239 * pkin(7);
t234 = V_base(6) + qJD(1);
t251 = -V_base(4) * pkin(6) + t234 * t231 + V_base(2);
t250 = V_base(4) * t230 - t231 * V_base(5) + V_base(3);
t249 = qJD(4) * t182 + t232 * t176 + t259;
t233 = qJD(2) * t239 + V_base(4);
t248 = (-Icges(3,3) * t241 + t239 * t253) * t232 + (Icges(3,3) * t239 + t241 * t253) * t233 + (Icges(3,5) * t238 + Icges(3,6) * t240) * t234;
t247 = qJD(3) * t205 + t234 * t190 + t251;
t246 = qJD(4) * t180 + t234 * t158 + t247;
t244 = -qJD(3) * t273 + t233 * t189 + t250;
t243 = qJD(4) * t203 + t233 * t157 + t244;
t194 = -Icges(3,6) * t241 + t239 * t254;
t195 = Icges(3,6) * t239 + t241 * t254;
t196 = -Icges(3,5) * t241 + t239 * t255;
t197 = Icges(3,5) * t239 + t241 * t255;
t218 = Icges(3,2) * t240 + t276;
t221 = Icges(3,1) * t238 + t275;
t242 = (-t195 * t238 + t197 * t240) * t233 + (-t194 * t238 + t196 * t240) * t232 + (-t218 * t238 + t221 * t240) * t234;
t235 = Icges(2,4) * t241;
t226 = rSges(2,1) * t241 - t239 * rSges(2,2);
t225 = t239 * rSges(2,1) + rSges(2,2) * t241;
t224 = rSges(3,1) * t238 + rSges(3,2) * t240;
t223 = Icges(2,1) * t241 - t277;
t222 = Icges(2,1) * t239 + t235;
t220 = -Icges(2,2) * t239 + t235;
t219 = Icges(2,2) * t241 + t277;
t214 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t213 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t212 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t199 = t239 * rSges(3,3) + t241 * t256;
t198 = -rSges(3,3) * t241 + t239 * t256;
t187 = V_base(5) * rSges(2,3) - t225 * t234 + t266;
t186 = t226 * t234 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t185 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t175 = rSges(4,1) * t204 - rSges(4,2) * t203 - rSges(4,3) * t273;
t174 = -rSges(5,1) * t273 - rSges(5,2) * t204 + rSges(5,3) * t203;
t160 = t224 * t232 + (-t198 - t230) * t234 + t266;
t159 = t199 * t234 - t224 * t233 + t251;
t154 = t198 * t233 - t199 * t232 + t250;
t153 = rSges(5,1) * t206 - rSges(5,2) * t183 + rSges(5,3) * t182;
t151 = rSges(5,1) * t205 - rSges(5,2) * t181 + rSges(5,3) * t180;
t149 = rSges(4,1) * t183 - rSges(4,2) * t182 + rSges(4,3) * t206;
t148 = rSges(4,1) * t181 - rSges(4,2) * t180 + rSges(4,3) * t205;
t129 = t175 * t232 + (-t148 + t267) * t234 + t259;
t128 = t149 * t234 + (-t175 - t208) * t233 + t247;
t127 = t148 * t233 + (-t149 - t190) * t232 + t244;
t126 = t174 * t232 + (-t151 + t263) * t234 + t249;
t125 = t153 * t234 + (-t174 + t268) * t233 + t246;
t124 = t151 * t233 + (-t153 + t270) * t232 + t243;
t123 = qJD(5) * t183 + t269 * t232 + (t263 - t272) * t234 + t249;
t122 = qJD(5) * t181 + t271 * t234 + (t268 - t269) * t233 + t246;
t121 = qJD(5) * t204 + t272 * t233 + (t270 - t271) * t232 + t243;
t1 = m(1) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(2) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(3) * (t154 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(4) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(5) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(6) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + ((-t239 * t219 + t222 * t241 + Icges(1,4)) * V_base(5) + (-t239 * t220 + t223 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t219 * t241 + t239 * t222 + Icges(1,2)) * V_base(5) + (t220 * t241 + t239 * t223 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t239 * t242 - t248 * t241 + (t286 * t180 + t287 * t181 + t285 * t205) * t234 + (t290 * t180 + t292 * t181 + t288 * t205) * t233 + (t291 * t180 + t293 * t181 + t289 * t205) * t232) * t232 / 0.2e1 + (t239 * t248 + t241 * t242 + (t286 * t182 + t287 * t183 + t285 * t206) * t234 + (t290 * t182 + t292 * t183 + t288 * t206) * t233 + (t291 * t182 + t293 * t183 + t289 * t206) * t232) * t233 / 0.2e1 + ((t195 * t240 + t197 * t238 + t290 * t203 + t292 * t204 - t288 * t273) * t233 + (t194 * t240 + t196 * t238 + t291 * t203 + t293 * t204 - t289 * t273) * t232 + (t286 * t203 + t287 * t204 + t218 * t240 + t221 * t238 - t285 * t273 + Icges(2,3)) * t234) * t234 / 0.2e1 + t234 * V_base(4) * (Icges(2,5) * t241 - Icges(2,6) * t239) + V_base(5) * t234 * (Icges(2,5) * t239 + Icges(2,6) * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
