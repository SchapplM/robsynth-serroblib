% Calculate kinetic energy for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:34:58
% EndTime: 2019-03-09 04:35:00
% DurationCPUTime: 2.63s
% Computational Cost: add. (1816->269), mult. (2006->366), div. (0->0), fcn. (1924->8), ass. (0->131)
t293 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t292 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t291 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t290 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t289 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t288 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t287 = rSges(7,1) + pkin(5);
t286 = rSges(7,3) + qJ(6);
t226 = qJ(1) + pkin(9);
t220 = sin(t226);
t221 = cos(t226);
t230 = cos(qJ(4));
t227 = sin(qJ(4));
t231 = cos(qJ(3));
t261 = t227 * t231;
t162 = t220 * t261 + t221 * t230;
t259 = t230 * t231;
t163 = t220 * t259 - t221 * t227;
t228 = sin(qJ(3));
t263 = t220 * t228;
t283 = t292 * t162 + t293 * t163 - t291 * t263;
t164 = -t220 * t230 + t221 * t261;
t165 = t220 * t227 + t221 * t259;
t262 = t221 * t228;
t282 = t292 * t164 + t293 * t165 - t291 * t262;
t281 = t290 * t162 + t292 * t163 - t289 * t263;
t280 = t290 * t164 + t292 * t165 - t289 * t262;
t279 = -t289 * t162 - t291 * t163 - t288 * t263;
t278 = -t289 * t164 - t291 * t165 - t288 * t262;
t277 = t288 * t231 + (-t289 * t227 - t291 * t230) * t228;
t276 = t289 * t231 + (t290 * t227 + t292 * t230) * t228;
t275 = t291 * t231 + (t292 * t227 + t293 * t230) * t228;
t229 = sin(qJ(1));
t270 = pkin(1) * t229;
t232 = cos(qJ(1));
t269 = pkin(1) * t232;
t268 = -pkin(6) - qJ(2);
t267 = Icges(2,4) * t229;
t266 = Icges(3,4) * t220;
t265 = Icges(4,4) * t228;
t264 = Icges(4,4) * t231;
t260 = t228 * t230;
t258 = rSges(7,2) * t162 + t286 * t163 + t287 * t263;
t257 = rSges(7,2) * t164 + t286 * t165 + t287 * t262;
t256 = (rSges(7,2) * t227 + rSges(7,3) * t230) * t228 + qJ(6) * t260 - t287 * t231;
t255 = qJD(4) * t228;
t222 = V_base(6) + qJD(1);
t254 = t222 * t269 + V_base(2);
t253 = V_base(5) * pkin(6) + V_base(1);
t198 = qJD(3) * t220 + V_base(4);
t191 = pkin(2) * t220 - pkin(7) * t221;
t250 = -t191 - t270;
t249 = V_base(5) * qJ(2) + t253;
t248 = V_base(4) * t270 + qJD(2) + V_base(3);
t247 = pkin(3) * t231 + pkin(8) * t228;
t197 = -qJD(3) * t221 + V_base(5);
t246 = rSges(4,1) * t231 - rSges(4,2) * t228;
t245 = Icges(4,1) * t231 - t265;
t244 = -Icges(4,2) * t228 + t264;
t243 = Icges(4,5) * t231 - Icges(4,6) * t228;
t242 = (-Icges(4,3) * t221 + t220 * t243) * t197 + (Icges(4,3) * t220 + t221 * t243) * t198 + (Icges(4,5) * t228 + Icges(4,6) * t231) * t222;
t192 = pkin(2) * t221 + pkin(7) * t220;
t241 = t222 * t192 + t268 * V_base(4) + t254;
t179 = t247 * t220;
t212 = pkin(3) * t228 - pkin(8) * t231;
t240 = t197 * t212 + (-t179 + t250) * t222 + t249;
t239 = V_base(4) * t191 + (-t192 - t269) * V_base(5) + t248;
t180 = t247 * t221;
t238 = t222 * t180 - t198 * t212 + t241;
t160 = t220 * t255 + t197;
t182 = (pkin(4) * t230 + qJ(5) * t227) * t228;
t237 = qJD(5) * t164 + t160 * t182 + t240;
t140 = pkin(4) * t165 + qJ(5) * t164;
t208 = -qJD(4) * t231 + t222;
t236 = qJD(5) * t162 + t208 * t140 + t238;
t235 = t198 * t179 - t197 * t180 + t239;
t139 = pkin(4) * t163 + qJ(5) * t162;
t161 = t221 * t255 + t198;
t234 = qJD(5) * t228 * t227 + t161 * t139 + t235;
t148 = -Icges(4,6) * t221 + t220 * t244;
t149 = Icges(4,6) * t220 + t221 * t244;
t150 = -Icges(4,5) * t221 + t220 * t245;
t151 = Icges(4,5) * t220 + t221 * t245;
t202 = Icges(4,2) * t231 + t265;
t205 = Icges(4,1) * t228 + t264;
t233 = (-t149 * t228 + t151 * t231) * t198 + (-t148 * t228 + t150 * t231) * t197 + (-t202 * t228 + t205 * t231) * t222;
t224 = Icges(2,4) * t232;
t219 = Icges(3,4) * t221;
t211 = rSges(2,1) * t232 - t229 * rSges(2,2);
t210 = t229 * rSges(2,1) + rSges(2,2) * t232;
t209 = rSges(4,1) * t228 + rSges(4,2) * t231;
t207 = Icges(2,1) * t232 - t267;
t206 = Icges(2,1) * t229 + t224;
t204 = -Icges(2,2) * t229 + t224;
t203 = Icges(2,2) * t232 + t267;
t196 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t195 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t194 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t190 = rSges(3,1) * t221 - rSges(3,2) * t220;
t189 = rSges(3,1) * t220 + rSges(3,2) * t221;
t188 = Icges(3,1) * t221 - t266;
t187 = Icges(3,1) * t220 + t219;
t186 = -Icges(3,2) * t220 + t219;
t185 = Icges(3,2) * t221 + t266;
t177 = -rSges(6,1) * t231 + (-rSges(6,2) * t230 + rSges(6,3) * t227) * t228;
t175 = -rSges(5,3) * t231 + (rSges(5,1) * t230 - rSges(5,2) * t227) * t228;
t155 = rSges(4,3) * t220 + t221 * t246;
t154 = -rSges(4,3) * t221 + t220 * t246;
t153 = V_base(5) * rSges(2,3) - t210 * t222 + t253;
t152 = t211 * t222 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t210 * V_base(4) - t211 * V_base(5) + V_base(3);
t138 = V_base(5) * rSges(3,3) + (-t189 - t270) * t222 + t249;
t137 = t190 * t222 + (-rSges(3,3) + t268) * V_base(4) + t254;
t136 = V_base(4) * t189 + (-t190 - t269) * V_base(5) + t248;
t135 = rSges(5,1) * t165 - rSges(5,2) * t164 + rSges(5,3) * t262;
t134 = rSges(5,1) * t163 - rSges(5,2) * t162 + rSges(5,3) * t263;
t133 = rSges(6,1) * t262 - rSges(6,2) * t165 + rSges(6,3) * t164;
t131 = rSges(6,1) * t263 - rSges(6,2) * t163 + rSges(6,3) * t162;
t109 = t197 * t209 + (-t154 + t250) * t222 + t249;
t108 = t155 * t222 - t198 * t209 + t241;
t107 = t198 * t154 - t197 * t155 + t239;
t106 = -t134 * t208 + t160 * t175 + t240;
t105 = t135 * t208 - t161 * t175 + t238;
t104 = t161 * t134 - t160 * t135 + t235;
t103 = t160 * t177 + (-t131 - t139) * t208 + t237;
t102 = t133 * t208 + (-t177 - t182) * t161 + t236;
t101 = qJD(6) * t165 + t256 * t160 + (-t139 - t258) * t208 + t237;
t100 = qJD(6) * t163 + t257 * t208 + (-t182 - t256) * t161 + t236;
t99 = t161 * t131 + (-t133 - t140) * t160 + t234;
t98 = qJD(6) * t260 + t258 * t161 + (-t140 - t257) * t160 + t234;
t1 = t198 * (t242 * t220 + t233 * t221) / 0.2e1 + t197 * (t233 * t220 - t242 * t221) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(2) * (t145 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(1) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + ((t162 * t276 + t163 * t275 + t263 * t277) * t208 + (t162 * t280 + t163 * t282 + t263 * t278) * t161 + (t281 * t162 + t283 * t163 + t279 * t263) * t160) * t160 / 0.2e1 + ((t164 * t276 + t165 * t275 + t262 * t277) * t208 + (t280 * t164 + t282 * t165 + t278 * t262) * t161 + (t164 * t281 + t165 * t283 + t262 * t279) * t160) * t161 / 0.2e1 + ((-t160 * t279 - t161 * t278 - t208 * t277) * t231 + ((t227 * t276 + t230 * t275) * t208 + (t227 * t280 + t230 * t282) * t161 + (t227 * t281 + t230 * t283) * t160) * t228) * t208 / 0.2e1 + ((t149 * t231 + t151 * t228) * t198 + (t148 * t231 + t150 * t228) * t197 + (t231 * t202 + t228 * t205 + Icges(2,3) + Icges(3,3)) * t222) * t222 / 0.2e1 + ((-t185 * t220 + t187 * t221 - t229 * t203 + t206 * t232 + Icges(1,4)) * V_base(5) + (-t220 * t186 + t221 * t188 - t229 * t204 + t232 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t221 * t185 + t220 * t187 + t232 * t203 + t229 * t206 + Icges(1,2)) * V_base(5) + (t186 * t221 + t188 * t220 + t204 * t232 + t229 * t207 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t222 * (Icges(2,5) * t229 + Icges(3,5) * t220 + Icges(2,6) * t232 + Icges(3,6) * t221) + V_base(4) * t222 * (Icges(2,5) * t232 + Icges(3,5) * t221 - Icges(2,6) * t229 - Icges(3,6) * t220) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
