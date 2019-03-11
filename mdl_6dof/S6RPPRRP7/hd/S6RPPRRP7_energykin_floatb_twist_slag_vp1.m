% Calculate kinetic energy for
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:31
% EndTime: 2019-03-09 02:12:34
% DurationCPUTime: 2.68s
% Computational Cost: add. (1369->294), mult. (1711->390), div. (0->0), fcn. (1549->8), ass. (0->154)
t317 = Icges(2,4) + Icges(3,6);
t316 = Icges(2,1) + Icges(3,2);
t315 = Icges(6,1) + Icges(7,1);
t314 = -Icges(3,4) + Icges(2,5);
t313 = Icges(6,4) + Icges(7,4);
t312 = Icges(3,5) - Icges(2,6);
t311 = Icges(7,5) + Icges(6,5);
t310 = Icges(2,2) + Icges(3,3);
t309 = Icges(6,2) + Icges(7,2);
t308 = Icges(7,6) + Icges(6,6);
t307 = Icges(7,3) + Icges(6,3);
t220 = cos(qJ(1));
t306 = t317 * t220;
t218 = sin(qJ(1));
t305 = t317 * t218;
t212 = pkin(9) + qJ(4);
t201 = sin(t212);
t219 = cos(qJ(5));
t262 = t219 * t220;
t217 = sin(qJ(5));
t264 = t218 * t217;
t158 = -t201 * t264 + t262;
t263 = t218 * t219;
t266 = t217 * t220;
t159 = t201 * t263 + t266;
t202 = cos(t212);
t268 = t202 * t218;
t304 = t308 * t158 + t311 * t159 - t307 * t268;
t160 = t201 * t266 + t263;
t161 = -t201 * t262 + t264;
t267 = t202 * t220;
t303 = t308 * t160 + t311 * t161 + t307 * t267;
t302 = t309 * t158 + t313 * t159 - t308 * t268;
t301 = t309 * t160 + t313 * t161 + t308 * t267;
t300 = t313 * t158 + t315 * t159 - t311 * t268;
t299 = t313 * t160 + t315 * t161 + t311 * t267;
t298 = (-t308 * t217 + t311 * t219) * t202 + t307 * t201;
t297 = (-t309 * t217 + t313 * t219) * t202 + t308 * t201;
t296 = (-t313 * t217 + t315 * t219) * t202 + t311 * t201;
t295 = t316 * t218 + t306;
t294 = t316 * t220 - t305;
t293 = t314 * t218 - t312 * t220;
t292 = t312 * t218 + t314 * t220;
t203 = V_base(6) + qJD(1);
t269 = qJ(3) * t220;
t291 = qJD(3) * t218 + t203 * t269;
t214 = cos(pkin(9));
t213 = sin(pkin(9));
t275 = Icges(4,4) * t213;
t240 = Icges(4,2) * t214 + t275;
t147 = Icges(4,6) * t218 - t220 * t240;
t274 = Icges(4,4) * t214;
t242 = Icges(4,1) * t213 + t274;
t149 = Icges(4,5) * t218 - t220 * t242;
t290 = t147 * t214 + t149 * t213 - t310 * t220 - t305;
t146 = Icges(4,6) * t220 + t218 * t240;
t148 = Icges(4,5) * t220 + t218 * t242;
t289 = t146 * t214 + t148 * t213 + t310 * t218 - t306;
t279 = pkin(5) * t219;
t288 = -qJ(6) * t202 + t201 * t279;
t273 = Icges(5,4) * t201;
t239 = Icges(5,2) * t202 + t273;
t135 = Icges(5,6) * t220 + t218 * t239;
t136 = Icges(5,6) * t218 - t220 * t239;
t272 = Icges(5,4) * t202;
t241 = Icges(5,1) * t201 + t272;
t137 = Icges(5,5) * t220 + t218 * t241;
t138 = Icges(5,5) * t218 - t220 * t241;
t164 = -Icges(5,2) * t201 + t272;
t165 = Icges(5,1) * t202 - t273;
t196 = qJD(4) * t218 + V_base(5);
t197 = qJD(4) * t220 + V_base(4);
t287 = (t135 * t202 + t137 * t201) * t197 + (t136 * t202 + t138 * t201) * t196 + (t164 * t202 + t165 * t201) * t203;
t286 = -pkin(2) - pkin(6);
t281 = pkin(3) * t213;
t280 = pkin(3) * t214;
t277 = rSges(7,1) * t159 + rSges(7,2) * t158 - rSges(7,3) * t268 + pkin(5) * t266 + t288 * t218;
t265 = t218 * qJ(3);
t260 = t161 * rSges(7,1) + t160 * rSges(7,2) + rSges(7,3) * t267 + pkin(5) * t264 - t288 * t220;
t259 = (rSges(7,1) * t219 - rSges(7,2) * t217 + t279) * t202 + (qJ(6) + rSges(7,3)) * t201;
t258 = qJD(2) * t220;
t257 = qJD(5) * t202;
t256 = qJD(6) * t202;
t192 = pkin(1) * t220 + t218 * qJ(2);
t255 = t203 * t192 + V_base(2);
t189 = t218 * pkin(1) - qJ(2) * t220;
t254 = V_base(4) * t189 + V_base(3);
t253 = V_base(5) * pkin(6) + V_base(1);
t250 = V_base(4) * t265 + t254;
t249 = qJD(2) * t218 + t253;
t248 = -t189 - t265;
t247 = -t192 - t269;
t246 = pkin(4) * t201 - pkin(8) * t202;
t154 = pkin(7) * t218 - t220 * t281;
t245 = -t154 + t248;
t244 = rSges(4,1) * t213 + rSges(4,2) * t214;
t243 = rSges(5,1) * t201 + rSges(5,2) * t202;
t238 = Icges(4,5) * t213 + Icges(4,6) * t214;
t237 = Icges(5,5) * t201 + Icges(5,6) * t202;
t174 = -Icges(4,2) * t213 + t274;
t175 = Icges(4,1) * t214 - t275;
t231 = t174 * t214 + t175 * t213;
t230 = t255 - t258;
t229 = V_base(5) * pkin(2) + qJD(3) * t220 + t249;
t228 = V_base(5) * t280 + t229;
t227 = (Icges(5,3) * t220 + t218 * t237) * t197 + (Icges(5,3) * t218 - t220 * t237) * t196 + (Icges(5,5) * t202 - Icges(5,6) * t201) * t203;
t226 = (Icges(4,3) * t220 + t218 * t238) * V_base(4) + (Icges(4,3) * t218 - t220 * t238) * V_base(5) + (Icges(4,5) * t214 - Icges(4,6) * t213) * t203;
t155 = pkin(7) * t220 + t218 * t281;
t225 = V_base(4) * t154 + (-t155 + t247) * V_base(5) + t250;
t224 = t203 * t155 + (-t280 + t286) * V_base(4) + t255 + t291;
t153 = t246 * t220;
t167 = pkin(4) * t202 + pkin(8) * t201;
t223 = t196 * t167 + (t153 + t245) * t203 + t228;
t152 = t246 * t218;
t222 = -t196 * t152 - t197 * t153 + t225;
t221 = t203 * t152 - t197 * t167 + t224;
t194 = rSges(2,1) * t220 - t218 * rSges(2,2);
t193 = -rSges(3,2) * t220 + t218 * rSges(3,3);
t191 = t218 * rSges(2,1) + rSges(2,2) * t220;
t190 = -t218 * rSges(3,2) - rSges(3,3) * t220;
t176 = rSges(4,1) * t214 - rSges(4,2) * t213;
t172 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t171 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t170 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t169 = qJD(5) * t201 + t203;
t166 = rSges(5,1) * t202 - rSges(5,2) * t201;
t157 = -t218 * t257 + t197;
t156 = t220 * t257 + t196;
t151 = t218 * rSges(4,3) - t220 * t244;
t150 = rSges(4,3) * t220 + t218 * t244;
t141 = t218 * rSges(5,3) - t220 * t243;
t140 = rSges(5,3) * t220 + t218 * t243;
t131 = rSges(6,3) * t201 + (rSges(6,1) * t219 - rSges(6,2) * t217) * t202;
t129 = V_base(5) * rSges(2,3) - t191 * t203 + t253;
t128 = t194 * t203 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t120 = t191 * V_base(4) - t194 * V_base(5) + V_base(3);
t118 = V_base(5) * rSges(3,1) + (-t189 - t190) * t203 + t249;
t117 = t203 * t193 + (-rSges(3,1) - pkin(6)) * V_base(4) + t230;
t116 = t161 * rSges(6,1) + t160 * rSges(6,2) + rSges(6,3) * t267;
t114 = rSges(6,1) * t159 + rSges(6,2) * t158 - rSges(6,3) * t268;
t98 = t190 * V_base(4) + (-t192 - t193) * V_base(5) + t254;
t97 = t176 * V_base(5) + (-t151 + t248) * t203 + t229;
t96 = t203 * t150 + (-t176 + t286) * V_base(4) + t230 + t291;
t95 = V_base(4) * t151 + (-t150 + t247) * V_base(5) + t250;
t94 = t166 * t196 + (-t141 + t245) * t203 + t228;
t93 = t203 * t140 - t197 * t166 + t224 - t258;
t92 = -t196 * t140 + t197 * t141 + t225;
t91 = -t116 * t169 + t131 * t156 + t223;
t90 = t169 * t114 - t157 * t131 + t221 - t258;
t89 = -t156 * t114 + t157 * t116 + t222;
t88 = t156 * t259 - t169 * t260 - t218 * t256 + t223;
t87 = (-qJD(2) + t256) * t220 + t277 * t169 - t259 * t157 + t221;
t86 = qJD(6) * t201 - t156 * t277 + t157 * t260 + t222;
t1 = t197 * (t287 * t218 + t227 * t220) / 0.2e1 + t196 * (t227 * t218 - t287 * t220) / 0.2e1 + m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t120 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t118 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + ((t297 * t160 + t296 * t161 + t298 * t267) * t169 + (t302 * t160 + t300 * t161 + t304 * t267) * t157 + (t301 * t160 + t299 * t161 + t303 * t267) * t156) * t156 / 0.2e1 + ((t297 * t158 + t296 * t159 - t298 * t268) * t169 + (t302 * t158 + t300 * t159 - t304 * t268) * t157 + (t301 * t158 + t299 * t159 - t303 * t268) * t156) * t157 / 0.2e1 + (((-t297 * t217 + t296 * t219) * t169 + (-t302 * t217 + t300 * t219) * t157 + (-t301 * t217 + t299 * t219) * t156) * t202 + (t303 * t156 + t304 * t157 + t298 * t169) * t201) * t169 / 0.2e1 + ((-t135 * t201 + t137 * t202) * t197 + (-t136 * t201 + t138 * t202) * t196 + (-t147 * t213 + t149 * t214 + t293) * V_base(5) + (-t146 * t213 + t148 * t214 + t292) * V_base(4) + (-t164 * t201 + t165 * t202 - t174 * t213 + t175 * t214 + Icges(3,1) + Icges(2,3)) * t203) * t203 / 0.2e1 + (t226 * t220 + (t231 * t218 + t292) * t203 + (t290 * t218 + t295 * t220 + Icges(1,4)) * V_base(5) + (t289 * t218 + t294 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t226 * t218 + (-t231 * t220 + t293) * t203 + (t295 * t218 - t290 * t220 + Icges(1,2)) * V_base(5) + (t294 * t218 - t289 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
