% Calculate kinetic energy for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:13
% EndTime: 2019-12-31 21:43:16
% DurationCPUTime: 3.15s
% Computational Cost: add. (1747->314), mult. (3964->451), div. (0->0), fcn. (4646->10), ass. (0->140)
t301 = Icges(4,1) + Icges(5,2);
t300 = Icges(5,1) + Icges(4,3);
t299 = -Icges(4,4) - Icges(5,6);
t298 = Icges(5,4) - Icges(4,5);
t297 = Icges(5,5) - Icges(4,6);
t296 = Icges(4,2) + Icges(5,3);
t251 = cos(pkin(5));
t254 = sin(qJ(1));
t256 = cos(qJ(2));
t274 = t254 * t256;
t253 = sin(qJ(2));
t257 = cos(qJ(1));
t275 = t253 * t257;
t219 = t251 * t275 + t274;
t250 = sin(pkin(5));
t283 = cos(qJ(3));
t268 = t250 * t283;
t282 = sin(qJ(3));
t198 = t219 * t282 + t257 * t268;
t267 = t250 * t282;
t199 = t219 * t283 - t257 * t267;
t273 = t256 * t257;
t276 = t253 * t254;
t218 = -t251 * t273 + t276;
t295 = t296 * t198 + t299 * t199 + t297 * t218;
t221 = -t251 * t276 + t273;
t200 = t221 * t282 - t254 * t268;
t201 = t221 * t283 + t254 * t267;
t220 = t251 * t274 + t275;
t294 = t296 * t200 + t299 * t201 + t297 * t220;
t293 = t297 * t198 - t298 * t199 + t300 * t218;
t292 = t297 * t200 - t298 * t201 + t300 * t220;
t291 = t299 * t198 + t301 * t199 - t298 * t218;
t290 = t299 * t200 + t301 * t201 - t298 * t220;
t216 = -t251 * t283 + t253 * t267;
t217 = t251 * t282 + t253 * t268;
t278 = t250 * t256;
t289 = t296 * t216 + t299 * t217 - t297 * t278;
t288 = t299 * t216 + t301 * t217 + t298 * t278;
t287 = t297 * t216 - t298 * t217 - t300 * t278;
t281 = pkin(7) * t251;
t280 = Icges(2,4) * t254;
t279 = t250 * t254;
t277 = t250 * t257;
t272 = qJD(2) * t250;
t271 = V_base(5) * pkin(6) + V_base(1);
t230 = t254 * t272 + V_base(4);
t247 = V_base(6) + qJD(1);
t197 = qJD(3) * t220 + t230;
t231 = qJD(2) * t251 + t247;
t229 = -t257 * t272 + V_base(5);
t224 = t254 * pkin(1) - pkin(7) * t277;
t266 = -t224 * t247 + V_base(5) * t281 + t271;
t225 = pkin(1) * t257 + pkin(7) * t279;
t265 = V_base(4) * t224 - t225 * V_base(5) + V_base(3);
t196 = qJD(3) * t218 + t229;
t214 = -qJD(3) * t278 + t231;
t264 = t247 * t225 + V_base(2) + (-pkin(6) - t281) * V_base(4);
t190 = pkin(2) * t219 + pkin(8) * t218;
t223 = (pkin(2) * t253 - pkin(8) * t256) * t250;
t263 = -t190 * t231 + t229 * t223 + t266;
t191 = pkin(2) * t221 + pkin(8) * t220;
t262 = t230 * t190 - t191 * t229 + t265;
t188 = pkin(3) * t217 + qJ(4) * t216;
t261 = qJD(4) * t200 + t196 * t188 + t263;
t160 = pkin(3) * t199 + qJ(4) * t198;
t260 = qJD(4) * t216 + t197 * t160 + t262;
t259 = t231 * t191 - t223 * t230 + t264;
t161 = pkin(3) * t201 + qJ(4) * t200;
t258 = qJD(4) * t198 + t214 * t161 + t259;
t255 = cos(qJ(5));
t252 = sin(qJ(5));
t248 = Icges(2,4) * t257;
t239 = rSges(2,1) * t257 - t254 * rSges(2,2);
t238 = t254 * rSges(2,1) + rSges(2,2) * t257;
t237 = Icges(2,1) * t257 - t280;
t236 = Icges(2,1) * t254 + t248;
t235 = -Icges(2,2) * t254 + t248;
t234 = Icges(2,2) * t257 + t280;
t228 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t227 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t226 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t210 = rSges(3,3) * t251 + (rSges(3,1) * t253 + rSges(3,2) * t256) * t250;
t209 = Icges(3,5) * t251 + (Icges(3,1) * t253 + Icges(3,4) * t256) * t250;
t208 = Icges(3,6) * t251 + (Icges(3,4) * t253 + Icges(3,2) * t256) * t250;
t207 = Icges(3,3) * t251 + (Icges(3,5) * t253 + Icges(3,6) * t256) * t250;
t206 = -pkin(4) * t278 + pkin(9) * t217;
t205 = V_base(5) * rSges(2,3) - t238 * t247 + t271;
t204 = t239 * t247 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t202 = t238 * V_base(4) - t239 * V_base(5) + V_base(3);
t195 = t216 * t252 - t255 * t278;
t194 = t216 * t255 + t252 * t278;
t189 = qJD(5) * t217 + t214;
t187 = rSges(3,1) * t221 - rSges(3,2) * t220 + rSges(3,3) * t279;
t186 = t219 * rSges(3,1) - t218 * rSges(3,2) - rSges(3,3) * t277;
t185 = Icges(3,1) * t221 - Icges(3,4) * t220 + Icges(3,5) * t279;
t184 = Icges(3,1) * t219 - Icges(3,4) * t218 - Icges(3,5) * t277;
t183 = Icges(3,4) * t221 - Icges(3,2) * t220 + Icges(3,6) * t279;
t182 = Icges(3,4) * t219 - Icges(3,2) * t218 - Icges(3,6) * t277;
t181 = Icges(3,5) * t221 - Icges(3,6) * t220 + Icges(3,3) * t279;
t180 = Icges(3,5) * t219 - Icges(3,6) * t218 - Icges(3,3) * t277;
t179 = rSges(4,1) * t217 - rSges(4,2) * t216 - rSges(4,3) * t278;
t178 = -rSges(5,1) * t278 - rSges(5,2) * t217 + rSges(5,3) * t216;
t169 = pkin(4) * t220 + pkin(9) * t201;
t168 = pkin(4) * t218 + pkin(9) * t199;
t167 = t200 * t252 + t220 * t255;
t166 = t200 * t255 - t220 * t252;
t165 = t198 * t252 + t218 * t255;
t164 = t198 * t255 - t218 * t252;
t163 = qJD(5) * t201 + t197;
t162 = qJD(5) * t199 + t196;
t158 = rSges(4,1) * t201 - rSges(4,2) * t200 + rSges(4,3) * t220;
t157 = rSges(4,1) * t199 - rSges(4,2) * t198 + rSges(4,3) * t218;
t156 = rSges(5,1) * t220 - rSges(5,2) * t201 + rSges(5,3) * t200;
t155 = rSges(5,1) * t218 - rSges(5,2) * t199 + rSges(5,3) * t198;
t141 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t217;
t140 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t217;
t139 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t217;
t138 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t217;
t136 = -t186 * t231 + t210 * t229 + t266;
t135 = t187 * t231 - t210 * t230 + t264;
t134 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t201;
t133 = rSges(6,1) * t165 + rSges(6,2) * t164 + rSges(6,3) * t199;
t132 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t201;
t131 = Icges(6,1) * t165 + Icges(6,4) * t164 + Icges(6,5) * t199;
t130 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t201;
t129 = Icges(6,4) * t165 + Icges(6,2) * t164 + Icges(6,6) * t199;
t128 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t201;
t127 = Icges(6,5) * t165 + Icges(6,6) * t164 + Icges(6,3) * t199;
t126 = t186 * t230 - t187 * t229 + t265;
t125 = -t157 * t214 + t179 * t196 + t263;
t124 = t158 * t214 - t179 * t197 + t259;
t123 = t157 * t197 - t158 * t196 + t262;
t122 = t178 * t196 + (-t155 - t160) * t214 + t261;
t121 = t156 * t214 + (-t178 - t188) * t197 + t258;
t120 = t155 * t197 + (-t156 - t161) * t196 + t260;
t119 = -t133 * t189 + t141 * t162 + t196 * t206 + (-t160 - t168) * t214 + t261;
t118 = t134 * t189 - t141 * t163 + t169 * t214 + (-t188 - t206) * t197 + t258;
t117 = t133 * t163 - t134 * t162 + t168 * t197 + (-t161 - t169) * t196 + t260;
t1 = m(1) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + m(2) * (t202 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(3) * (t126 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + t230 * ((t181 * t279 - t220 * t183 + t221 * t185) * t230 + (t180 * t279 - t182 * t220 + t184 * t221) * t229 + (t207 * t279 - t208 * t220 + t209 * t221) * t231) / 0.2e1 + t229 * ((-t181 * t277 - t218 * t183 + t219 * t185) * t230 + (-t180 * t277 - t218 * t182 + t219 * t184) * t229 + (-t207 * t277 - t218 * t208 + t219 * t209) * t231) / 0.2e1 + t231 * ((t180 * t229 + t181 * t230 + t207 * t231) * t251 + ((t183 * t256 + t185 * t253) * t230 + (t182 * t256 + t184 * t253) * t229 + (t208 * t256 + t209 * t253) * t231) * t250) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + t163 * ((t201 * t128 + t166 * t130 + t167 * t132) * t163 + (t127 * t201 + t129 * t166 + t131 * t167) * t162 + (t138 * t201 + t139 * t166 + t140 * t167) * t189) / 0.2e1 + t162 * ((t128 * t199 + t130 * t164 + t132 * t165) * t163 + (t199 * t127 + t164 * t129 + t165 * t131) * t162 + (t138 * t199 + t139 * t164 + t140 * t165) * t189) / 0.2e1 + t189 * ((t128 * t217 + t130 * t194 + t132 * t195) * t163 + (t127 * t217 + t129 * t194 + t131 * t195) * t162 + (t217 * t138 + t194 * t139 + t195 * t140) * t189) / 0.2e1 + ((t198 * t289 + t199 * t288 + t218 * t287) * t214 + (t198 * t294 + t199 * t290 + t218 * t292) * t197 + (t295 * t198 + t291 * t199 + t293 * t218) * t196) * t196 / 0.2e1 + ((t200 * t289 + t201 * t288 + t220 * t287) * t214 + (t200 * t294 + t290 * t201 + t292 * t220) * t197 + (t200 * t295 + t291 * t201 + t293 * t220) * t196) * t197 / 0.2e1 + ((t216 * t289 + t217 * t288 - t278 * t287) * t214 + (t216 * t294 + t217 * t290 - t278 * t292) * t197 + (t216 * t295 + t291 * t217 - t293 * t278) * t196) * t214 / 0.2e1 + ((-t254 * t234 + t236 * t257 + Icges(1,4)) * V_base(5) + (-t254 * t235 + t257 * t237 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t257 * t234 + t254 * t236 + Icges(1,2)) * V_base(5) + (t235 * t257 + t254 * t237 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t254 + Icges(2,6) * t257) * V_base(5) + (Icges(2,5) * t257 - Icges(2,6) * t254) * V_base(4) + Icges(2,3) * t247 / 0.2e1) * t247;
T = t1;
