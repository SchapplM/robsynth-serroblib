% Calculate kinetic energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:03
% EndTime: 2019-03-09 02:33:06
% DurationCPUTime: 3.15s
% Computational Cost: add. (1446->313), mult. (1486->440), div. (0->0), fcn. (1268->10), ass. (0->167)
t310 = Icges(2,4) + Icges(3,6);
t309 = Icges(2,1) + Icges(3,2);
t308 = -Icges(3,4) + Icges(2,5);
t307 = Icges(3,5) - Icges(2,6);
t306 = Icges(2,2) + Icges(3,3);
t225 = cos(qJ(1));
t305 = t310 * t225;
t223 = sin(qJ(1));
t304 = t310 * t223;
t303 = t309 * t223 + t305;
t302 = t309 * t225 - t304;
t301 = t308 * t223 - t307 * t225;
t300 = t307 * t223 + t308 * t225;
t220 = cos(pkin(10));
t219 = sin(pkin(10));
t289 = Icges(4,4) * t219;
t253 = Icges(4,2) * t220 + t289;
t144 = Icges(4,6) * t223 - t225 * t253;
t288 = Icges(4,4) * t220;
t256 = Icges(4,1) * t219 + t288;
t146 = Icges(4,5) * t223 - t225 * t256;
t299 = t144 * t220 + t146 * t219 - t306 * t225 - t304;
t143 = Icges(4,6) * t225 + t223 * t253;
t145 = Icges(4,5) * t225 + t223 * t256;
t298 = t143 * t220 + t145 * t219 + t306 * t223 - t305;
t218 = pkin(10) + qJ(4);
t205 = qJ(5) + t218;
t202 = cos(t205);
t201 = sin(t205);
t285 = Icges(6,4) * t201;
t251 = Icges(6,2) * t202 + t285;
t119 = Icges(6,6) * t225 + t223 * t251;
t120 = Icges(6,6) * t223 - t225 * t251;
t284 = Icges(6,4) * t202;
t254 = Icges(6,1) * t201 + t284;
t121 = Icges(6,5) * t225 + t223 * t254;
t122 = Icges(6,5) * t223 - t225 * t254;
t156 = -Icges(6,2) * t201 + t284;
t157 = Icges(6,1) * t202 - t285;
t197 = qJD(4) * t223 + V_base(5);
t167 = qJD(5) * t223 + t197;
t198 = qJD(4) * t225 + V_base(4);
t168 = qJD(5) * t225 + t198;
t206 = V_base(6) + qJD(1);
t297 = (t119 * t202 + t121 * t201) * t168 + (t120 * t202 + t122 * t201) * t167 + (t156 * t202 + t157 * t201) * t206;
t204 = cos(t218);
t203 = sin(t218);
t287 = Icges(5,4) * t203;
t252 = Icges(5,2) * t204 + t287;
t129 = Icges(5,6) * t225 + t223 * t252;
t130 = Icges(5,6) * t223 - t225 * t252;
t286 = Icges(5,4) * t204;
t255 = Icges(5,1) * t203 + t286;
t131 = Icges(5,5) * t225 + t223 * t255;
t132 = Icges(5,5) * t223 - t225 * t255;
t163 = -Icges(5,2) * t203 + t286;
t164 = Icges(5,1) * t204 - t287;
t296 = (t129 * t204 + t131 * t203) * t198 + (t130 * t204 + t132 * t203) * t197 + (t163 * t204 + t164 * t203) * t206;
t295 = -pkin(2) - pkin(6);
t293 = pkin(3) * t219;
t292 = pkin(3) * t220;
t291 = pkin(4) * t204;
t281 = qJ(3) * t225;
t280 = t202 * t223;
t279 = t202 * t225;
t222 = sin(qJ(6));
t278 = t222 * t225;
t277 = t223 * qJ(3);
t276 = t223 * t222;
t224 = cos(qJ(6));
t275 = t223 * t224;
t274 = t224 * t225;
t271 = qJD(6) * t202;
t190 = t223 * pkin(1) - qJ(2) * t225;
t270 = V_base(4) * t190 + V_base(3);
t269 = V_base(5) * pkin(6) + V_base(1);
t266 = pkin(4) * t203;
t265 = V_base(4) * t277 + t270;
t264 = qJD(2) * t223 + t269;
t263 = -t190 - t277;
t193 = pkin(1) * t225 + t223 * qJ(2);
t262 = -t193 - t281;
t261 = pkin(5) * t201 - pkin(9) * t202;
t153 = pkin(7) * t223 - t225 * t293;
t260 = -t153 + t263;
t259 = rSges(4,1) * t219 + rSges(4,2) * t220;
t258 = rSges(5,1) * t203 + rSges(5,2) * t204;
t257 = rSges(6,1) * t201 + rSges(6,2) * t202;
t250 = Icges(4,5) * t219 + Icges(4,6) * t220;
t249 = Icges(5,5) * t203 + Icges(5,6) * t204;
t248 = Icges(6,5) * t201 + Icges(6,6) * t202;
t175 = -Icges(4,2) * t219 + t288;
t176 = Icges(4,1) * t220 - t289;
t239 = t175 * t220 + t176 * t219;
t238 = -qJD(2) * t225 + t206 * t193 + V_base(2);
t237 = V_base(5) * pkin(2) + qJD(3) * t225 + t264;
t110 = pkin(8) * t223 - t225 * t266;
t236 = -t110 + t260;
t235 = V_base(5) * t292 + t237;
t234 = qJD(3) * t223 + t206 * t281 + t238;
t233 = t197 * t291 + t235;
t232 = (Icges(6,3) * t225 + t223 * t248) * t168 + (Icges(6,3) * t223 - t225 * t248) * t167 + (Icges(6,5) * t202 - Icges(6,6) * t201) * t206;
t231 = (Icges(5,3) * t225 + t223 * t249) * t198 + (Icges(5,3) * t223 - t225 * t249) * t197 + (Icges(5,5) * t204 - Icges(5,6) * t203) * t206;
t230 = (Icges(4,3) * t225 + t223 * t250) * V_base(4) + (Icges(4,3) * t223 - t225 * t250) * V_base(5) + (Icges(4,5) * t220 - Icges(4,6) * t219) * t206;
t154 = pkin(7) * t225 + t223 * t293;
t229 = V_base(4) * t153 + (-t154 + t262) * V_base(5) + t265;
t228 = t206 * t154 + (-t292 + t295) * V_base(4) + t234;
t111 = pkin(8) * t225 + t223 * t266;
t227 = t198 * t110 - t197 * t111 + t229;
t226 = t206 * t111 - t198 * t291 + t228;
t195 = rSges(2,1) * t225 - t223 * rSges(2,2);
t194 = -rSges(3,2) * t225 + t223 * rSges(3,3);
t192 = t223 * rSges(2,1) + rSges(2,2) * t225;
t191 = -t223 * rSges(3,2) - rSges(3,3) * t225;
t177 = rSges(4,1) * t220 - rSges(4,2) * t219;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t166 = qJD(6) * t201 + t206;
t165 = rSges(5,1) * t204 - rSges(5,2) * t203;
t159 = pkin(5) * t202 + pkin(9) * t201;
t158 = rSges(6,1) * t202 - rSges(6,2) * t201;
t152 = -t201 * t274 + t276;
t151 = t201 * t278 + t275;
t150 = t201 * t275 + t278;
t149 = -t201 * t276 + t274;
t148 = t223 * rSges(4,3) - t225 * t259;
t147 = rSges(4,3) * t225 + t223 * t259;
t140 = t261 * t225;
t139 = t261 * t223;
t137 = t223 * rSges(5,3) - t225 * t258;
t136 = rSges(5,3) * t225 + t223 * t258;
t135 = -t223 * t271 + t168;
t134 = t225 * t271 + t167;
t126 = t223 * rSges(6,3) - t225 * t257;
t125 = rSges(6,3) * t225 + t223 * t257;
t124 = V_base(5) * rSges(2,3) - t192 * t206 + t269;
t123 = t195 * t206 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t116 = t192 * V_base(4) - t195 * V_base(5) + V_base(3);
t115 = rSges(7,3) * t201 + (rSges(7,1) * t224 - rSges(7,2) * t222) * t202;
t114 = Icges(7,5) * t201 + (Icges(7,1) * t224 - Icges(7,4) * t222) * t202;
t113 = Icges(7,6) * t201 + (Icges(7,4) * t224 - Icges(7,2) * t222) * t202;
t112 = Icges(7,3) * t201 + (Icges(7,5) * t224 - Icges(7,6) * t222) * t202;
t107 = V_base(5) * rSges(3,1) + (-t190 - t191) * t206 + t264;
t106 = t206 * t194 + (-rSges(3,1) - pkin(6)) * V_base(4) + t238;
t105 = t152 * rSges(7,1) + t151 * rSges(7,2) + rSges(7,3) * t279;
t104 = rSges(7,1) * t150 + rSges(7,2) * t149 - rSges(7,3) * t280;
t103 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t279;
t102 = Icges(7,1) * t150 + Icges(7,4) * t149 - Icges(7,5) * t280;
t101 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t279;
t100 = Icges(7,4) * t150 + Icges(7,2) * t149 - Icges(7,6) * t280;
t99 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t279;
t98 = Icges(7,5) * t150 + Icges(7,6) * t149 - Icges(7,3) * t280;
t97 = t191 * V_base(4) + (-t193 - t194) * V_base(5) + t270;
t96 = t177 * V_base(5) + (-t148 + t263) * t206 + t237;
t95 = t206 * t147 + (-t177 + t295) * V_base(4) + t234;
t94 = V_base(4) * t148 + (-t147 + t262) * V_base(5) + t265;
t93 = t165 * t197 + (-t137 + t260) * t206 + t235;
t92 = t206 * t136 - t198 * t165 + t228;
t91 = -t197 * t136 + t198 * t137 + t229;
t90 = t158 * t167 + (-t126 + t236) * t206 + t233;
t89 = t206 * t125 - t168 * t158 + t226;
t88 = -t167 * t125 + t168 * t126 + t227;
t87 = -t105 * t166 + t115 * t134 + t159 * t167 + (t140 + t236) * t206 + t233;
t86 = t166 * t104 - t135 * t115 + t206 * t139 - t168 * t159 + t226;
t85 = -t134 * t104 + t135 * t105 - t167 * t139 - t168 * t140 + t227;
t1 = t166 * ((t112 * t166 + t99 * t134 + t98 * t135) * t201 + ((-t100 * t222 + t102 * t224) * t135 + (-t101 * t222 + t103 * t224) * t134 + (-t113 * t222 + t114 * t224) * t166) * t202) / 0.2e1 + m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t198 * (t296 * t223 + t231 * t225) / 0.2e1 + t197 * (t231 * t223 - t296 * t225) / 0.2e1 + t168 * (t297 * t223 + t232 * t225) / 0.2e1 + t167 * (t232 * t223 - t297 * t225) / 0.2e1 + t134 * ((t151 * t100 + t152 * t102 + t279 * t98) * t135 + (t151 * t101 + t152 * t103 + t99 * t279) * t134 + (t112 * t279 + t151 * t113 + t152 * t114) * t166) / 0.2e1 + t135 * ((t149 * t100 + t150 * t102 - t98 * t280) * t135 + (t101 * t149 + t103 * t150 - t280 * t99) * t134 + (-t112 * t280 + t113 * t149 + t114 * t150) * t166) / 0.2e1 + (t230 * t225 + (t223 * t239 + t300) * t206 + (t299 * t223 + t303 * t225 + Icges(1,4)) * V_base(5) + (t298 * t223 + t302 * t225 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t230 * t223 + (-t225 * t239 + t301) * t206 + (t303 * t223 - t299 * t225 + Icges(1,2)) * V_base(5) + (t302 * t223 - t298 * t225 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t129 * t203 + t131 * t204) * t198 + (-t130 * t203 + t132 * t204) * t197 + (-t119 * t201 + t121 * t202) * t168 + (-t120 * t201 + t122 * t202) * t167 + (-t144 * t219 + t146 * t220 + t301) * V_base(5) + (-t143 * t219 + t145 * t220 + t300) * V_base(4) + (-t201 * t156 + t202 * t157 - t203 * t163 + t204 * t164 - t219 * t175 + t220 * t176 + Icges(3,1) + Icges(2,3)) * t206) * t206 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
