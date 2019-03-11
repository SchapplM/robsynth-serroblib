% Calculate kinetic energy for
% S6RPRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:54:53
% EndTime: 2019-03-09 03:54:55
% DurationCPUTime: 2.74s
% Computational Cost: add. (1468->307), mult. (1508->429), div. (0->0), fcn. (1290->10), ass. (0->161)
t316 = Icges(2,4) + Icges(3,6);
t315 = Icges(2,1) + Icges(3,2);
t314 = -Icges(3,4) + Icges(2,5);
t313 = Icges(3,5) - Icges(2,6);
t312 = Icges(2,2) + Icges(3,3);
t311 = Icges(4,3) + Icges(5,3);
t218 = qJ(3) + pkin(10);
t203 = sin(t218);
t204 = cos(t218);
t221 = sin(qJ(3));
t224 = cos(qJ(3));
t310 = Icges(4,5) * t221 + Icges(5,5) * t203 + Icges(4,6) * t224 + Icges(5,6) * t204;
t225 = cos(qJ(1));
t309 = t316 * t225;
t222 = sin(qJ(1));
t308 = t316 * t222;
t283 = Icges(5,4) * t203;
t252 = Icges(5,2) * t204 + t283;
t130 = Icges(5,6) * t225 + t222 * t252;
t131 = Icges(5,6) * t222 - t225 * t252;
t282 = Icges(5,4) * t204;
t255 = Icges(5,1) * t203 + t282;
t132 = Icges(5,5) * t225 + t222 * t255;
t133 = Icges(5,5) * t222 - t225 * t255;
t285 = Icges(4,4) * t221;
t253 = Icges(4,2) * t224 + t285;
t143 = Icges(4,6) * t225 + t222 * t253;
t144 = Icges(4,6) * t222 - t225 * t253;
t284 = Icges(4,4) * t224;
t256 = Icges(4,1) * t221 + t284;
t145 = Icges(4,5) * t225 + t222 * t256;
t146 = Icges(4,5) * t222 - t225 * t256;
t163 = -Icges(5,2) * t203 + t282;
t164 = Icges(5,1) * t204 - t283;
t182 = -Icges(4,2) * t221 + t284;
t187 = Icges(4,1) * t224 - t285;
t198 = qJD(3) * t222 + V_base(5);
t199 = qJD(3) * t225 + V_base(4);
t206 = V_base(6) + qJD(1);
t307 = t198 * (t131 * t204 + t133 * t203 + t144 * t224 + t146 * t221) + t199 * (t130 * t204 + t132 * t203 + t143 * t224 + t145 * t221) + t206 * (t163 * t204 + t164 * t203 + t182 * t224 + t187 * t221);
t306 = -t312 * t225 - t308;
t305 = t312 * t222 - t309;
t304 = t315 * t222 + t309;
t303 = t315 * t225 - t308;
t300 = (Icges(4,5) * t224 + Icges(5,5) * t204 - Icges(4,6) * t221 - Icges(5,6) * t203) * t206 + (t310 * t222 + t311 * t225) * t199 + (t311 * t222 - t310 * t225) * t198;
t205 = qJ(5) + t218;
t202 = cos(t205);
t201 = sin(t205);
t281 = Icges(6,4) * t201;
t251 = Icges(6,2) * t202 + t281;
t119 = Icges(6,6) * t225 + t222 * t251;
t120 = Icges(6,6) * t222 - t225 * t251;
t280 = Icges(6,4) * t202;
t254 = Icges(6,1) * t201 + t280;
t121 = Icges(6,5) * t225 + t222 * t254;
t122 = Icges(6,5) * t222 - t225 * t254;
t156 = -Icges(6,2) * t201 + t280;
t157 = Icges(6,1) * t202 - t281;
t167 = qJD(5) * t222 + t198;
t168 = qJD(5) * t225 + t199;
t296 = (t119 * t202 + t121 * t201) * t168 + (t120 * t202 + t122 * t201) * t167 + (t156 * t202 + t157 * t201) * t206;
t292 = pkin(3) * t221;
t291 = pkin(3) * t224;
t290 = pkin(4) * t204;
t289 = pkin(7) * t225;
t288 = t222 * pkin(7);
t277 = t202 * t222;
t276 = t202 * t225;
t220 = sin(qJ(6));
t275 = t220 * t225;
t274 = t222 * t220;
t223 = cos(qJ(6));
t273 = t222 * t223;
t272 = t223 * t225;
t270 = qJD(6) * t202;
t190 = t222 * pkin(1) - qJ(2) * t225;
t269 = V_base(4) * t190 + V_base(3);
t268 = V_base(5) * pkin(6) + V_base(1);
t265 = pkin(4) * t203;
t264 = -t190 - t288;
t263 = qJD(2) * t222 + t268;
t153 = qJ(4) * t222 - t225 * t292;
t262 = -t153 + t264;
t261 = V_base(5) * pkin(2) + t263;
t260 = pkin(5) * t201 - pkin(9) * t202;
t259 = rSges(4,1) * t221 + rSges(4,2) * t224;
t258 = rSges(5,1) * t203 + rSges(5,2) * t204;
t257 = rSges(6,1) * t201 + rSges(6,2) * t202;
t248 = Icges(6,5) * t201 + Icges(6,6) * t202;
t194 = pkin(1) * t225 + t222 * qJ(2);
t238 = -qJD(2) * t225 + t206 * t194 + V_base(2);
t110 = pkin(8) * t222 - t225 * t265;
t237 = -t110 + t262;
t236 = qJD(4) * t225 + t198 * t291 + t261;
t235 = t198 * t290 + t236;
t234 = (Icges(6,3) * t225 + t222 * t248) * t168 + (Icges(6,3) * t222 - t225 * t248) * t167 + (Icges(6,5) * t202 - Icges(6,6) * t201) * t206;
t231 = V_base(4) * t288 + (-t194 - t289) * V_base(5) + t269;
t230 = t206 * t289 + (-pkin(2) - pkin(6)) * V_base(4) + t238;
t229 = t199 * t153 + t231;
t154 = qJ(4) * t225 + t222 * t292;
t228 = qJD(4) * t222 + t206 * t154 + t230;
t111 = pkin(8) * t225 + t222 * t265;
t227 = t199 * t110 + (-t111 - t154) * t198 + t229;
t226 = t206 * t111 + (-t290 - t291) * t199 + t228;
t196 = rSges(2,1) * t225 - t222 * rSges(2,2);
t195 = -rSges(3,2) * t225 + t222 * rSges(3,3);
t193 = rSges(4,1) * t224 - rSges(4,2) * t221;
t192 = t222 * rSges(2,1) + rSges(2,2) * t225;
t191 = -t222 * rSges(3,2) - rSges(3,3) * t225;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t166 = qJD(6) * t201 + t206;
t165 = rSges(5,1) * t204 - rSges(5,2) * t203;
t159 = pkin(5) * t202 + pkin(9) * t201;
t158 = rSges(6,1) * t202 - rSges(6,2) * t201;
t152 = -t201 * t272 + t274;
t151 = t201 * t275 + t273;
t150 = t201 * t273 + t275;
t149 = -t201 * t274 + t272;
t148 = t222 * rSges(4,3) - t225 * t259;
t147 = rSges(4,3) * t225 + t222 * t259;
t140 = t260 * t225;
t139 = t260 * t222;
t137 = t222 * rSges(5,3) - t225 * t258;
t136 = rSges(5,3) * t225 + t222 * t258;
t135 = -t222 * t270 + t168;
t134 = t225 * t270 + t167;
t127 = t222 * rSges(6,3) - t225 * t257;
t126 = rSges(6,3) * t225 + t222 * t257;
t125 = V_base(5) * rSges(2,3) - t192 * t206 + t268;
t124 = t196 * t206 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t116 = t192 * V_base(4) - t196 * V_base(5) + V_base(3);
t115 = rSges(7,3) * t201 + (rSges(7,1) * t223 - rSges(7,2) * t220) * t202;
t114 = Icges(7,5) * t201 + (Icges(7,1) * t223 - Icges(7,4) * t220) * t202;
t113 = Icges(7,6) * t201 + (Icges(7,4) * t223 - Icges(7,2) * t220) * t202;
t112 = Icges(7,3) * t201 + (Icges(7,5) * t223 - Icges(7,6) * t220) * t202;
t107 = V_base(5) * rSges(3,1) + (-t190 - t191) * t206 + t263;
t106 = t206 * t195 + (-rSges(3,1) - pkin(6)) * V_base(4) + t238;
t105 = t152 * rSges(7,1) + t151 * rSges(7,2) + rSges(7,3) * t276;
t104 = rSges(7,1) * t150 + rSges(7,2) * t149 - rSges(7,3) * t277;
t103 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t276;
t102 = Icges(7,1) * t150 + Icges(7,4) * t149 - Icges(7,5) * t277;
t101 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t276;
t100 = Icges(7,4) * t150 + Icges(7,2) * t149 - Icges(7,6) * t277;
t99 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t276;
t98 = Icges(7,5) * t150 + Icges(7,6) * t149 - Icges(7,3) * t277;
t97 = t191 * V_base(4) + (-t194 - t195) * V_base(5) + t269;
t96 = t193 * t198 + (-t148 + t264) * t206 + t261;
t95 = t206 * t147 - t199 * t193 + t230;
t94 = -t198 * t147 + t199 * t148 + t231;
t93 = t165 * t198 + (-t137 + t262) * t206 + t236;
t92 = t206 * t136 + (-t165 - t291) * t199 + t228;
t91 = t199 * t137 + (-t136 - t154) * t198 + t229;
t90 = t158 * t167 + (-t127 + t237) * t206 + t235;
t89 = t206 * t126 - t168 * t158 + t226;
t88 = -t167 * t126 + t168 * t127 + t227;
t87 = -t105 * t166 + t115 * t134 + t159 * t167 + (t140 + t237) * t206 + t235;
t86 = t166 * t104 - t135 * t115 + t206 * t139 - t168 * t159 + t226;
t85 = -t134 * t104 + t135 * t105 - t167 * t139 - t168 * t140 + t227;
t1 = t166 * ((t112 * t166 + t99 * t134 + t98 * t135) * t201 + ((-t100 * t220 + t102 * t223) * t135 + (-t101 * t220 + t103 * t223) * t134 + (-t113 * t220 + t114 * t223) * t166) * t202) / 0.2e1 + m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t134 * ((t151 * t100 + t152 * t102 + t98 * t276) * t135 + (t151 * t101 + t152 * t103 + t99 * t276) * t134 + (t112 * t276 + t151 * t113 + t152 * t114) * t166) / 0.2e1 + t135 * ((t149 * t100 + t150 * t102 - t98 * t277) * t135 + (t101 * t149 + t103 * t150 - t99 * t277) * t134 + (-t112 * t277 + t113 * t149 + t114 * t150) * t166) / 0.2e1 + t168 * (t296 * t222 + t234 * t225) / 0.2e1 + t167 * (t234 * t222 - t296 * t225) / 0.2e1 + (t300 * t222 - t307 * t225) * t198 / 0.2e1 + (t307 * t222 + t300 * t225) * t199 / 0.2e1 + ((t222 * t306 + t304 * t225 + Icges(1,4)) * V_base(5) + (t305 * t222 + t303 * t225 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t304 * t222 - t306 * t225 + Icges(1,2)) * V_base(5) + (t222 * t303 - t225 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t119 * t201 + t121 * t202) * t168 + (-t120 * t201 + t122 * t202) * t167 + (-t130 * t203 + t132 * t204 - t143 * t221 + t145 * t224) * t199 + (-t131 * t203 + t133 * t204 - t144 * t221 + t146 * t224) * t198 + (-t201 * t156 + t202 * t157 - t203 * t163 + t204 * t164 - t221 * t182 + t224 * t187 + Icges(3,1) + Icges(2,3)) * t206) * t206 / 0.2e1 + t206 * V_base(5) * (t314 * t222 - t313 * t225) + t206 * V_base(4) * (t313 * t222 + t314 * t225) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
