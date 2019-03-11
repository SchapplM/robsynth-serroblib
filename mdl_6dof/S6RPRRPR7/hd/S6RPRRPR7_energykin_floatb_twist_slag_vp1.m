% Calculate kinetic energy for
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:47
% EndTime: 2019-03-09 05:19:50
% DurationCPUTime: 2.80s
% Computational Cost: add. (1486->307), mult. (1526->429), div. (0->0), fcn. (1308->10), ass. (0->161)
t316 = Icges(2,4) + Icges(3,6);
t315 = Icges(2,1) + Icges(3,2);
t314 = -Icges(3,4) + Icges(2,5);
t313 = Icges(3,5) - Icges(2,6);
t312 = Icges(2,2) + Icges(3,3);
t311 = Icges(5,3) + Icges(6,3);
t218 = qJ(3) + qJ(4);
t203 = pkin(10) + t218;
t201 = sin(t203);
t202 = cos(t203);
t212 = sin(t218);
t213 = cos(t218);
t310 = Icges(5,5) * t212 + Icges(6,5) * t201 + Icges(5,6) * t213 + Icges(6,6) * t202;
t224 = cos(qJ(1));
t309 = t316 * t224;
t221 = sin(qJ(1));
t308 = t316 * t221;
t281 = Icges(6,4) * t201;
t251 = Icges(6,2) * t202 + t281;
t119 = Icges(6,6) * t224 + t221 * t251;
t120 = Icges(6,6) * t221 - t224 * t251;
t280 = Icges(6,4) * t202;
t254 = Icges(6,1) * t201 + t280;
t121 = Icges(6,5) * t224 + t221 * t254;
t122 = Icges(6,5) * t221 - t224 * t254;
t283 = Icges(5,4) * t212;
t252 = Icges(5,2) * t213 + t283;
t132 = Icges(5,6) * t224 + t221 * t252;
t133 = Icges(5,6) * t221 - t224 * t252;
t282 = Icges(5,4) * t213;
t255 = Icges(5,1) * t212 + t282;
t134 = Icges(5,5) * t224 + t221 * t255;
t135 = Icges(5,5) * t221 - t224 * t255;
t157 = -Icges(6,2) * t201 + t280;
t158 = Icges(6,1) * t202 - t281;
t163 = -Icges(5,2) * t212 + t282;
t164 = Icges(5,1) * t213 - t283;
t198 = qJD(3) * t221 + V_base(5);
t167 = qJD(4) * t221 + t198;
t199 = qJD(3) * t224 + V_base(4);
t168 = qJD(4) * t224 + t199;
t204 = V_base(6) + qJD(1);
t307 = t167 * (t120 * t202 + t122 * t201 + t133 * t213 + t135 * t212) + t168 * (t119 * t202 + t121 * t201 + t132 * t213 + t134 * t212) + t204 * (t157 * t202 + t158 * t201 + t163 * t213 + t164 * t212);
t306 = -t312 * t224 - t308;
t305 = t312 * t221 - t309;
t304 = t315 * t221 + t309;
t303 = t315 * t224 - t308;
t300 = (Icges(5,5) * t213 + Icges(6,5) * t202 - Icges(5,6) * t212 - Icges(6,6) * t201) * t204 + (t310 * t221 + t311 * t224) * t168 + (t311 * t221 - t310 * t224) * t167;
t223 = cos(qJ(3));
t220 = sin(qJ(3));
t285 = Icges(4,4) * t220;
t253 = Icges(4,2) * t223 + t285;
t143 = Icges(4,6) * t224 + t221 * t253;
t144 = Icges(4,6) * t221 - t224 * t253;
t284 = Icges(4,4) * t223;
t256 = Icges(4,1) * t220 + t284;
t145 = Icges(4,5) * t224 + t221 * t256;
t146 = Icges(4,5) * t221 - t224 * t256;
t182 = -Icges(4,2) * t220 + t284;
t187 = Icges(4,1) * t223 - t285;
t296 = (t143 * t223 + t145 * t220) * t199 + (t144 * t223 + t146 * t220) * t198 + (t182 * t223 + t187 * t220) * t204;
t292 = pkin(3) * t220;
t291 = pkin(3) * t223;
t290 = pkin(4) * t213;
t289 = t221 * pkin(7);
t288 = t224 * pkin(7);
t277 = t202 * t221;
t276 = t202 * t224;
t219 = sin(qJ(6));
t275 = t219 * t221;
t274 = t219 * t224;
t222 = cos(qJ(6));
t273 = t221 * t222;
t272 = t222 * t224;
t270 = qJD(6) * t202;
t190 = pkin(1) * t221 - qJ(2) * t224;
t269 = V_base(4) * t190 + V_base(3);
t268 = V_base(5) * pkin(6) + V_base(1);
t265 = pkin(4) * t212;
t264 = -t190 - t289;
t263 = qJD(2) * t221 + t268;
t154 = pkin(8) * t221 - t224 * t292;
t262 = -t154 + t264;
t261 = V_base(5) * pkin(2) + t263;
t260 = pkin(5) * t201 - pkin(9) * t202;
t259 = rSges(4,1) * t220 + rSges(4,2) * t223;
t258 = rSges(5,1) * t212 + rSges(5,2) * t213;
t257 = rSges(6,1) * t201 + rSges(6,2) * t202;
t250 = Icges(4,5) * t220 + Icges(4,6) * t223;
t194 = pkin(1) * t224 + qJ(2) * t221;
t238 = -qJD(2) * t224 + t204 * t194 + V_base(2);
t110 = qJ(5) * t221 - t224 * t265;
t237 = -t110 + t262;
t236 = t198 * t291 + t261;
t235 = qJD(5) * t224 + t167 * t290 + t236;
t232 = (Icges(4,3) * t224 + t221 * t250) * t199 + (Icges(4,3) * t221 - t224 * t250) * t198 + (Icges(4,5) * t223 - Icges(4,6) * t220) * t204;
t231 = V_base(4) * t289 + (-t194 - t288) * V_base(5) + t269;
t230 = t204 * t288 + (-pkin(2) - pkin(6)) * V_base(4) + t238;
t155 = pkin(8) * t224 + t221 * t292;
t229 = t199 * t154 - t155 * t198 + t231;
t228 = t168 * t110 + t229;
t227 = t204 * t155 - t199 * t291 + t230;
t111 = qJ(5) * t224 + t221 * t265;
t226 = qJD(5) * t221 + t204 * t111 + t227;
t196 = rSges(2,1) * t224 - rSges(2,2) * t221;
t195 = -rSges(3,2) * t224 + rSges(3,3) * t221;
t193 = rSges(4,1) * t223 - rSges(4,2) * t220;
t192 = rSges(2,1) * t221 + rSges(2,2) * t224;
t191 = -rSges(3,2) * t221 - rSges(3,3) * t224;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t166 = qJD(6) * t201 + t204;
t165 = rSges(5,1) * t213 - rSges(5,2) * t212;
t160 = pkin(5) * t202 + pkin(9) * t201;
t159 = rSges(6,1) * t202 - rSges(6,2) * t201;
t153 = -t201 * t272 + t275;
t152 = t201 * t274 + t273;
t151 = t201 * t273 + t274;
t150 = -t201 * t275 + t272;
t148 = rSges(4,3) * t221 - t224 * t259;
t147 = rSges(4,3) * t224 + t221 * t259;
t140 = t260 * t224;
t139 = t260 * t221;
t137 = rSges(5,3) * t221 - t224 * t258;
t136 = rSges(5,3) * t224 + t221 * t258;
t129 = -t221 * t270 + t168;
t128 = t224 * t270 + t167;
t127 = rSges(6,3) * t221 - t224 * t257;
t126 = rSges(6,3) * t224 + t221 * t257;
t125 = V_base(5) * rSges(2,3) - t192 * t204 + t268;
t124 = t196 * t204 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t116 = t192 * V_base(4) - t196 * V_base(5) + V_base(3);
t115 = rSges(7,3) * t201 + (rSges(7,1) * t222 - rSges(7,2) * t219) * t202;
t114 = Icges(7,5) * t201 + (Icges(7,1) * t222 - Icges(7,4) * t219) * t202;
t113 = Icges(7,6) * t201 + (Icges(7,4) * t222 - Icges(7,2) * t219) * t202;
t112 = Icges(7,3) * t201 + (Icges(7,5) * t222 - Icges(7,6) * t219) * t202;
t108 = V_base(5) * rSges(3,1) + (-t190 - t191) * t204 + t263;
t107 = t195 * t204 + (-rSges(3,1) - pkin(6)) * V_base(4) + t238;
t105 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t276;
t104 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t277;
t103 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t276;
t102 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t277;
t101 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t276;
t100 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t277;
t99 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t276;
t98 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t277;
t97 = t191 * V_base(4) + (-t194 - t195) * V_base(5) + t269;
t96 = t193 * t198 + (-t148 + t264) * t204 + t261;
t95 = t147 * t204 - t193 * t199 + t230;
t94 = -t147 * t198 + t148 * t199 + t231;
t93 = t165 * t167 + (-t137 + t262) * t204 + t236;
t92 = t136 * t204 - t165 * t168 + t227;
t91 = -t136 * t167 + t137 * t168 + t229;
t90 = t159 * t167 + (-t127 + t237) * t204 + t235;
t89 = t126 * t204 + (-t159 - t290) * t168 + t226;
t88 = t127 * t168 + (-t111 - t126) * t167 + t228;
t87 = -t105 * t166 + t115 * t128 + t160 * t167 + (t140 + t237) * t204 + t235;
t86 = t104 * t166 - t115 * t129 + t139 * t204 + (-t160 - t290) * t168 + t226;
t85 = -t104 * t128 + t105 * t129 - t140 * t168 + (-t111 - t139) * t167 + t228;
t1 = m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t107 ^ 2 + t108 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t129 * ((t150 * t100 + t151 * t102 - t98 * t277) * t129 + (t101 * t150 + t103 * t151 - t277 * t99) * t128 + (-t112 * t277 + t113 * t150 + t114 * t151) * t166) / 0.2e1 + t128 * ((t100 * t152 + t102 * t153 + t276 * t98) * t129 + (t152 * t101 + t153 * t103 + t99 * t276) * t128 + (t112 * t276 + t113 * t152 + t114 * t153) * t166) / 0.2e1 + t199 * (t296 * t221 + t232 * t224) / 0.2e1 + t198 * (t232 * t221 - t296 * t224) / 0.2e1 + t166 * ((t112 * t166 + t99 * t128 + t98 * t129) * t201 + ((-t100 * t219 + t102 * t222) * t129 + (-t101 * t219 + t103 * t222) * t128 + (-t113 * t219 + t114 * t222) * t166) * t202) / 0.2e1 + (t300 * t221 - t307 * t224) * t167 / 0.2e1 + (t307 * t221 + t300 * t224) * t168 / 0.2e1 + ((t221 * t306 + t304 * t224 + Icges(1,4)) * V_base(5) + (t305 * t221 + t303 * t224 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t304 * t221 - t306 * t224 + Icges(1,2)) * V_base(5) + (t221 * t303 - t224 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t143 * t220 + t145 * t223) * t199 + (-t144 * t220 + t146 * t223) * t198 + (-t119 * t201 + t121 * t202 - t132 * t212 + t134 * t213) * t168 + (-t120 * t201 + t122 * t202 - t133 * t212 + t135 * t213) * t167 + (-t201 * t157 + t202 * t158 - t212 * t163 + t213 * t164 - t220 * t182 + t223 * t187 + Icges(3,1) + Icges(2,3)) * t204) * t204 / 0.2e1 + t204 * V_base(5) * (t314 * t221 - t313 * t224) + t204 * V_base(4) * (t313 * t221 + t314 * t224) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
