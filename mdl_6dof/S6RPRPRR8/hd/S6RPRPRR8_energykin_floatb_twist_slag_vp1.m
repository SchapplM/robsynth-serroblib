% Calculate kinetic energy for
% S6RPRPRR8
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:40
% EndTime: 2019-03-09 03:57:43
% DurationCPUTime: 3.33s
% Computational Cost: add. (1561->332), mult. (1758->468), div. (0->0), fcn. (1596->10), ass. (0->165)
t314 = Icges(2,4) + Icges(3,6);
t313 = Icges(2,1) + Icges(3,2);
t312 = -Icges(3,4) + Icges(2,5);
t311 = Icges(3,5) - Icges(2,6);
t310 = Icges(2,2) + Icges(3,3);
t309 = Icges(4,3) + Icges(5,3);
t222 = qJ(3) + pkin(10);
t209 = sin(t222);
t210 = cos(t222);
t226 = sin(qJ(3));
t229 = cos(qJ(3));
t308 = Icges(4,5) * t226 + Icges(5,5) * t209 + Icges(4,6) * t229 + Icges(5,6) * t210;
t230 = cos(qJ(1));
t307 = t314 * t230;
t227 = sin(qJ(1));
t306 = t314 * t227;
t280 = Icges(5,4) * t209;
t251 = Icges(5,2) * t210 + t280;
t137 = Icges(5,6) * t230 + t227 * t251;
t138 = Icges(5,6) * t227 - t230 * t251;
t279 = Icges(5,4) * t210;
t253 = Icges(5,1) * t209 + t279;
t139 = Icges(5,5) * t230 + t227 * t253;
t140 = Icges(5,5) * t227 - t230 * t253;
t282 = Icges(4,4) * t226;
t252 = Icges(4,2) * t229 + t282;
t151 = Icges(4,6) * t230 + t227 * t252;
t152 = Icges(4,6) * t227 - t230 * t252;
t281 = Icges(4,4) * t229;
t254 = Icges(4,1) * t226 + t281;
t153 = Icges(4,5) * t230 + t227 * t254;
t154 = Icges(4,5) * t227 - t230 * t254;
t170 = -Icges(5,2) * t209 + t279;
t171 = Icges(5,1) * t210 - t280;
t187 = -Icges(4,2) * t226 + t281;
t192 = Icges(4,1) * t229 - t282;
t204 = qJD(3) * t227 + V_base(5);
t205 = qJD(3) * t230 + V_base(4);
t211 = V_base(6) + qJD(1);
t305 = t204 * (t138 * t210 + t140 * t209 + t152 * t229 + t154 * t226) + t205 * (t137 * t210 + t139 * t209 + t151 * t229 + t153 * t226) + t211 * (t170 * t210 + t171 * t209 + t187 * t229 + t192 * t226);
t304 = -t230 * t310 - t306;
t303 = t227 * t310 - t307;
t302 = t227 * t313 + t307;
t301 = t230 * t313 - t306;
t298 = (Icges(4,5) * t229 + Icges(5,5) * t210 - Icges(4,6) * t226 - Icges(5,6) * t209) * t211 + (t227 * t308 + t230 * t309) * t205 + (t227 * t309 - t230 * t308) * t204;
t228 = cos(qJ(5));
t286 = pkin(5) * t228;
t294 = -pkin(9) * t210 + t209 * t286;
t290 = pkin(3) * t226;
t289 = pkin(3) * t229;
t288 = pkin(7) * t227;
t287 = pkin(7) * t230;
t276 = t210 * t227;
t275 = t210 * t230;
t223 = qJ(5) + qJ(6);
t217 = sin(t223);
t274 = t217 * t227;
t273 = t217 * t230;
t218 = cos(t223);
t272 = t218 * t227;
t271 = t218 * t230;
t225 = sin(qJ(5));
t270 = t225 * t227;
t269 = t225 * t230;
t268 = t227 * t228;
t267 = t228 * t230;
t266 = qJD(5) * t210;
t196 = pkin(1) * t227 - qJ(2) * t230;
t265 = V_base(4) * t196 + V_base(3);
t264 = V_base(5) * pkin(6) + V_base(1);
t261 = -t196 - t288;
t260 = qJD(2) * t227 + t264;
t162 = t230 * t266 + t204;
t176 = qJD(5) * t209 + t211;
t160 = qJ(4) * t227 - t230 * t290;
t259 = -t160 + t261;
t258 = V_base(5) * pkin(2) + t260;
t257 = pkin(4) * t209 - pkin(8) * t210;
t256 = rSges(4,1) * t226 + rSges(4,2) * t229;
t255 = rSges(5,1) * t209 + rSges(5,2) * t210;
t200 = pkin(1) * t230 + qJ(2) * t227;
t242 = -qJD(2) * t230 + t211 * t200 + V_base(2);
t241 = qJD(4) * t230 + t204 * t289 + t258;
t238 = V_base(4) * t288 + (-t200 - t287) * V_base(5) + t265;
t237 = t211 * t287 + (-pkin(2) - pkin(6)) * V_base(4) + t242;
t236 = t205 * t160 + t238;
t161 = qJ(4) * t230 + t227 * t290;
t235 = qJD(4) * t227 + t211 * t161 + t237;
t158 = t257 * t230;
t173 = t210 * pkin(4) + t209 * pkin(8);
t234 = t204 * t173 + (t158 + t259) * t211 + t241;
t157 = t257 * t227;
t233 = -t205 * t158 + (-t157 - t161) * t204 + t236;
t232 = t211 * t157 + (-t173 - t289) * t205 + t235;
t202 = rSges(2,1) * t230 - rSges(2,2) * t227;
t201 = -rSges(3,2) * t230 + rSges(3,3) * t227;
t199 = rSges(4,1) * t229 - rSges(4,2) * t226;
t198 = rSges(2,1) * t227 + rSges(2,2) * t230;
t197 = -rSges(3,2) * t227 - rSges(3,3) * t230;
t179 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t178 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t177 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = rSges(5,1) * t210 - rSges(5,2) * t209;
t167 = -t209 * t267 + t270;
t166 = t209 * t269 + t268;
t165 = t209 * t268 + t269;
t164 = -t209 * t270 + t267;
t163 = -t227 * t266 + t205;
t159 = qJD(6) * t209 + t176;
t156 = rSges(4,3) * t227 - t230 * t256;
t155 = rSges(4,3) * t230 + t227 * t256;
t148 = -t209 * t271 + t274;
t147 = t209 * t273 + t272;
t146 = t209 * t272 + t273;
t145 = -t209 * t274 + t271;
t142 = rSges(5,3) * t227 - t230 * t255;
t141 = rSges(5,3) * t230 + t227 * t255;
t133 = rSges(6,3) * t209 + (rSges(6,1) * t228 - rSges(6,2) * t225) * t210;
t132 = V_base(5) * rSges(2,3) - t198 * t211 + t264;
t131 = t202 * t211 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t130 = Icges(6,5) * t209 + (Icges(6,1) * t228 - Icges(6,4) * t225) * t210;
t129 = Icges(6,6) * t209 + (Icges(6,4) * t228 - Icges(6,2) * t225) * t210;
t128 = Icges(6,3) * t209 + (Icges(6,5) * t228 - Icges(6,6) * t225) * t210;
t126 = (-qJD(5) - qJD(6)) * t276 + t205;
t125 = qJD(6) * t275 + t162;
t123 = t198 * V_base(4) - t202 * V_base(5) + V_base(3);
t122 = rSges(7,3) * t209 + (rSges(7,1) * t218 - rSges(7,2) * t217) * t210;
t121 = Icges(7,5) * t209 + (Icges(7,1) * t218 - Icges(7,4) * t217) * t210;
t120 = Icges(7,6) * t209 + (Icges(7,4) * t218 - Icges(7,2) * t217) * t210;
t119 = Icges(7,3) * t209 + (Icges(7,5) * t218 - Icges(7,6) * t217) * t210;
t118 = pkin(9) * t209 + t210 * t286;
t117 = V_base(5) * rSges(3,1) + (-t196 - t197) * t211 + t260;
t116 = t201 * t211 + (-rSges(3,1) - pkin(6)) * V_base(4) + t242;
t115 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t275;
t114 = rSges(6,1) * t165 + rSges(6,2) * t164 - rSges(6,3) * t276;
t113 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t275;
t112 = Icges(6,1) * t165 + Icges(6,4) * t164 - Icges(6,5) * t276;
t111 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t275;
t110 = Icges(6,4) * t165 + Icges(6,2) * t164 - Icges(6,6) * t276;
t109 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t275;
t108 = Icges(6,5) * t165 + Icges(6,6) * t164 - Icges(6,3) * t276;
t107 = pkin(5) * t270 - t230 * t294;
t106 = pkin(5) * t269 + t227 * t294;
t105 = t197 * V_base(4) + (-t200 - t201) * V_base(5) + t265;
t104 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t275;
t103 = rSges(7,1) * t146 + rSges(7,2) * t145 - rSges(7,3) * t276;
t102 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t275;
t101 = Icges(7,1) * t146 + Icges(7,4) * t145 - Icges(7,5) * t276;
t100 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t275;
t99 = Icges(7,4) * t146 + Icges(7,2) * t145 - Icges(7,6) * t276;
t98 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t275;
t97 = Icges(7,5) * t146 + Icges(7,6) * t145 - Icges(7,3) * t276;
t96 = t199 * t204 + (-t156 + t261) * t211 + t258;
t95 = t155 * t211 - t199 * t205 + t237;
t94 = -t155 * t204 + t156 * t205 + t238;
t93 = t172 * t204 + (-t142 + t259) * t211 + t241;
t92 = t141 * t211 + (-t172 - t289) * t205 + t235;
t91 = t142 * t205 + (-t141 - t161) * t204 + t236;
t90 = -t115 * t176 + t133 * t162 + t234;
t89 = t114 * t176 - t133 * t163 + t232;
t88 = -t114 * t162 + t115 * t163 + t233;
t87 = -t104 * t159 - t107 * t176 + t118 * t162 + t122 * t125 + t234;
t86 = t103 * t159 + t106 * t176 - t118 * t163 - t122 * t126 + t232;
t85 = -t103 * t125 + t104 * t126 - t106 * t162 + t107 * t163 + t233;
t1 = t163 * ((-t108 * t276 + t110 * t164 + t112 * t165) * t163 + (-t109 * t276 + t111 * t164 + t113 * t165) * t162 + (-t128 * t276 + t129 * t164 + t130 * t165) * t176) / 0.2e1 + t126 * ((t101 * t146 + t145 * t99 - t97 * t276) * t126 + (t100 * t145 + t102 * t146 - t276 * t98) * t125 + (-t119 * t276 + t120 * t145 + t121 * t146) * t159) / 0.2e1 + t162 * ((t108 * t275 + t110 * t166 + t112 * t167) * t163 + (t109 * t275 + t111 * t166 + t113 * t167) * t162 + (t128 * t275 + t129 * t166 + t130 * t167) * t176) / 0.2e1 + t125 * ((t101 * t148 + t147 * t99 + t275 * t97) * t126 + (t100 * t147 + t102 * t148 + t98 * t275) * t125 + (t119 * t275 + t120 * t147 + t121 * t148) * t159) / 0.2e1 + t176 * ((t108 * t163 + t109 * t162 + t128 * t176) * t209 + ((-t110 * t225 + t112 * t228) * t163 + (-t111 * t225 + t113 * t228) * t162 + (-t129 * t225 + t130 * t228) * t176) * t210) / 0.2e1 + t159 * ((t119 * t159 + t98 * t125 + t97 * t126) * t209 + ((t101 * t218 - t217 * t99) * t126 + (-t100 * t217 + t102 * t218) * t125 + (-t120 * t217 + t121 * t218) * t159) * t210) / 0.2e1 + m(1) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + (t298 * t227 - t305 * t230) * t204 / 0.2e1 + (t305 * t227 + t298 * t230) * t205 / 0.2e1 + ((t227 * t304 + t302 * t230 + Icges(1,4)) * V_base(5) + (t227 * t303 + t230 * t301 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t302 * t227 - t230 * t304 + Icges(1,2)) * V_base(5) + (t227 * t301 - t230 * t303 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t137 * t209 + t139 * t210 - t151 * t226 + t153 * t229) * t205 + (-t138 * t209 + t140 * t210 - t152 * t226 + t154 * t229) * t204 + (-t170 * t209 + t171 * t210 - t187 * t226 + t192 * t229 + Icges(3,1) + Icges(2,3)) * t211) * t211 / 0.2e1 + t211 * V_base(5) * (t227 * t312 - t230 * t311) + t211 * V_base(4) * (t227 * t311 + t312 * t230) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
