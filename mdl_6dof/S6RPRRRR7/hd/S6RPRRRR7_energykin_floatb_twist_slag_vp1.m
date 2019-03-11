% Calculate kinetic energy for
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:52
% EndTime: 2019-03-09 07:16:55
% DurationCPUTime: 2.80s
% Computational Cost: add. (1510->309), mult. (1550->451), div. (0->0), fcn. (1332->10), ass. (0->165)
t311 = Icges(2,4) + Icges(3,6);
t310 = Icges(2,1) + Icges(3,2);
t309 = -Icges(3,4) + Icges(2,5);
t308 = Icges(3,5) - Icges(2,6);
t307 = Icges(2,2) + Icges(3,3);
t226 = cos(qJ(1));
t306 = t311 * t226;
t223 = sin(qJ(1));
t305 = t311 * t223;
t304 = -t307 * t226 - t305;
t303 = t307 * t223 - t306;
t302 = t310 * t223 + t306;
t301 = t310 * t226 - t305;
t220 = qJ(3) + qJ(4);
t216 = qJ(5) + t220;
t204 = cos(t216);
t203 = sin(t216);
t283 = Icges(6,4) * t203;
t253 = Icges(6,2) * t204 + t283;
t124 = Icges(6,6) * t226 + t223 * t253;
t125 = Icges(6,6) * t223 - t226 * t253;
t282 = Icges(6,4) * t204;
t256 = Icges(6,1) * t203 + t282;
t126 = Icges(6,5) * t226 + t223 * t256;
t127 = Icges(6,5) * t223 - t226 * t256;
t200 = qJD(3) * t223 + V_base(5);
t168 = qJD(4) * t223 + t200;
t154 = qJD(5) * t223 + t168;
t201 = qJD(3) * t226 + V_base(4);
t169 = qJD(4) * t226 + t201;
t155 = qJD(5) * t226 + t169;
t159 = -Icges(6,2) * t203 + t282;
t160 = Icges(6,1) * t204 - t283;
t205 = V_base(6) + qJD(1);
t298 = (t124 * t204 + t126 * t203) * t155 + (t125 * t204 + t127 * t203) * t154 + (t159 * t204 + t160 * t203) * t205;
t214 = cos(t220);
t213 = sin(t220);
t285 = Icges(5,4) * t213;
t254 = Icges(5,2) * t214 + t285;
t132 = Icges(5,6) * t226 + t223 * t254;
t133 = Icges(5,6) * t223 - t226 * t254;
t284 = Icges(5,4) * t214;
t257 = Icges(5,1) * t213 + t284;
t134 = Icges(5,5) * t226 + t223 * t257;
t135 = Icges(5,5) * t223 - t226 * t257;
t165 = -Icges(5,2) * t213 + t284;
t166 = Icges(5,1) * t214 - t285;
t297 = (t132 * t214 + t134 * t213) * t169 + (t133 * t214 + t135 * t213) * t168 + (t165 * t214 + t166 * t213) * t205;
t225 = cos(qJ(3));
t222 = sin(qJ(3));
t287 = Icges(4,4) * t222;
t255 = Icges(4,2) * t225 + t287;
t143 = Icges(4,6) * t226 + t223 * t255;
t144 = Icges(4,6) * t223 - t226 * t255;
t286 = Icges(4,4) * t225;
t258 = Icges(4,1) * t222 + t286;
t145 = Icges(4,5) * t226 + t223 * t258;
t146 = Icges(4,5) * t223 - t226 * t258;
t184 = -Icges(4,2) * t222 + t286;
t189 = Icges(4,1) * t225 - t287;
t296 = (t143 * t225 + t145 * t222) * t201 + (t144 * t225 + t146 * t222) * t200 + (t184 * t225 + t189 * t222) * t205;
t294 = pkin(3) * t222;
t293 = pkin(3) * t225;
t292 = pkin(4) * t214;
t291 = t223 * pkin(7);
t290 = t226 * pkin(7);
t279 = t204 * t223;
t278 = t204 * t226;
t221 = sin(qJ(6));
t277 = t221 * t223;
t276 = t221 * t226;
t224 = cos(qJ(6));
t275 = t223 * t224;
t274 = t224 * t226;
t272 = qJD(6) * t204;
t192 = pkin(1) * t223 - qJ(2) * t226;
t271 = V_base(4) * t192 + V_base(3);
t270 = V_base(5) * pkin(6) + V_base(1);
t267 = pkin(4) * t213;
t266 = -t192 - t291;
t265 = qJD(2) * t223 + t270;
t156 = pkin(8) * t223 - t226 * t294;
t264 = -t156 + t266;
t263 = V_base(5) * pkin(2) + t265;
t262 = pkin(5) * t203 - pkin(10) * t204;
t261 = rSges(4,1) * t222 + rSges(4,2) * t225;
t260 = rSges(5,1) * t213 + rSges(5,2) * t214;
t259 = rSges(6,1) * t203 + rSges(6,2) * t204;
t252 = Icges(4,5) * t222 + Icges(4,6) * t225;
t251 = Icges(5,5) * t213 + Icges(5,6) * t214;
t250 = Icges(6,5) * t203 + Icges(6,6) * t204;
t196 = pkin(1) * t226 + qJ(2) * t223;
t240 = -qJD(2) * t226 + t205 * t196 + V_base(2);
t110 = pkin(9) * t223 - t226 * t267;
t239 = -t110 + t264;
t238 = t200 * t293 + t263;
t237 = t168 * t292 + t238;
t236 = (Icges(6,3) * t226 + t223 * t250) * t155 + (Icges(6,3) * t223 - t226 * t250) * t154 + (Icges(6,5) * t204 - Icges(6,6) * t203) * t205;
t235 = (Icges(5,3) * t226 + t223 * t251) * t169 + (Icges(5,3) * t223 - t226 * t251) * t168 + (Icges(5,5) * t214 - Icges(5,6) * t213) * t205;
t234 = (Icges(4,3) * t226 + t223 * t252) * t201 + (Icges(4,3) * t223 - t226 * t252) * t200 + (Icges(4,5) * t225 - Icges(4,6) * t222) * t205;
t233 = V_base(4) * t291 + (-t196 - t290) * V_base(5) + t271;
t232 = t205 * t290 + (-pkin(2) - pkin(6)) * V_base(4) + t240;
t157 = pkin(8) * t226 + t223 * t294;
t231 = t201 * t156 - t157 * t200 + t233;
t230 = t205 * t157 - t201 * t293 + t232;
t111 = pkin(9) * t226 + t223 * t267;
t229 = t169 * t110 - t111 * t168 + t231;
t228 = t205 * t111 - t169 * t292 + t230;
t198 = rSges(2,1) * t226 - rSges(2,2) * t223;
t197 = -rSges(3,2) * t226 + rSges(3,3) * t223;
t195 = rSges(4,1) * t225 - rSges(4,2) * t222;
t194 = rSges(2,1) * t223 + rSges(2,2) * t226;
t193 = -rSges(3,2) * t223 - rSges(3,3) * t226;
t175 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t174 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t173 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t171 = qJD(6) * t203 + t205;
t167 = rSges(5,1) * t214 - rSges(5,2) * t213;
t162 = pkin(5) * t204 + pkin(10) * t203;
t161 = rSges(6,1) * t204 - rSges(6,2) * t203;
t153 = -t203 * t274 + t277;
t152 = t203 * t276 + t275;
t151 = t203 * t275 + t276;
t150 = -t203 * t277 + t274;
t148 = rSges(4,3) * t223 - t226 * t261;
t147 = rSges(4,3) * t226 + t223 * t261;
t140 = t262 * t226;
t139 = t262 * t223;
t137 = rSges(5,3) * t223 - t226 * t260;
t136 = rSges(5,3) * t226 + t223 * t260;
t129 = rSges(6,3) * t223 - t226 * t259;
t128 = rSges(6,3) * t226 + t223 * t259;
t121 = V_base(5) * rSges(2,3) - t194 * t205 + t270;
t120 = t198 * t205 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t118 = -t223 * t272 + t155;
t117 = t226 * t272 + t154;
t116 = t194 * V_base(4) - t198 * V_base(5) + V_base(3);
t115 = rSges(7,3) * t203 + (rSges(7,1) * t224 - rSges(7,2) * t221) * t204;
t114 = Icges(7,5) * t203 + (Icges(7,1) * t224 - Icges(7,4) * t221) * t204;
t113 = Icges(7,6) * t203 + (Icges(7,4) * t224 - Icges(7,2) * t221) * t204;
t112 = Icges(7,3) * t203 + (Icges(7,5) * t224 - Icges(7,6) * t221) * t204;
t108 = V_base(5) * rSges(3,1) + (-t192 - t193) * t205 + t265;
t107 = t197 * t205 + (-rSges(3,1) - pkin(6)) * V_base(4) + t240;
t105 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t278;
t104 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t279;
t103 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t278;
t102 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t279;
t101 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t278;
t100 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t279;
t99 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t278;
t98 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t279;
t97 = t193 * V_base(4) + (-t196 - t197) * V_base(5) + t271;
t96 = t195 * t200 + (-t148 + t266) * t205 + t263;
t95 = t147 * t205 - t195 * t201 + t232;
t94 = -t147 * t200 + t148 * t201 + t233;
t93 = t167 * t168 + (-t137 + t264) * t205 + t238;
t92 = t136 * t205 - t167 * t169 + t230;
t91 = -t136 * t168 + t137 * t169 + t231;
t90 = t154 * t161 + (-t129 + t239) * t205 + t237;
t89 = t128 * t205 - t155 * t161 + t228;
t88 = -t128 * t154 + t129 * t155 + t229;
t87 = -t105 * t171 + t115 * t117 + t154 * t162 + (t140 + t239) * t205 + t237;
t86 = t104 * t171 - t115 * t118 + t139 * t205 - t155 * t162 + t228;
t85 = -t104 * t117 + t105 * t118 - t139 * t154 - t140 * t155 + t229;
t1 = m(1) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t107 ^ 2 + t108 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + t118 * ((t150 * t100 + t151 * t102 - t98 * t279) * t118 + (t101 * t150 + t103 * t151 - t99 * t279) * t117 + (-t112 * t279 + t113 * t150 + t114 * t151) * t171) / 0.2e1 + t117 * ((t100 * t152 + t102 * t153 + t278 * t98) * t118 + (t152 * t101 + t153 * t103 + t99 * t278) * t117 + (t112 * t278 + t113 * t152 + t114 * t153) * t171) / 0.2e1 + t200 * (t234 * t223 - t296 * t226) / 0.2e1 + t169 * (t297 * t223 + t235 * t226) / 0.2e1 + t168 * (t235 * t223 - t297 * t226) / 0.2e1 + t155 * (t298 * t223 + t236 * t226) / 0.2e1 + t154 * (t236 * t223 - t298 * t226) / 0.2e1 + t201 * (t296 * t223 + t234 * t226) / 0.2e1 + t171 * ((t112 * t171 + t99 * t117 + t98 * t118) * t203 + ((-t100 * t221 + t102 * t224) * t118 + (-t101 * t221 + t103 * t224) * t117 + (-t113 * t221 + t114 * t224) * t171) * t204) / 0.2e1 + ((t223 * t304 + t302 * t226 + Icges(1,4)) * V_base(5) + (t303 * t223 + t301 * t226 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t302 * t223 - t304 * t226 + Icges(1,2)) * V_base(5) + (t223 * t301 - t226 * t303 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t143 * t222 + t145 * t225) * t201 + (-t144 * t222 + t146 * t225) * t200 + (-t132 * t213 + t134 * t214) * t169 + (-t133 * t213 + t135 * t214) * t168 + (-t124 * t203 + t126 * t204) * t155 + (-t125 * t203 + t127 * t204) * t154 + (-t203 * t159 + t204 * t160 - t213 * t165 + t214 * t166 - t222 * t184 + t225 * t189 + Icges(3,1) + Icges(2,3)) * t205) * t205 / 0.2e1 + t205 * V_base(5) * (t309 * t223 - t308 * t226) + t205 * V_base(4) * (t308 * t223 + t309 * t226) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
