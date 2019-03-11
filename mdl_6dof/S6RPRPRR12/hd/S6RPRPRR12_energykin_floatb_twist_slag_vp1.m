% Calculate kinetic energy for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:06
% EndTime: 2019-03-09 04:18:10
% DurationCPUTime: 3.79s
% Computational Cost: add. (1100->313), mult. (1755->442), div. (0->0), fcn. (1593->8), ass. (0->152)
t317 = -Icges(4,4) - Icges(5,6);
t316 = Icges(4,1) + Icges(5,2);
t315 = Icges(4,2) + Icges(5,3);
t222 = cos(qJ(3));
t314 = t317 * t222;
t219 = sin(qJ(3));
t313 = t317 * t219;
t312 = -Icges(5,4) + Icges(4,5);
t311 = Icges(5,5) - Icges(4,6);
t310 = t222 * t315 - t313;
t309 = t219 * t316 - t314;
t308 = Icges(2,4) + Icges(3,6);
t307 = Icges(2,1) + Icges(3,2);
t306 = Icges(5,1) + Icges(4,3);
t305 = -Icges(3,4) + Icges(2,5);
t304 = Icges(3,5) - Icges(2,6);
t303 = Icges(2,2) + Icges(3,3);
t220 = sin(qJ(1));
t223 = cos(qJ(1));
t302 = -t310 * t220 + t311 * t223;
t301 = t311 * t220 + t310 * t223;
t300 = t309 * t220 + t312 * t223;
t299 = t312 * t220 - t309 * t223;
t298 = t219 * t315 + t314;
t297 = t222 * t316 + t313;
t296 = t312 * t219 - t311 * t222;
t295 = t308 * t223;
t294 = t308 * t220;
t200 = qJD(3) * t220 + V_base(5);
t201 = qJD(3) * t223 + V_base(4);
t206 = V_base(6) + qJD(1);
t293 = t200 * (-t299 * t219 + t301 * t222) + t201 * (-t300 * t219 + t302 * t222) - t206 * (t297 * t219 - t298 * t222);
t292 = -t303 * t223 - t294;
t291 = t303 * t220 - t295;
t290 = t307 * t220 + t295;
t289 = t307 * t223 - t294;
t286 = (t311 * t219 + t312 * t222) * t206 + (t296 * t220 + t306 * t223) * t201 + (t306 * t220 - t296 * t223) * t200;
t218 = sin(qJ(5));
t278 = pkin(5) * t218;
t282 = pkin(9) * t219 - t222 * t278;
t277 = pkin(7) * t220;
t276 = pkin(7) * t223;
t275 = t222 * pkin(8);
t221 = cos(qJ(5));
t274 = pkin(5) * t221;
t265 = t219 * t220;
t264 = t219 * t223;
t263 = t220 * t222;
t262 = t222 * t223;
t261 = qJD(4) * t222;
t260 = qJD(5) * t219;
t190 = pkin(1) * t220 - qJ(2) * t223;
t259 = V_base(4) * t190 + V_base(3);
t258 = V_base(5) * pkin(6) + V_base(1);
t254 = -t190 - t277;
t253 = qJD(2) * t220 + t258;
t155 = t220 * t260 + t201;
t189 = qJD(5) * t222 + t206;
t248 = pkin(3) * t219 - qJ(4) * t222;
t162 = t248 * t223;
t252 = t162 + t254;
t251 = V_base(5) * pkin(2) + t253;
t250 = rSges(4,1) * t219 + rSges(4,2) * t222;
t249 = rSges(5,2) * t219 + rSges(5,3) * t222;
t196 = pkin(1) * t223 + qJ(2) * t220;
t235 = -qJD(2) * t223 + t206 * t196 + V_base(2);
t193 = pkin(3) * t222 + qJ(4) * t219;
t232 = t200 * t193 - t220 * t261 + t251;
t231 = V_base(4) * t277 + (-t196 - t276) * V_base(5) + t259;
t230 = t206 * t276 + (-pkin(2) - pkin(6)) * V_base(4) + t235;
t229 = qJD(4) * t219 - t201 * t162 + t231;
t161 = t248 * t220;
t228 = t206 * t161 + t223 * t261 + t230;
t166 = t220 * pkin(4) - pkin(8) * t264;
t167 = t223 * pkin(4) + pkin(8) * t265;
t227 = t201 * t166 + (-t161 - t167) * t200 + t229;
t226 = t200 * t275 + (-t166 + t252) * t206 + t232;
t225 = t206 * t167 + (-t193 - t275) * t201 + t228;
t217 = qJ(5) + qJ(6);
t213 = cos(t217);
t212 = sin(t217);
t198 = rSges(2,1) * t223 - rSges(2,2) * t220;
t197 = -rSges(3,2) * t223 + rSges(3,3) * t220;
t195 = rSges(4,1) * t222 - rSges(4,2) * t219;
t194 = -rSges(5,2) * t222 + rSges(5,3) * t219;
t192 = rSges(2,1) * t220 + rSges(2,2) * t223;
t191 = -rSges(3,2) * t220 - rSges(3,3) * t223;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = qJD(6) * t222 + t189;
t159 = -t218 * t263 + t221 * t223;
t158 = -t218 * t223 - t221 * t263;
t157 = t218 * t262 + t220 * t221;
t156 = -t218 * t220 + t221 * t262;
t154 = -t223 * t260 + t200;
t152 = pkin(9) * t222 + t219 * t278;
t150 = -t212 * t263 + t213 * t223;
t149 = -t212 * t223 - t213 * t263;
t148 = t212 * t262 + t213 * t220;
t147 = -t212 * t220 + t213 * t262;
t146 = rSges(5,1) * t223 - t220 * t249;
t145 = rSges(5,1) * t220 + t223 * t249;
t144 = rSges(4,3) * t220 - t223 * t250;
t143 = rSges(4,3) * t223 + t220 * t250;
t142 = rSges(6,3) * t222 + (rSges(6,1) * t218 + rSges(6,2) * t221) * t219;
t133 = Icges(6,5) * t222 + (Icges(6,1) * t218 + Icges(6,4) * t221) * t219;
t130 = Icges(6,6) * t222 + (Icges(6,4) * t218 + Icges(6,2) * t221) * t219;
t127 = Icges(6,3) * t222 + (Icges(6,5) * t218 + Icges(6,6) * t221) * t219;
t124 = qJD(6) * t265 + t155;
t123 = (-qJD(5) - qJD(6)) * t264 + t200;
t122 = rSges(7,3) * t222 + (rSges(7,1) * t212 + rSges(7,2) * t213) * t219;
t120 = Icges(7,5) * t222 + (Icges(7,1) * t212 + Icges(7,4) * t213) * t219;
t119 = Icges(7,6) * t222 + (Icges(7,4) * t212 + Icges(7,2) * t213) * t219;
t118 = Icges(7,3) * t222 + (Icges(7,5) * t212 + Icges(7,6) * t213) * t219;
t117 = V_base(5) * rSges(2,3) - t192 * t206 + t258;
t116 = t198 * t206 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t115 = t192 * V_base(4) - t198 * V_base(5) + V_base(3);
t114 = t220 * t282 + t274 * t223;
t113 = t274 * t220 - t223 * t282;
t112 = rSges(6,1) * t159 + rSges(6,2) * t158 + rSges(6,3) * t265;
t111 = rSges(6,1) * t157 + rSges(6,2) * t156 - rSges(6,3) * t264;
t110 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t265;
t109 = Icges(6,1) * t157 + Icges(6,4) * t156 - Icges(6,5) * t264;
t108 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t265;
t107 = Icges(6,4) * t157 + Icges(6,2) * t156 - Icges(6,6) * t264;
t106 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t265;
t105 = Icges(6,5) * t157 + Icges(6,6) * t156 - Icges(6,3) * t264;
t104 = V_base(5) * rSges(3,1) + (-t190 - t191) * t206 + t253;
t103 = t197 * t206 + (-rSges(3,1) - pkin(6)) * V_base(4) + t235;
t102 = rSges(7,1) * t150 + rSges(7,2) * t149 + rSges(7,3) * t265;
t101 = rSges(7,1) * t148 + rSges(7,2) * t147 - rSges(7,3) * t264;
t100 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t265;
t99 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t264;
t98 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t265;
t97 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t264;
t96 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t265;
t95 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t264;
t94 = t191 * V_base(4) + (-t196 - t197) * V_base(5) + t259;
t93 = t195 * t200 + (-t144 + t254) * t206 + t251;
t92 = t143 * t206 - t195 * t201 + t230;
t91 = -t143 * t200 + t144 * t201 + t231;
t90 = t194 * t200 + (-t145 + t252) * t206 + t232;
t89 = t146 * t206 + (-t193 - t194) * t201 + t228;
t88 = t145 * t201 + (-t146 - t161) * t200 + t229;
t87 = -t111 * t189 + t142 * t154 + t226;
t86 = t112 * t189 - t142 * t155 + t225;
t85 = t111 * t155 - t112 * t154 + t227;
t84 = -t101 * t163 - t113 * t189 + t122 * t123 + t152 * t154 + t226;
t83 = t102 * t163 + t114 * t189 - t122 * t124 - t152 * t155 + t225;
t82 = t101 * t124 - t102 * t123 + t113 * t155 - t114 * t154 + t227;
t1 = t155 * ((t106 * t265 + t158 * t108 + t159 * t110) * t155 + (t105 * t265 + t107 * t158 + t109 * t159) * t154 + (t127 * t265 + t130 * t158 + t133 * t159) * t189) / 0.2e1 + t124 * ((t150 * t100 + t149 * t98 + t96 * t265) * t124 + (t149 * t97 + t150 * t99 + t265 * t95) * t123 + (t118 * t265 + t119 * t149 + t120 * t150) * t163) / 0.2e1 + t154 * ((-t106 * t264 + t108 * t156 + t110 * t157) * t155 + (-t105 * t264 + t156 * t107 + t157 * t109) * t154 + (-t127 * t264 + t130 * t156 + t133 * t157) * t189) / 0.2e1 + t123 * ((t100 * t148 + t147 * t98 - t264 * t96) * t124 + (t147 * t97 + t148 * t99 - t95 * t264) * t123 + (-t118 * t264 + t119 * t147 + t120 * t148) * t163) / 0.2e1 + m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t103 ^ 2 + t104 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t189 * ((t105 * t154 + t106 * t155 + t127 * t189) * t222 + ((t108 * t221 + t110 * t218) * t155 + (t107 * t221 + t109 * t218) * t154 + (t130 * t221 + t133 * t218) * t189) * t219) / 0.2e1 + t163 * ((t118 * t163 + t95 * t123 + t96 * t124) * t222 + ((t100 * t212 + t213 * t98) * t124 + (t212 * t99 + t213 * t97) * t123 + (t119 * t213 + t120 * t212) * t163) * t219) / 0.2e1 + (t286 * t220 + t293 * t223) * t200 / 0.2e1 + (-t293 * t220 + t286 * t223) * t201 / 0.2e1 + ((t220 * t292 + t290 * t223 + Icges(1,4)) * V_base(5) + (t291 * t220 + t289 * t223 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t290 * t220 - t292 * t223 + Icges(1,2)) * V_base(5) + (t220 * t289 - t223 * t291 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t302 * t219 + t300 * t222) * t201 + (t301 * t219 + t299 * t222) * t200 + (t298 * t219 + t297 * t222 + Icges(3,1) + Icges(2,3)) * t206) * t206 / 0.2e1 + t206 * V_base(5) * (t305 * t220 - t304 * t223) + t206 * V_base(4) * (t304 * t220 + t305 * t223) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
