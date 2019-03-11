% Calculate kinetic energy for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:44
% EndTime: 2019-03-09 02:22:47
% DurationCPUTime: 3.11s
% Computational Cost: add. (1627->319), mult. (1566->450), div. (0->0), fcn. (1402->10), ass. (0->156)
t296 = Icges(3,4) + Icges(4,6);
t295 = Icges(3,1) + Icges(4,2);
t294 = -Icges(4,4) + Icges(3,5);
t293 = Icges(4,5) - Icges(3,6);
t292 = Icges(3,2) + Icges(4,3);
t223 = qJ(1) + pkin(10);
t214 = cos(t223);
t291 = t296 * t214;
t213 = sin(t223);
t290 = t296 * t213;
t289 = -t214 * t292 - t290;
t288 = t213 * t292 - t291;
t287 = t213 * t295 + t291;
t286 = t214 * t295 - t290;
t226 = sin(qJ(4));
t229 = cos(qJ(4));
t228 = cos(qJ(5));
t276 = pkin(5) * t228;
t283 = -pkin(9) * t229 + t226 * t276;
t271 = Icges(5,4) * t226;
t245 = Icges(5,2) * t229 + t271;
t132 = Icges(5,6) * t214 + t213 * t245;
t133 = Icges(5,6) * t213 - t214 * t245;
t270 = Icges(5,4) * t229;
t246 = Icges(5,1) * t226 + t270;
t134 = Icges(5,5) * t214 + t213 * t246;
t135 = Icges(5,5) * t213 - t214 * t246;
t188 = qJD(4) * t213 + V_base(5);
t189 = qJD(4) * t214 + V_base(4);
t193 = -Icges(5,2) * t226 + t270;
t196 = Icges(5,1) * t229 - t271;
t215 = V_base(6) + qJD(1);
t282 = (t132 * t229 + t134 * t226) * t189 + (t133 * t229 + t135 * t226) * t188 + (t193 * t229 + t196 * t226) * t215;
t227 = sin(qJ(1));
t280 = pkin(1) * t227;
t230 = cos(qJ(1));
t279 = pkin(1) * t230;
t278 = pkin(7) * t213;
t277 = pkin(7) * t214;
t275 = -pkin(6) - qJ(2);
t273 = Icges(2,4) * t227;
t225 = sin(qJ(5));
t267 = t213 * t225;
t266 = t213 * t229;
t265 = t214 * t225;
t264 = t214 * t229;
t224 = qJ(5) + qJ(6);
t217 = sin(t224);
t263 = t217 * t226;
t218 = cos(t224);
t262 = t218 * t226;
t261 = t225 * t226;
t260 = t226 * t228;
t259 = qJD(5) * t229;
t258 = t215 * t279 + V_base(2);
t257 = V_base(5) * pkin(6) + V_base(1);
t178 = pkin(2) * t213 - qJ(3) * t214;
t254 = -t178 - t280;
t181 = pkin(2) * t214 + qJ(3) * t213;
t253 = -t181 - t279;
t252 = V_base(5) * qJ(2) + t257;
t251 = V_base(4) * t280 + qJD(2) + V_base(3);
t151 = t214 * t259 + t188;
t200 = qJD(5) * t226 + t215;
t250 = qJD(3) * t213 + t252;
t249 = pkin(4) * t226 - pkin(8) * t229;
t248 = V_base(4) * t178 + t251;
t247 = rSges(5,1) * t226 + rSges(5,2) * t229;
t244 = Icges(5,5) * t226 + Icges(5,6) * t229;
t240 = V_base(5) * pkin(3) + t250;
t239 = t254 - t278;
t238 = -qJD(3) * t214 + t215 * t181 + t258;
t237 = (Icges(5,3) * t214 + t213 * t244) * t189 + (Icges(5,3) * t213 - t214 * t244) * t188 + (Icges(5,5) * t229 - Icges(5,6) * t226) * t215;
t236 = t215 * t277 + (-pkin(3) + t275) * V_base(4) + t238;
t235 = V_base(4) * t278 + (t253 - t277) * V_base(5) + t248;
t163 = t249 * t214;
t204 = t229 * pkin(4) + t226 * pkin(8);
t234 = t188 * t204 + (t163 + t239) * t215 + t240;
t162 = t249 * t213;
t233 = t215 * t162 - t189 * t204 + t236;
t232 = -t162 * t188 - t189 * t163 + t235;
t220 = Icges(2,4) * t230;
t203 = rSges(2,1) * t230 - rSges(2,2) * t227;
t202 = rSges(5,1) * t229 - rSges(5,2) * t226;
t201 = rSges(2,1) * t227 + rSges(2,2) * t230;
t198 = Icges(2,1) * t230 - t273;
t197 = Icges(2,1) * t227 + t220;
t195 = -Icges(2,2) * t227 + t220;
t194 = Icges(2,2) * t230 + t273;
t186 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t185 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t184 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t183 = rSges(3,1) * t214 - rSges(3,2) * t213;
t182 = -rSges(4,2) * t214 + rSges(4,3) * t213;
t180 = rSges(3,1) * t213 + rSges(3,2) * t214;
t179 = -rSges(4,2) * t213 - rSges(4,3) * t214;
t177 = qJD(6) * t226 + t200;
t161 = rSges(6,3) * t226 + (rSges(6,1) * t228 - rSges(6,2) * t225) * t229;
t159 = Icges(6,5) * t226 + (Icges(6,1) * t228 - Icges(6,4) * t225) * t229;
t158 = Icges(6,6) * t226 + (Icges(6,4) * t228 - Icges(6,2) * t225) * t229;
t157 = Icges(6,3) * t226 + (Icges(6,5) * t228 - Icges(6,6) * t225) * t229;
t156 = -t214 * t260 + t267;
t155 = t213 * t228 + t214 * t261;
t154 = t213 * t260 + t265;
t153 = -t213 * t261 + t214 * t228;
t152 = -t213 * t259 + t189;
t149 = rSges(7,3) * t226 + (rSges(7,1) * t218 - rSges(7,2) * t217) * t229;
t148 = Icges(7,5) * t226 + (Icges(7,1) * t218 - Icges(7,4) * t217) * t229;
t147 = Icges(7,6) * t226 + (Icges(7,4) * t218 - Icges(7,2) * t217) * t229;
t146 = Icges(7,3) * t226 + (Icges(7,5) * t218 - Icges(7,6) * t217) * t229;
t145 = t213 * t217 - t214 * t262;
t144 = t213 * t218 + t214 * t263;
t143 = t213 * t262 + t214 * t217;
t142 = -t213 * t263 + t214 * t218;
t141 = pkin(9) * t226 + t229 * t276;
t139 = rSges(5,3) * t213 - t214 * t247;
t138 = rSges(5,3) * t214 + t213 * t247;
t137 = V_base(5) * rSges(2,3) - t201 * t215 + t257;
t136 = t203 * t215 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t129 = t201 * V_base(4) - t203 * V_base(5) + V_base(3);
t128 = (-qJD(5) - qJD(6)) * t266 + t189;
t127 = qJD(6) * t264 + t151;
t125 = V_base(5) * rSges(3,3) + (-t180 - t280) * t215 + t252;
t124 = t183 * t215 + (-rSges(3,3) + t275) * V_base(4) + t258;
t123 = t180 * V_base(4) + (-t183 - t279) * V_base(5) + t251;
t122 = pkin(5) * t267 - t214 * t283;
t121 = pkin(5) * t265 + t213 * t283;
t120 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t264;
t119 = rSges(6,1) * t154 + rSges(6,2) * t153 - rSges(6,3) * t266;
t118 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t264;
t117 = Icges(6,1) * t154 + Icges(6,4) * t153 - Icges(6,5) * t266;
t116 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t264;
t115 = Icges(6,4) * t154 + Icges(6,2) * t153 - Icges(6,6) * t266;
t114 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t264;
t113 = Icges(6,5) * t154 + Icges(6,6) * t153 - Icges(6,3) * t266;
t112 = rSges(7,1) * t145 + rSges(7,2) * t144 + rSges(7,3) * t264;
t111 = rSges(7,1) * t143 + rSges(7,2) * t142 - rSges(7,3) * t266;
t110 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t264;
t109 = Icges(7,1) * t143 + Icges(7,4) * t142 - Icges(7,5) * t266;
t108 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t264;
t107 = Icges(7,4) * t143 + Icges(7,2) * t142 - Icges(7,6) * t266;
t106 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t264;
t105 = Icges(7,5) * t143 + Icges(7,6) * t142 - Icges(7,3) * t266;
t104 = V_base(5) * rSges(4,1) + (-t179 + t254) * t215 + t250;
t103 = t182 * t215 + (-rSges(4,1) + t275) * V_base(4) + t238;
t102 = t179 * V_base(4) + (-t182 + t253) * V_base(5) + t248;
t101 = t188 * t202 + (-t139 + t239) * t215 + t240;
t100 = t138 * t215 - t189 * t202 + t236;
t99 = -t138 * t188 + t139 * t189 + t235;
t98 = -t120 * t200 + t151 * t161 + t234;
t97 = t119 * t200 - t152 * t161 + t233;
t96 = -t119 * t151 + t120 * t152 + t232;
t95 = -t112 * t177 - t122 * t200 + t127 * t149 + t141 * t151 + t234;
t94 = t111 * t177 + t121 * t200 - t128 * t149 - t141 * t152 + t233;
t93 = -t111 * t127 + t112 * t128 - t121 * t151 + t122 * t152 + t232;
t1 = t152 * ((-t113 * t266 + t153 * t115 + t154 * t117) * t152 + (-t114 * t266 + t116 * t153 + t118 * t154) * t151 + (t153 * t158 + t154 * t159 - t157 * t266) * t200) / 0.2e1 + t128 * ((-t105 * t266 + t142 * t107 + t143 * t109) * t128 + (-t106 * t266 + t108 * t142 + t110 * t143) * t127 + (t142 * t147 + t143 * t148 - t146 * t266) * t177) / 0.2e1 + t151 * ((t113 * t264 + t115 * t155 + t117 * t156) * t152 + (t114 * t264 + t155 * t116 + t156 * t118) * t151 + (t155 * t158 + t156 * t159 + t157 * t264) * t200) / 0.2e1 + t127 * ((t105 * t264 + t107 * t144 + t109 * t145) * t128 + (t106 * t264 + t144 * t108 + t145 * t110) * t127 + (t144 * t147 + t145 * t148 + t146 * t264) * t177) / 0.2e1 + t189 * (t282 * t213 + t237 * t214) / 0.2e1 + t188 * (t237 * t213 - t282 * t214) / 0.2e1 + t200 * ((t113 * t152 + t114 * t151 + t157 * t200) * t226 + ((-t115 * t225 + t117 * t228) * t152 + (-t116 * t225 + t118 * t228) * t151 + (-t158 * t225 + t159 * t228) * t200) * t229) / 0.2e1 + m(7) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(2) * (t129 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(1) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + t177 * ((t105 * t128 + t106 * t127 + t146 * t177) * t226 + ((-t107 * t217 + t109 * t218) * t128 + (-t108 * t217 + t110 * t218) * t127 + (-t147 * t217 + t148 * t218) * t177) * t229) / 0.2e1 + ((-t132 * t226 + t134 * t229) * t189 + (-t133 * t226 + t135 * t229) * t188 + (-t226 * t193 + t229 * t196 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t215) * t215 / 0.2e1 + ((-t194 * t227 + t197 * t230 + t213 * t289 + t287 * t214 + Icges(1,4)) * V_base(5) + (-t227 * t195 + t230 * t198 + t288 * t213 + t286 * t214 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t230 * t194 + t227 * t197 + t287 * t213 - t289 * t214 + Icges(1,2)) * V_base(5) + (t195 * t230 + t198 * t227 + t213 * t286 - t214 * t288 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t215 * (Icges(2,5) * t227 + Icges(2,6) * t230 + t213 * t294 - t214 * t293) + V_base(4) * t215 * (Icges(2,5) * t230 - Icges(2,6) * t227 + t213 * t293 + t214 * t294) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
