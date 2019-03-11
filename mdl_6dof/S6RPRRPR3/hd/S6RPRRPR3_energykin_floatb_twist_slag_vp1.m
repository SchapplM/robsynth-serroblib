% Calculate kinetic energy for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:21
% EndTime: 2019-03-09 05:04:24
% DurationCPUTime: 3.24s
% Computational Cost: add. (2089->327), mult. (2390->460), div. (0->0), fcn. (2446->10), ass. (0->152)
t298 = Icges(5,1) + Icges(6,1);
t297 = -Icges(5,4) + Icges(6,5);
t296 = Icges(6,4) + Icges(5,5);
t295 = Icges(5,2) + Icges(6,3);
t294 = -Icges(6,6) + Icges(5,6);
t293 = -Icges(5,3) - Icges(6,2);
t234 = qJ(1) + pkin(10);
t228 = sin(t234);
t229 = cos(t234);
t240 = cos(qJ(4));
t236 = sin(qJ(4));
t241 = cos(qJ(3));
t268 = t236 * t241;
t170 = t228 * t268 + t229 * t240;
t267 = t240 * t241;
t171 = t228 * t267 - t229 * t236;
t237 = sin(qJ(3));
t270 = t228 * t237;
t292 = t295 * t170 + t297 * t171 - t294 * t270;
t172 = -t228 * t240 + t229 * t268;
t173 = t228 * t236 + t229 * t267;
t269 = t229 * t237;
t291 = t295 * t172 + t297 * t173 - t294 * t269;
t290 = -t294 * t170 + t296 * t171 - t293 * t270;
t289 = -t294 * t172 + t296 * t173 - t293 * t269;
t288 = t297 * t170 + t298 * t171 + t296 * t270;
t287 = t297 * t172 + t298 * t173 + t296 * t269;
t286 = t294 * t241 + (t295 * t236 + t297 * t240) * t237;
t285 = t293 * t241 + (-t294 * t236 + t296 * t240) * t237;
t284 = -t296 * t241 + (t297 * t236 + t298 * t240) * t237;
t238 = sin(qJ(1));
t277 = pkin(1) * t238;
t242 = cos(qJ(1));
t276 = pkin(1) * t242;
t275 = -pkin(6) - qJ(2);
t274 = Icges(2,4) * t238;
t273 = Icges(3,4) * t228;
t272 = Icges(4,4) * t237;
t271 = Icges(4,4) * t241;
t266 = qJD(4) * t237;
t265 = qJD(6) * t237;
t230 = V_base(6) + qJD(1);
t264 = t230 * t276 + V_base(2);
t263 = V_base(5) * pkin(6) + V_base(1);
t205 = qJD(3) * t228 + V_base(4);
t198 = pkin(2) * t228 - pkin(7) * t229;
t260 = -t198 - t277;
t259 = V_base(5) * qJ(2) + t263;
t258 = V_base(4) * t277 + qJD(2) + V_base(3);
t169 = t229 * t266 + t205;
t257 = pkin(3) * t241 + pkin(8) * t237;
t204 = -qJD(3) * t229 + V_base(5);
t256 = rSges(4,1) * t241 - rSges(4,2) * t237;
t255 = Icges(4,1) * t241 - t272;
t254 = -Icges(4,2) * t237 + t271;
t253 = Icges(4,5) * t241 - Icges(4,6) * t237;
t168 = t228 * t266 + t204;
t252 = (-Icges(4,3) * t229 + t228 * t253) * t204 + (Icges(4,3) * t228 + t229 * t253) * t205 + (Icges(4,5) * t237 + Icges(4,6) * t241) * t230;
t199 = pkin(2) * t229 + pkin(7) * t228;
t251 = t230 * t199 + t275 * V_base(4) + t264;
t183 = t257 * t228;
t221 = pkin(3) * t237 - pkin(8) * t241;
t250 = t204 * t221 + (-t183 + t260) * t230 + t259;
t249 = V_base(4) * t198 + (-t199 - t276) * V_base(5) + t258;
t184 = t257 * t229;
t248 = t230 * t184 - t205 * t221 + t251;
t188 = (pkin(4) * t240 + qJ(5) * t236) * t237;
t247 = qJD(5) * t172 + t168 * t188 + t250;
t146 = pkin(4) * t173 + qJ(5) * t172;
t217 = -qJD(4) * t241 + t230;
t246 = qJD(5) * t170 + t217 * t146 + t248;
t245 = t205 * t183 - t204 * t184 + t249;
t145 = pkin(4) * t171 + qJ(5) * t170;
t244 = qJD(5) * t237 * t236 + t169 * t145 + t245;
t156 = -Icges(4,6) * t229 + t228 * t254;
t157 = Icges(4,6) * t228 + t229 * t254;
t158 = -Icges(4,5) * t229 + t228 * t255;
t159 = Icges(4,5) * t228 + t229 * t255;
t209 = Icges(4,2) * t241 + t272;
t212 = Icges(4,1) * t237 + t271;
t243 = (-t157 * t237 + t159 * t241) * t205 + (-t156 * t237 + t158 * t241) * t204 + (-t209 * t237 + t212 * t241) * t230;
t239 = cos(qJ(6));
t235 = sin(qJ(6));
t232 = Icges(2,4) * t242;
t227 = Icges(3,4) * t229;
t220 = rSges(2,1) * t242 - t238 * rSges(2,2);
t219 = t238 * rSges(2,1) + rSges(2,2) * t242;
t218 = rSges(4,1) * t237 + rSges(4,2) * t241;
t214 = Icges(2,1) * t242 - t274;
t213 = Icges(2,1) * t238 + t232;
t211 = -Icges(2,2) * t238 + t232;
t210 = Icges(2,2) * t242 + t274;
t203 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t202 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t201 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t200 = pkin(5) * t237 * t240 + pkin(9) * t241;
t197 = rSges(3,1) * t229 - rSges(3,2) * t228;
t196 = rSges(3,1) * t228 + rSges(3,2) * t229;
t195 = (-qJD(4) + qJD(6)) * t241 + t230;
t194 = Icges(3,1) * t229 - t273;
t193 = Icges(3,1) * t228 + t227;
t192 = -Icges(3,2) * t228 + t227;
t191 = Icges(3,2) * t229 + t273;
t186 = (t235 * t236 + t239 * t240) * t237;
t185 = (-t235 * t240 + t236 * t239) * t237;
t181 = -rSges(5,3) * t241 + (rSges(5,1) * t240 - rSges(5,2) * t236) * t237;
t180 = -rSges(6,2) * t241 + (rSges(6,1) * t240 + rSges(6,3) * t236) * t237;
t163 = rSges(4,3) * t228 + t229 * t256;
t162 = -rSges(4,3) * t229 + t228 * t256;
t161 = V_base(5) * rSges(2,3) - t219 * t230 + t263;
t160 = t220 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t153 = t219 * V_base(4) - t220 * V_base(5) + V_base(3);
t152 = -t229 * t265 + t169;
t151 = -t228 * t265 + t168;
t149 = pkin(5) * t173 - pkin(9) * t269;
t148 = pkin(5) * t171 - pkin(9) * t270;
t144 = V_base(5) * rSges(3,3) + (-t196 - t277) * t230 + t259;
t143 = t197 * t230 + (-rSges(3,3) + t275) * V_base(4) + t264;
t142 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t241;
t141 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t241;
t140 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t241;
t139 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t241;
t138 = t172 * t235 + t173 * t239;
t137 = t172 * t239 - t173 * t235;
t136 = t170 * t235 + t171 * t239;
t135 = t170 * t239 - t171 * t235;
t134 = V_base(4) * t196 + (-t197 - t276) * V_base(5) + t258;
t133 = rSges(5,1) * t173 - rSges(5,2) * t172 + rSges(5,3) * t269;
t132 = rSges(6,1) * t173 + rSges(6,2) * t269 + rSges(6,3) * t172;
t131 = rSges(5,1) * t171 - rSges(5,2) * t170 + rSges(5,3) * t270;
t130 = rSges(6,1) * t171 + rSges(6,2) * t270 + rSges(6,3) * t170;
t115 = t204 * t218 + (-t162 + t260) * t230 + t259;
t114 = t163 * t230 - t205 * t218 + t251;
t113 = rSges(7,1) * t138 + rSges(7,2) * t137 - rSges(7,3) * t269;
t112 = rSges(7,1) * t136 + rSges(7,2) * t135 - rSges(7,3) * t270;
t111 = Icges(7,1) * t138 + Icges(7,4) * t137 - Icges(7,5) * t269;
t110 = Icges(7,1) * t136 + Icges(7,4) * t135 - Icges(7,5) * t270;
t109 = Icges(7,4) * t138 + Icges(7,2) * t137 - Icges(7,6) * t269;
t108 = Icges(7,4) * t136 + Icges(7,2) * t135 - Icges(7,6) * t270;
t107 = Icges(7,5) * t138 + Icges(7,6) * t137 - Icges(7,3) * t269;
t106 = Icges(7,5) * t136 + Icges(7,6) * t135 - Icges(7,3) * t270;
t105 = t205 * t162 - t204 * t163 + t249;
t104 = -t131 * t217 + t168 * t181 + t250;
t103 = t133 * t217 - t169 * t181 + t248;
t102 = t169 * t131 - t168 * t133 + t245;
t101 = t168 * t180 + (-t130 - t145) * t217 + t247;
t100 = t132 * t217 + (-t180 - t188) * t169 + t246;
t99 = t169 * t130 + (-t132 - t146) * t168 + t244;
t98 = -t112 * t195 + t142 * t151 + t168 * t200 + (-t145 - t148) * t217 + t247;
t97 = t113 * t195 - t142 * t152 + t149 * t217 + (-t188 - t200) * t169 + t246;
t96 = t152 * t112 - t151 * t113 + t169 * t148 + (-t146 - t149) * t168 + t244;
t1 = t151 * ((-t107 * t270 + t109 * t135 + t111 * t136) * t152 + (-t106 * t270 + t135 * t108 + t136 * t110) * t151 + (t135 * t140 + t136 * t141 - t139 * t270) * t195) / 0.2e1 + t195 * ((t107 * t241 + t109 * t185 + t111 * t186) * t152 + (t106 * t241 + t108 * t185 + t110 * t186) * t151 + (t241 * t139 + t185 * t140 + t186 * t141) * t195) / 0.2e1 + m(1) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(2) * (t153 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + t205 * (t252 * t228 + t243 * t229) / 0.2e1 + t204 * (t243 * t228 - t252 * t229) / 0.2e1 + t152 * ((-t107 * t269 + t137 * t109 + t138 * t111) * t152 + (-t106 * t269 + t108 * t137 + t110 * t138) * t151 + (t137 * t140 + t138 * t141 - t139 * t269) * t195) / 0.2e1 + ((t170 * t286 + t171 * t284 + t270 * t285) * t217 + (t170 * t291 + t171 * t287 + t270 * t289) * t169 + (t292 * t170 + t288 * t171 + t290 * t270) * t168) * t168 / 0.2e1 + ((t172 * t286 + t173 * t284 + t269 * t285) * t217 + (t291 * t172 + t287 * t173 + t289 * t269) * t169 + (t172 * t292 + t288 * t173 + t290 * t269) * t168) * t169 / 0.2e1 + ((-t168 * t290 - t169 * t289 - t217 * t285) * t241 + ((t236 * t286 + t240 * t284) * t217 + (t236 * t291 + t240 * t287) * t169 + (t236 * t292 + t288 * t240) * t168) * t237) * t217 / 0.2e1 + ((t157 * t241 + t159 * t237) * t205 + (t156 * t241 + t158 * t237) * t204 + (t241 * t209 + t237 * t212 + Icges(2,3) + Icges(3,3)) * t230) * t230 / 0.2e1 + ((-t191 * t228 + t193 * t229 - t238 * t210 + t213 * t242 + Icges(1,4)) * V_base(5) + (-t192 * t228 + t194 * t229 - t238 * t211 + t214 * t242 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t191 * t229 + t193 * t228 + t210 * t242 + t238 * t213 + Icges(1,2)) * V_base(5) + (t192 * t229 + t194 * t228 + t211 * t242 + t238 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t230 * (Icges(2,5) * t238 + Icges(3,5) * t228 + Icges(2,6) * t242 + Icges(3,6) * t229) + V_base(4) * t230 * (Icges(2,5) * t242 + Icges(3,5) * t229 - Icges(2,6) * t238 - Icges(3,6) * t228) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
