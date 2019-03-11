% Calculate kinetic energy for
% S6RPPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:27
% EndTime: 2019-03-09 02:30:30
% DurationCPUTime: 3.00s
% Computational Cost: add. (1020->304), mult. (1578->432), div. (0->0), fcn. (1416->8), ass. (0->148)
t279 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t278 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t277 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t276 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t275 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t212 = sin(qJ(1));
t274 = t279 * t212;
t215 = cos(qJ(1));
t273 = t279 * t215;
t272 = t278 * t215 - t274;
t271 = t278 * t212 + t273;
t270 = -t275 * t215 - t274;
t269 = t275 * t212 - t273;
t266 = -pkin(2) - pkin(6);
t264 = pkin(7) * t212;
t263 = pkin(7) * t215;
t213 = cos(qJ(5));
t262 = pkin(5) * t213;
t211 = sin(qJ(4));
t259 = Icges(5,4) * t211;
t214 = cos(qJ(4));
t258 = Icges(5,4) * t214;
t254 = qJ(3) * t212;
t253 = qJ(3) * t215;
t210 = sin(qJ(5));
t252 = t210 * t212;
t251 = t210 * t215;
t250 = t211 * t212;
t249 = t211 * t215;
t248 = t212 * t213;
t247 = t212 * t214;
t246 = t213 * t215;
t245 = t214 * t215;
t244 = qJD(5) * t214;
t179 = pkin(1) * t212 - qJ(2) * t215;
t243 = V_base(4) * t179 + V_base(3);
t242 = V_base(5) * pkin(6) + V_base(1);
t190 = qJD(4) * t215 + V_base(5);
t196 = V_base(6) + qJD(1);
t239 = V_base(4) * t254 + t243;
t238 = qJD(2) * t212 + t242;
t237 = -t179 - t254;
t184 = pkin(1) * t215 + qJ(2) * t212;
t236 = -t184 - t253;
t235 = (-qJD(5) - qJD(6)) * t214;
t178 = qJD(5) * t211 + t196;
t234 = pkin(4) * t211 - pkin(8) * t214;
t191 = -qJD(4) * t212 + V_base(4);
t233 = rSges(5,1) * t211 + rSges(5,2) * t214;
t232 = Icges(5,1) * t211 + t258;
t231 = Icges(5,2) * t214 + t259;
t230 = Icges(5,5) * t211 + Icges(5,6) * t214;
t229 = -qJD(2) * t215 + t196 * t184 + V_base(2);
t228 = V_base(5) * pkin(2) + qJD(3) * t215 + t238;
t227 = t237 - t263;
t226 = V_base(5) * pkin(3) + t228;
t225 = qJD(3) * t212 + t196 * t253 + t229;
t224 = -pkin(9) * t214 + t211 * t262;
t223 = (Icges(5,3) * t215 + t212 * t230) * t190 + (-Icges(5,3) * t212 + t215 * t230) * t191 + (Icges(5,5) * t214 - Icges(5,6) * t211) * t196;
t222 = V_base(4) * t263 + t239 + (t236 + t264) * V_base(5);
t221 = (-pkin(3) + t266) * V_base(4) + t225;
t150 = t234 * t212;
t151 = t234 * t215;
t220 = t191 * t150 - t151 * t190 + t222;
t189 = t214 * pkin(4) + t211 * pkin(8);
t219 = t190 * t189 + (-t150 + t227) * t196 + t226;
t218 = -t189 * t191 + t221 + (t151 - t264) * t196;
t130 = Icges(5,6) * t215 + t212 * t231;
t131 = -Icges(5,6) * t212 + t215 * t231;
t133 = Icges(5,5) * t215 + t212 * t232;
t134 = -Icges(5,5) * t212 + t215 * t232;
t168 = -Icges(5,2) * t211 + t258;
t175 = Icges(5,1) * t214 - t259;
t217 = (t131 * t214 + t134 * t211) * t191 + (t130 * t214 + t133 * t211) * t190 + (t168 * t214 + t175 * t211) * t196;
t209 = qJ(5) + qJ(6);
t203 = cos(t209);
t202 = sin(t209);
t187 = rSges(2,1) * t215 - rSges(2,2) * t212;
t186 = -rSges(3,2) * t215 + rSges(3,3) * t212;
t185 = -rSges(4,2) * t215 + rSges(4,3) * t212;
t183 = rSges(5,1) * t214 - rSges(5,2) * t211;
t182 = rSges(2,1) * t212 + rSges(2,2) * t215;
t181 = -rSges(3,2) * t212 - rSges(3,3) * t215;
t180 = rSges(4,2) * t212 + rSges(4,3) * t215;
t156 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t155 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t154 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t152 = qJD(6) * t211 + t178;
t148 = t211 * t246 - t252;
t147 = -t210 * t249 - t248;
t146 = t211 * t248 + t251;
t145 = -t210 * t250 + t246;
t144 = -t215 * t244 + t191;
t143 = -t212 * t244 + t190;
t141 = -t202 * t212 + t203 * t249;
t140 = -t202 * t249 - t203 * t212;
t139 = t202 * t215 + t203 * t250;
t138 = -t202 * t250 + t203 * t215;
t137 = -rSges(5,3) * t212 + t215 * t233;
t136 = rSges(6,3) * t211 + (rSges(6,1) * t213 - rSges(6,2) * t210) * t214;
t135 = rSges(5,3) * t215 + t212 * t233;
t132 = Icges(6,5) * t211 + (Icges(6,1) * t213 - Icges(6,4) * t210) * t214;
t129 = Icges(6,6) * t211 + (Icges(6,4) * t213 - Icges(6,2) * t210) * t214;
t126 = Icges(6,3) * t211 + (Icges(6,5) * t213 - Icges(6,6) * t210) * t214;
t124 = t215 * t235 + t191;
t123 = t212 * t235 + t190;
t122 = rSges(7,3) * t211 + (rSges(7,1) * t203 - rSges(7,2) * t202) * t214;
t120 = Icges(7,5) * t211 + (Icges(7,1) * t203 - Icges(7,4) * t202) * t214;
t119 = Icges(7,6) * t211 + (Icges(7,4) * t203 - Icges(7,2) * t202) * t214;
t118 = Icges(7,3) * t211 + (Icges(7,5) * t203 - Icges(7,6) * t202) * t214;
t117 = pkin(9) * t211 + t214 * t262;
t116 = V_base(5) * rSges(2,3) - t182 * t196 + t242;
t115 = t187 * t196 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t114 = t182 * V_base(4) - t187 * V_base(5) + V_base(3);
t113 = -pkin(5) * t252 + t215 * t224;
t112 = pkin(5) * t251 + t212 * t224;
t111 = rSges(6,1) * t148 + rSges(6,2) * t147 - rSges(6,3) * t245;
t110 = rSges(6,1) * t146 + rSges(6,2) * t145 - rSges(6,3) * t247;
t109 = Icges(6,1) * t148 + Icges(6,4) * t147 - Icges(6,5) * t245;
t108 = Icges(6,1) * t146 + Icges(6,4) * t145 - Icges(6,5) * t247;
t107 = Icges(6,4) * t148 + Icges(6,2) * t147 - Icges(6,6) * t245;
t106 = Icges(6,4) * t146 + Icges(6,2) * t145 - Icges(6,6) * t247;
t105 = Icges(6,5) * t148 + Icges(6,6) * t147 - Icges(6,3) * t245;
t104 = Icges(6,5) * t146 + Icges(6,6) * t145 - Icges(6,3) * t247;
t103 = V_base(5) * rSges(3,1) + (-t179 - t181) * t196 + t238;
t102 = t186 * t196 + (-rSges(3,1) - pkin(6)) * V_base(4) + t229;
t101 = rSges(7,1) * t141 + rSges(7,2) * t140 - rSges(7,3) * t245;
t100 = rSges(7,1) * t139 + rSges(7,2) * t138 - rSges(7,3) * t247;
t99 = Icges(7,1) * t141 + Icges(7,4) * t140 - Icges(7,5) * t245;
t98 = Icges(7,1) * t139 + Icges(7,4) * t138 - Icges(7,5) * t247;
t97 = Icges(7,4) * t141 + Icges(7,2) * t140 - Icges(7,6) * t245;
t96 = Icges(7,4) * t139 + Icges(7,2) * t138 - Icges(7,6) * t247;
t95 = Icges(7,5) * t141 + Icges(7,6) * t140 - Icges(7,3) * t245;
t94 = Icges(7,5) * t139 + Icges(7,6) * t138 - Icges(7,3) * t247;
t93 = t181 * V_base(4) + (-t184 - t186) * V_base(5) + t243;
t92 = V_base(5) * rSges(4,1) + (-t185 + t237) * t196 + t228;
t91 = t180 * t196 + (-rSges(4,1) + t266) * V_base(4) + t225;
t90 = t185 * V_base(4) + (-t180 + t236) * V_base(5) + t239;
t89 = t183 * t190 + (-t135 + t227) * t196 + t226;
t88 = -t183 * t191 + (t137 - t264) * t196 + t221;
t87 = t135 * t191 - t137 * t190 + t222;
t86 = -t110 * t178 + t136 * t143 + t219;
t85 = t111 * t178 - t136 * t144 + t218;
t84 = t110 * t144 - t111 * t143 + t220;
t83 = -t100 * t152 - t112 * t178 + t117 * t143 + t122 * t123 + t219;
t82 = t101 * t152 + t113 * t178 - t117 * t144 - t122 * t124 + t218;
t81 = t100 * t124 - t101 * t123 + t112 * t144 - t113 * t143 + t220;
t1 = t143 * ((-t105 * t247 + t107 * t145 + t109 * t146) * t144 + (-t104 * t247 + t145 * t106 + t146 * t108) * t143 + (-t126 * t247 + t129 * t145 + t132 * t146) * t178) / 0.2e1 + t123 * ((t138 * t97 + t139 * t99 - t247 * t95) * t124 + (t138 * t96 + t139 * t98 - t94 * t247) * t123 + (-t118 * t247 + t119 * t138 + t120 * t139) * t152) / 0.2e1 + t144 * ((-t105 * t245 + t147 * t107 + t148 * t109) * t144 + (-t104 * t245 + t106 * t147 + t108 * t148) * t143 + (-t126 * t245 + t129 * t147 + t132 * t148) * t178) / 0.2e1 + t124 * ((t140 * t97 + t141 * t99 - t95 * t245) * t124 + (t140 * t96 + t141 * t98 - t245 * t94) * t123 + (-t118 * t245 + t119 * t140 + t120 * t141) * t152) / 0.2e1 + t190 * (t217 * t212 + t223 * t215) / 0.2e1 + t191 * (-t223 * t212 + t217 * t215) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t178 * ((t104 * t143 + t105 * t144 + t126 * t178) * t211 + ((-t107 * t210 + t109 * t213) * t144 + (-t106 * t210 + t108 * t213) * t143 + (-t129 * t210 + t132 * t213) * t178) * t214) / 0.2e1 + t152 * ((t118 * t152 + t94 * t123 + t95 * t124) * t211 + ((-t202 * t97 + t203 * t99) * t124 + (-t202 * t96 + t203 * t98) * t123 + (-t119 * t202 + t120 * t203) * t152) * t214) / 0.2e1 + m(1) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(2) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((-t131 * t211 + t134 * t214) * t191 + (-t130 * t211 + t133 * t214) * t190 + (-t168 * t211 + t175 * t214 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t196) * t196 / 0.2e1 + ((t270 * t212 + t271 * t215 + Icges(1,4)) * V_base(5) + (t269 * t212 + t272 * t215 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t271 * t212 - t270 * t215 + Icges(1,2)) * V_base(5) + (t212 * t272 - t269 * t215 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t196 * (t276 * t212 + t277 * t215) + V_base(5) * t196 * (t277 * t212 - t276 * t215) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
