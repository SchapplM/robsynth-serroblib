% Calculate kinetic energy for
% S6RPPRRR5
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:14
% DurationCPUTime: 2.44s
% Computational Cost: add. (1009->281), mult. (1304->390), div. (0->0), fcn. (1086->8), ass. (0->145)
t278 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t277 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t276 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t275 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t274 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t205 = sin(qJ(1));
t273 = t278 * t205;
t208 = cos(qJ(1));
t272 = t278 * t208;
t271 = t277 * t208 - t273;
t270 = t277 * t205 + t272;
t269 = -t274 * t208 - t273;
t268 = t274 * t205 - t272;
t265 = -pkin(2) - pkin(6);
t204 = sin(qJ(4));
t263 = pkin(4) * t204;
t207 = cos(qJ(4));
t262 = pkin(4) * t207;
t261 = t205 * pkin(7);
t260 = t208 * pkin(7);
t257 = Icges(5,4) * t204;
t256 = Icges(5,4) * t207;
t202 = qJ(4) + qJ(5);
t195 = sin(t202);
t255 = Icges(6,4) * t195;
t196 = cos(t202);
t254 = Icges(6,4) * t196;
t250 = qJ(3) * t205;
t249 = qJ(3) * t208;
t248 = t196 * t205;
t247 = t196 * t208;
t203 = sin(qJ(6));
t246 = t203 * t205;
t245 = t203 * t208;
t206 = cos(qJ(6));
t244 = t205 * t206;
t243 = t206 * t208;
t242 = qJD(6) * t196;
t174 = pkin(1) * t205 - qJ(2) * t208;
t241 = V_base(4) * t174 + V_base(3);
t240 = V_base(5) * pkin(6) + V_base(1);
t184 = qJD(4) * t208 + V_base(5);
t189 = V_base(6) + qJD(1);
t237 = V_base(4) * t250 + t241;
t236 = qJD(2) * t205 + t240;
t235 = -t174 - t250;
t179 = pkin(1) * t208 + qJ(2) * t205;
t234 = -t179 - t249;
t145 = qJD(5) * t208 + t184;
t233 = pkin(5) * t195 - pkin(9) * t196;
t232 = rSges(5,1) * t204 + rSges(5,2) * t207;
t231 = rSges(6,1) * t195 + rSges(6,2) * t196;
t230 = Icges(5,1) * t204 + t256;
t229 = Icges(6,1) * t195 + t254;
t228 = Icges(5,2) * t207 + t257;
t227 = Icges(6,2) * t196 + t255;
t226 = Icges(5,5) * t204 + Icges(5,6) * t207;
t225 = Icges(6,5) * t195 + Icges(6,6) * t196;
t224 = -qJD(2) * t208 + t189 * t179 + V_base(2);
t223 = V_base(5) * pkin(2) + qJD(3) * t208 + t236;
t222 = t235 - t260;
t146 = V_base(4) + (-qJD(4) - qJD(5)) * t205;
t221 = V_base(5) * pkin(3) + t223;
t134 = pkin(8) * t208 + t205 * t263;
t220 = -t134 + t222;
t219 = qJD(3) * t205 + t189 * t249 + t224;
t218 = t184 * t262 + t221;
t217 = (Icges(6,3) * t208 + t205 * t225) * t145 + (-Icges(6,3) * t205 + t208 * t225) * t146 + (Icges(6,5) * t196 - Icges(6,6) * t195) * t189;
t185 = -qJD(4) * t205 + V_base(4);
t216 = (Icges(5,3) * t208 + t205 * t226) * t184 + (-Icges(5,3) * t205 + t208 * t226) * t185 + (Icges(5,5) * t207 - Icges(5,6) * t204) * t189;
t215 = V_base(4) * t260 + t237 + (t234 + t261) * V_base(5);
t214 = (-pkin(3) + t265) * V_base(4) + t219;
t133 = -pkin(8) * t205 + t208 * t263;
t213 = -t133 * t184 + t185 * t134 + t215;
t212 = t189 * t133 - t185 * t262 + t214;
t114 = Icges(6,6) * t208 + t205 * t227;
t115 = -Icges(6,6) * t205 + t208 * t227;
t116 = Icges(6,5) * t208 + t205 * t229;
t117 = -Icges(6,5) * t205 + t208 * t229;
t141 = -Icges(6,2) * t195 + t254;
t142 = Icges(6,1) * t196 - t255;
t211 = (t115 * t196 + t117 * t195) * t146 + (t114 * t196 + t116 * t195) * t145 + (t141 * t196 + t142 * t195) * t189;
t125 = Icges(5,6) * t208 + t205 * t228;
t126 = -Icges(5,6) * t205 + t208 * t228;
t127 = Icges(5,5) * t208 + t205 * t230;
t128 = -Icges(5,5) * t205 + t208 * t230;
t164 = -Icges(5,2) * t204 + t256;
t171 = Icges(5,1) * t207 - t257;
t210 = (t126 * t207 + t128 * t204) * t185 + (t125 * t207 + t127 * t204) * t184 + (t164 * t207 + t171 * t204) * t189;
t182 = rSges(2,1) * t208 - rSges(2,2) * t205;
t181 = -rSges(3,2) * t208 + rSges(3,3) * t205;
t180 = -rSges(4,2) * t208 + rSges(4,3) * t205;
t178 = rSges(5,1) * t207 - rSges(5,2) * t204;
t177 = rSges(2,1) * t205 + rSges(2,2) * t208;
t176 = -rSges(3,2) * t205 - rSges(3,3) * t208;
t175 = rSges(4,2) * t205 + rSges(4,3) * t208;
t152 = qJD(6) * t195 + t189;
t151 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t150 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t149 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t144 = pkin(5) * t196 + pkin(9) * t195;
t143 = rSges(6,1) * t196 - rSges(6,2) * t195;
t138 = t195 * t243 - t246;
t137 = -t195 * t245 - t244;
t136 = t195 * t244 + t245;
t135 = -t195 * t246 + t243;
t132 = t233 * t208;
t131 = t233 * t205;
t130 = -rSges(5,3) * t205 + t208 * t232;
t129 = rSges(5,3) * t208 + t205 * t232;
t122 = -t208 * t242 + t146;
t121 = -t205 * t242 + t145;
t119 = -rSges(6,3) * t205 + t208 * t231;
t118 = rSges(6,3) * t208 + t205 * t231;
t111 = rSges(7,3) * t195 + (rSges(7,1) * t206 - rSges(7,2) * t203) * t196;
t110 = Icges(7,5) * t195 + (Icges(7,1) * t206 - Icges(7,4) * t203) * t196;
t109 = Icges(7,6) * t195 + (Icges(7,4) * t206 - Icges(7,2) * t203) * t196;
t108 = Icges(7,3) * t195 + (Icges(7,5) * t206 - Icges(7,6) * t203) * t196;
t107 = V_base(5) * rSges(2,3) - t177 * t189 + t240;
t106 = t182 * t189 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t104 = t177 * V_base(4) - t182 * V_base(5) + V_base(3);
t103 = V_base(5) * rSges(3,1) + (-t174 - t176) * t189 + t236;
t102 = t181 * t189 + (-rSges(3,1) - pkin(6)) * V_base(4) + t224;
t101 = rSges(7,1) * t138 + rSges(7,2) * t137 - rSges(7,3) * t247;
t100 = rSges(7,1) * t136 + rSges(7,2) * t135 - rSges(7,3) * t248;
t99 = Icges(7,1) * t138 + Icges(7,4) * t137 - Icges(7,5) * t247;
t98 = Icges(7,1) * t136 + Icges(7,4) * t135 - Icges(7,5) * t248;
t97 = Icges(7,4) * t138 + Icges(7,2) * t137 - Icges(7,6) * t247;
t96 = Icges(7,4) * t136 + Icges(7,2) * t135 - Icges(7,6) * t248;
t95 = Icges(7,5) * t138 + Icges(7,6) * t137 - Icges(7,3) * t247;
t94 = Icges(7,5) * t136 + Icges(7,6) * t135 - Icges(7,3) * t248;
t93 = t176 * V_base(4) + (-t179 - t181) * V_base(5) + t241;
t92 = V_base(5) * rSges(4,1) + (-t180 + t235) * t189 + t223;
t91 = t175 * t189 + (-rSges(4,1) + t265) * V_base(4) + t219;
t90 = t180 * V_base(4) + (-t175 + t234) * V_base(5) + t237;
t89 = t178 * t184 + (-t129 + t222) * t189 + t221;
t88 = -t178 * t185 + (t130 - t261) * t189 + t214;
t87 = t129 * t185 - t130 * t184 + t215;
t86 = t143 * t145 + (-t118 + t220) * t189 + t218;
t85 = -t143 * t146 + (t119 - t261) * t189 + t212;
t84 = t118 * t146 - t119 * t145 + t213;
t83 = -t100 * t152 + t111 * t121 + t144 * t145 + (-t131 + t220) * t189 + t218;
t82 = t101 * t152 - t111 * t122 - t144 * t146 + (t132 - t261) * t189 + t212;
t81 = t100 * t122 - t101 * t121 + t131 * t146 - t132 * t145 + t213;
t1 = t122 * ((t137 * t97 + t138 * t99 - t95 * t247) * t122 + (t137 * t96 + t138 * t98 - t247 * t94) * t121 + (-t108 * t247 + t109 * t137 + t110 * t138) * t152) / 0.2e1 + t121 * ((t135 * t97 + t136 * t99 - t248 * t95) * t122 + (t135 * t96 + t136 * t98 - t94 * t248) * t121 + (-t108 * t248 + t109 * t135 + t110 * t136) * t152) / 0.2e1 + t184 * (t205 * t210 + t208 * t216) / 0.2e1 + t185 * (-t205 * t216 + t208 * t210) / 0.2e1 + t145 * (t211 * t205 + t217 * t208) / 0.2e1 + t146 * (-t217 * t205 + t211 * t208) / 0.2e1 + t152 * ((t108 * t152 + t94 * t121 + t95 * t122) * t195 + ((-t203 * t97 + t206 * t99) * t122 + (-t203 * t96 + t206 * t98) * t121 + (-t109 * t203 + t110 * t206) * t152) * t196) / 0.2e1 + m(1) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + ((t205 * t269 + t208 * t270 + Icges(1,4)) * V_base(5) + (t268 * t205 + t271 * t208 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t270 * t205 - t269 * t208 + Icges(1,2)) * V_base(5) + (t205 * t271 - t268 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t126 * t204 + t128 * t207) * t185 + (-t125 * t204 + t127 * t207) * t184 + (-t115 * t195 + t117 * t196) * t146 + (-t114 * t195 + t116 * t196) * t145 + (-t141 * t195 + t142 * t196 - t164 * t204 + t171 * t207 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t189) * t189 / 0.2e1 + t189 * V_base(4) * (t275 * t205 + t276 * t208) + t189 * V_base(5) * (t276 * t205 - t275 * t208) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
