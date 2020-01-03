% Calculate kinetic energy for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:09
% EndTime: 2019-12-31 21:22:12
% DurationCPUTime: 3.24s
% Computational Cost: add. (1490->325), mult. (1974->479), div. (0->0), fcn. (1924->10), ass. (0->146)
t261 = -Icges(5,3) - Icges(4,3);
t208 = qJ(3) + pkin(9);
t199 = sin(t208);
t200 = cos(t208);
t215 = cos(qJ(1));
t212 = sin(qJ(1));
t214 = cos(qJ(2));
t243 = t212 * t214;
t140 = -t199 * t243 - t200 * t215;
t141 = -t199 * t215 + t200 * t243;
t210 = sin(qJ(3));
t213 = cos(qJ(3));
t160 = -t210 * t243 - t213 * t215;
t247 = t210 * t215;
t161 = t213 * t243 - t247;
t211 = sin(qJ(2));
t246 = t211 * t212;
t260 = Icges(4,5) * t161 + Icges(5,5) * t141 + Icges(4,6) * t160 + Icges(5,6) * t140 - t261 * t246;
t242 = t214 * t215;
t142 = -t199 * t242 + t212 * t200;
t143 = t212 * t199 + t200 * t242;
t162 = -t210 * t242 + t212 * t213;
t244 = t212 * t210;
t163 = t213 * t242 + t244;
t245 = t211 * t215;
t259 = Icges(4,5) * t163 + Icges(5,5) * t143 + Icges(4,6) * t162 + Icges(5,6) * t142 - t261 * t245;
t258 = t261 * t214 + (Icges(4,5) * t213 + Icges(5,5) * t200 - Icges(4,6) * t210 - Icges(5,6) * t199) * t211;
t252 = t213 * pkin(3);
t250 = Icges(2,4) * t212;
t249 = Icges(3,4) * t211;
t248 = Icges(3,4) * t214;
t241 = pkin(4) * t200;
t239 = qJD(3) * t211;
t238 = qJD(4) * t211;
t237 = qJD(5) * t211;
t236 = V_base(5) * pkin(5) + V_base(1);
t191 = qJD(2) * t212 + V_base(4);
t202 = V_base(6) + qJD(1);
t233 = pkin(4) * t199;
t159 = t215 * t239 + t191;
t232 = pkin(2) * t214 + pkin(7) * t211;
t190 = -qJD(2) * t215 + V_base(5);
t231 = rSges(3,1) * t214 - rSges(3,2) * t211;
t230 = Icges(3,1) * t214 - t249;
t229 = -Icges(3,2) * t211 + t248;
t228 = Icges(3,5) * t214 - Icges(3,6) * t211;
t158 = t212 * t239 + t190;
t189 = pkin(1) * t215 + t212 * pkin(6);
t227 = -V_base(4) * pkin(5) + t202 * t189 + V_base(2);
t188 = t212 * pkin(1) - pkin(6) * t215;
t226 = V_base(4) * t188 - t189 * V_base(5) + V_base(3);
t225 = qJ(4) * t211 + t214 * t252;
t165 = t232 * t212;
t187 = pkin(2) * t211 - pkin(7) * t214;
t224 = t190 * t187 + (-t165 - t188) * t202 + t236;
t223 = (-Icges(3,3) * t215 + t212 * t228) * t190 + (Icges(3,3) * t212 + t215 * t228) * t191 + (Icges(3,5) * t211 + Icges(3,6) * t214) * t202;
t222 = pkin(8) * t211 + t214 * t241;
t166 = t232 * t215;
t221 = t202 * t166 - t187 * t191 + t227;
t128 = -qJ(4) * t214 + t211 * t252;
t220 = t158 * t128 + t215 * t238 + t224;
t219 = t191 * t165 - t166 * t190 + t226;
t116 = pkin(3) * t244 + t215 * t225;
t183 = -qJD(3) * t214 + t202;
t218 = t183 * t116 + t212 * t238 + t221;
t115 = -pkin(3) * t247 + t212 * t225;
t217 = -qJD(4) * t214 + t159 * t115 + t219;
t149 = -Icges(3,6) * t215 + t212 * t229;
t150 = Icges(3,6) * t212 + t215 * t229;
t152 = -Icges(3,5) * t215 + t212 * t230;
t153 = Icges(3,5) * t212 + t215 * t230;
t177 = Icges(3,2) * t214 + t249;
t180 = Icges(3,1) * t211 + t248;
t216 = (-t150 * t211 + t153 * t214) * t191 + (-t149 * t211 + t152 * t214) * t190 + (-t177 * t211 + t180 * t214) * t202;
t204 = Icges(2,4) * t215;
t201 = qJ(5) + t208;
t197 = cos(t201);
t196 = sin(t201);
t186 = rSges(2,1) * t215 - t212 * rSges(2,2);
t185 = t212 * rSges(2,1) + rSges(2,2) * t215;
t184 = rSges(3,1) * t211 + rSges(3,2) * t214;
t182 = Icges(2,1) * t215 - t250;
t181 = Icges(2,1) * t212 + t204;
t179 = -Icges(2,2) * t212 + t204;
t178 = Icges(2,2) * t215 + t250;
t172 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t171 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t170 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t167 = (-qJD(3) - qJD(5)) * t214 + t202;
t156 = t212 * rSges(3,3) + t215 * t231;
t155 = -rSges(3,3) * t215 + t212 * t231;
t154 = -rSges(4,3) * t214 + (rSges(4,1) * t213 - rSges(4,2) * t210) * t211;
t151 = -Icges(4,5) * t214 + (Icges(4,1) * t213 - Icges(4,4) * t210) * t211;
t148 = -Icges(4,6) * t214 + (Icges(4,4) * t213 - Icges(4,2) * t210) * t211;
t139 = t212 * t196 + t197 * t242;
t138 = -t196 * t242 + t212 * t197;
t137 = -t196 * t215 + t197 * t243;
t136 = -t196 * t243 - t197 * t215;
t135 = t215 * t237 + t159;
t134 = t212 * t237 + t158;
t132 = -rSges(5,3) * t214 + (rSges(5,1) * t200 - rSges(5,2) * t199) * t211;
t131 = -Icges(5,5) * t214 + (Icges(5,1) * t200 - Icges(5,4) * t199) * t211;
t130 = -Icges(5,6) * t214 + (Icges(5,4) * t200 - Icges(5,2) * t199) * t211;
t127 = -rSges(6,3) * t214 + (rSges(6,1) * t197 - rSges(6,2) * t196) * t211;
t126 = V_base(5) * rSges(2,3) - t185 * t202 + t236;
t125 = t186 * t202 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t124 = -Icges(6,5) * t214 + (Icges(6,1) * t197 - Icges(6,4) * t196) * t211;
t123 = -Icges(6,6) * t214 + (Icges(6,4) * t197 - Icges(6,2) * t196) * t211;
t122 = -Icges(6,3) * t214 + (Icges(6,5) * t197 - Icges(6,6) * t196) * t211;
t121 = t185 * V_base(4) - t186 * V_base(5) + V_base(3);
t120 = -pkin(8) * t214 + t211 * t241;
t118 = t163 * rSges(4,1) + t162 * rSges(4,2) + rSges(4,3) * t245;
t117 = rSges(4,1) * t161 + rSges(4,2) * t160 + rSges(4,3) * t246;
t114 = Icges(4,1) * t163 + Icges(4,4) * t162 + Icges(4,5) * t245;
t113 = Icges(4,1) * t161 + Icges(4,4) * t160 + Icges(4,5) * t246;
t112 = Icges(4,4) * t163 + Icges(4,2) * t162 + Icges(4,6) * t245;
t111 = Icges(4,4) * t161 + Icges(4,2) * t160 + Icges(4,6) * t246;
t108 = t143 * rSges(5,1) + t142 * rSges(5,2) + rSges(5,3) * t245;
t107 = rSges(5,1) * t141 + rSges(5,2) * t140 + rSges(5,3) * t246;
t106 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t245;
t105 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t246;
t104 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t245;
t103 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t246;
t99 = t139 * rSges(6,1) + t138 * rSges(6,2) + rSges(6,3) * t245;
t98 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t246;
t97 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t245;
t96 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t246;
t95 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t245;
t94 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t246;
t93 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t245;
t92 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t246;
t90 = t184 * t190 + (-t155 - t188) * t202 + t236;
t89 = t156 * t202 - t184 * t191 + t227;
t88 = t212 * t233 + t215 * t222;
t87 = t212 * t222 - t215 * t233;
t86 = t155 * t191 - t156 * t190 + t226;
t85 = -t117 * t183 + t154 * t158 + t224;
t84 = t118 * t183 - t154 * t159 + t221;
t83 = t117 * t159 - t118 * t158 + t219;
t82 = t132 * t158 + (-t107 - t115) * t183 + t220;
t81 = t108 * t183 + (-t128 - t132) * t159 + t218;
t80 = t107 * t159 + (-t108 - t116) * t158 + t217;
t79 = t120 * t158 + t127 * t134 - t167 * t98 + (-t115 - t87) * t183 + t220;
t78 = -t127 * t135 + t167 * t99 + t183 * t88 + (-t120 - t128) * t159 + t218;
t77 = -t134 * t99 + t135 * t98 + t159 * t87 + (-t116 - t88) * t158 + t217;
t1 = m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t191 * (t212 * t223 + t215 * t216) / 0.2e1 + t190 * (t212 * t216 - t223 * t215) / 0.2e1 + m(4) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t135 * ((t138 * t95 + t139 * t97 + t93 * t245) * t135 + (t138 * t94 + t139 * t96 + t245 * t92) * t134 + (t122 * t245 + t138 * t123 + t139 * t124) * t167) / 0.2e1 + t134 * ((t136 * t95 + t137 * t97 + t246 * t93) * t135 + (t136 * t94 + t137 * t96 + t92 * t246) * t134 + (t122 * t246 + t123 * t136 + t124 * t137) * t167) / 0.2e1 + t167 * ((-t122 * t167 - t92 * t134 - t93 * t135) * t214 + ((-t196 * t95 + t197 * t97) * t135 + (-t196 * t94 + t197 * t96) * t134 + (-t123 * t196 + t124 * t197) * t167) * t211) / 0.2e1 + ((t130 * t140 + t131 * t141 + t148 * t160 + t151 * t161 + t246 * t258) * t183 + (t104 * t140 + t106 * t141 + t112 * t160 + t114 * t161 + t246 * t259) * t159 + (t103 * t140 + t105 * t141 + t111 * t160 + t113 * t161 + t260 * t246) * t158) * t158 / 0.2e1 + ((t142 * t130 + t143 * t131 + t162 * t148 + t163 * t151 + t245 * t258) * t183 + (t142 * t104 + t143 * t106 + t162 * t112 + t163 * t114 + t259 * t245) * t159 + (t142 * t103 + t143 * t105 + t162 * t111 + t163 * t113 + t245 * t260) * t158) * t159 / 0.2e1 + ((-t158 * t260 - t259 * t159 - t258 * t183) * t214 + ((-t130 * t199 + t131 * t200 - t148 * t210 + t151 * t213) * t183 + (-t104 * t199 + t106 * t200 - t112 * t210 + t114 * t213) * t159 + (-t103 * t199 + t105 * t200 - t111 * t210 + t113 * t213) * t158) * t211) * t183 / 0.2e1 + ((t150 * t214 + t153 * t211) * t191 + (t149 * t214 + t152 * t211) * t190 + (t177 * t214 + t180 * t211 + Icges(2,3)) * t202) * t202 / 0.2e1 + ((-t212 * t178 + t181 * t215 + Icges(1,4)) * V_base(5) + (-t212 * t179 + t182 * t215 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t178 * t215 + t212 * t181 + Icges(1,2)) * V_base(5) + (t179 * t215 + t212 * t182 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t202 * (Icges(2,5) * t215 - Icges(2,6) * t212) + V_base(5) * t202 * (Icges(2,5) * t212 + Icges(2,6) * t215) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
