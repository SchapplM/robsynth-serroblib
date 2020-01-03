% Calculate kinetic energy for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:31
% EndTime: 2019-12-31 19:31:33
% DurationCPUTime: 2.57s
% Computational Cost: add. (1415->306), mult. (1580->431), div. (0->0), fcn. (1474->10), ass. (0->148)
t264 = Icges(3,3) + Icges(4,3);
t199 = qJ(2) + pkin(8);
t190 = sin(t199);
t192 = cos(t199);
t204 = sin(qJ(2));
t206 = cos(qJ(2));
t263 = Icges(3,5) * t206 + Icges(4,5) * t192 - Icges(3,6) * t204 - Icges(4,6) * t190;
t205 = sin(qJ(1));
t207 = cos(qJ(1));
t249 = Icges(4,4) * t192;
t222 = -Icges(4,2) * t190 + t249;
t130 = -Icges(4,6) * t207 + t205 * t222;
t131 = Icges(4,6) * t205 + t207 * t222;
t250 = Icges(4,4) * t190;
t224 = Icges(4,1) * t192 - t250;
t132 = -Icges(4,5) * t207 + t205 * t224;
t133 = Icges(4,5) * t205 + t207 * t224;
t251 = Icges(3,4) * t206;
t223 = -Icges(3,2) * t204 + t251;
t143 = -Icges(3,6) * t207 + t205 * t223;
t144 = Icges(3,6) * t205 + t207 * t223;
t252 = Icges(3,4) * t204;
t225 = Icges(3,1) * t206 - t252;
t145 = -Icges(3,5) * t207 + t205 * t225;
t146 = Icges(3,5) * t205 + t207 * t225;
t159 = Icges(4,2) * t192 + t250;
t160 = Icges(4,1) * t190 + t249;
t172 = Icges(3,2) * t206 + t252;
t175 = Icges(3,1) * t204 + t251;
t185 = -qJD(2) * t207 + V_base(5);
t186 = qJD(2) * t205 + V_base(4);
t193 = V_base(6) + qJD(1);
t262 = (-t159 * t190 + t160 * t192 - t172 * t204 + t175 * t206) * t193 + (-t131 * t190 + t133 * t192 - t144 * t204 + t146 * t206) * t186 + (-t130 * t190 + t132 * t192 - t143 * t204 + t145 * t206) * t185;
t261 = (Icges(3,5) * t204 + Icges(4,5) * t190 + Icges(3,6) * t206 + Icges(4,6) * t192) * t193 + (t264 * t205 + t263 * t207) * t186 + (t263 * t205 - t264 * t207) * t185;
t257 = pkin(2) * t204;
t256 = pkin(2) * t206;
t201 = cos(pkin(9));
t255 = pkin(4) * t201;
t253 = Icges(2,4) * t205;
t248 = t190 * t205;
t247 = t190 * t207;
t246 = t192 * t207;
t200 = sin(pkin(9));
t245 = t200 * t207;
t244 = t201 * t207;
t198 = pkin(9) + qJ(5);
t189 = sin(t198);
t243 = t205 * t189;
t191 = cos(t198);
t242 = t205 * t191;
t241 = t205 * t200;
t240 = t205 * t201;
t126 = -qJ(3) * t207 + t205 * t256;
t183 = t205 * pkin(1) - pkin(6) * t207;
t238 = -t126 - t183;
t127 = qJ(3) * t205 + t207 * t256;
t226 = pkin(3) * t192 + qJ(4) * t190;
t148 = t226 * t207;
t237 = -t127 - t148;
t236 = qJD(4) * t190;
t235 = qJD(5) * t190;
t234 = V_base(5) * pkin(5) + V_base(1);
t147 = t226 * t205;
t231 = -t147 + t238;
t161 = pkin(3) * t190 - qJ(4) * t192;
t230 = -t161 - t257;
t229 = qJD(3) * t205 + t185 * t257 + t234;
t228 = rSges(3,1) * t206 - rSges(3,2) * t204;
t227 = rSges(4,1) * t192 - rSges(4,2) * t190;
t184 = pkin(1) * t207 + t205 * pkin(6);
t219 = -V_base(4) * pkin(5) + t193 * t184 + V_base(2);
t218 = V_base(4) * t183 - t184 * V_base(5) + V_base(3);
t217 = t185 * t161 + t207 * t236 + t229;
t216 = t186 * t126 + t218;
t213 = pkin(7) * t190 + t192 * t255;
t212 = -qJD(3) * t207 + t193 * t127 + t219;
t211 = -qJD(4) * t192 + t186 * t147 + t216;
t210 = t193 * t148 + t205 * t236 + t212;
t196 = Icges(2,4) * t207;
t182 = rSges(2,1) * t207 - t205 * rSges(2,2);
t181 = t205 * rSges(2,1) + rSges(2,2) * t207;
t180 = rSges(3,1) * t204 + rSges(3,2) * t206;
t177 = Icges(2,1) * t207 - t253;
t176 = Icges(2,1) * t205 + t196;
t174 = -Icges(2,2) * t205 + t196;
t173 = Icges(2,2) * t207 + t253;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t165 = -qJD(5) * t192 + t193;
t162 = rSges(4,1) * t190 + rSges(4,2) * t192;
t156 = t207 * t235 + t186;
t155 = t205 * t235 + t185;
t154 = t192 * t244 + t241;
t153 = -t192 * t245 + t240;
t152 = t192 * t240 - t245;
t151 = -t192 * t241 - t244;
t150 = t205 * rSges(3,3) + t207 * t228;
t149 = -rSges(3,3) * t207 + t205 * t228;
t140 = t191 * t246 + t243;
t139 = -t189 * t246 + t242;
t138 = -t189 * t207 + t192 * t242;
t137 = -t191 * t207 - t192 * t243;
t135 = t205 * rSges(4,3) + t207 * t227;
t134 = -rSges(4,3) * t207 + t205 * t227;
t124 = V_base(5) * rSges(2,3) - t181 * t193 + t234;
t123 = t182 * t193 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t122 = -rSges(5,3) * t192 + (rSges(5,1) * t201 - rSges(5,2) * t200) * t190;
t121 = -Icges(5,5) * t192 + (Icges(5,1) * t201 - Icges(5,4) * t200) * t190;
t120 = -Icges(5,6) * t192 + (Icges(5,4) * t201 - Icges(5,2) * t200) * t190;
t119 = -Icges(5,3) * t192 + (Icges(5,5) * t201 - Icges(5,6) * t200) * t190;
t117 = t181 * V_base(4) - t182 * V_base(5) + V_base(3);
t115 = -rSges(6,3) * t192 + (rSges(6,1) * t191 - rSges(6,2) * t189) * t190;
t114 = -Icges(6,5) * t192 + (Icges(6,1) * t191 - Icges(6,4) * t189) * t190;
t113 = -Icges(6,6) * t192 + (Icges(6,4) * t191 - Icges(6,2) * t189) * t190;
t112 = -Icges(6,3) * t192 + (Icges(6,5) * t191 - Icges(6,6) * t189) * t190;
t110 = -pkin(7) * t192 + t190 * t255;
t109 = t154 * rSges(5,1) + t153 * rSges(5,2) + rSges(5,3) * t247;
t108 = rSges(5,1) * t152 + rSges(5,2) * t151 + rSges(5,3) * t248;
t107 = Icges(5,1) * t154 + Icges(5,4) * t153 + Icges(5,5) * t247;
t106 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t248;
t105 = Icges(5,4) * t154 + Icges(5,2) * t153 + Icges(5,6) * t247;
t104 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t248;
t103 = Icges(5,5) * t154 + Icges(5,6) * t153 + Icges(5,3) * t247;
t102 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t248;
t101 = pkin(4) * t241 + t207 * t213;
t100 = -pkin(4) * t245 + t205 * t213;
t99 = t140 * rSges(6,1) + t139 * rSges(6,2) + rSges(6,3) * t247;
t98 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t248;
t97 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t247;
t96 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t248;
t95 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t247;
t94 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t248;
t93 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t247;
t92 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t248;
t91 = t180 * t185 + (-t149 - t183) * t193 + t234;
t90 = t150 * t193 - t180 * t186 + t219;
t89 = t149 * t186 - t150 * t185 + t218;
t88 = t162 * t185 + (-t134 + t238) * t193 + t229;
t87 = t193 * t135 + (-t162 - t257) * t186 + t212;
t86 = t134 * t186 + (-t127 - t135) * t185 + t216;
t85 = t122 * t185 + (-t108 + t231) * t193 + t217;
t84 = t193 * t109 + (-t122 + t230) * t186 + t210;
t83 = t108 * t186 + (-t109 + t237) * t185 + t211;
t82 = t110 * t185 + t115 * t155 - t165 * t98 + (-t100 + t231) * t193 + t217;
t81 = t193 * t101 - t156 * t115 + t165 * t99 + (-t110 + t230) * t186 + t210;
t80 = t100 * t186 - t155 * t99 + t156 * t98 + (-t101 + t237) * t185 + t211;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t156 * ((t139 * t95 + t140 * t97 + t93 * t247) * t156 + (t139 * t94 + t140 * t96 + t247 * t92) * t155 + (t112 * t247 + t139 * t113 + t140 * t114) * t165) / 0.2e1 + t155 * ((t137 * t95 + t138 * t97 + t248 * t93) * t156 + (t137 * t94 + t138 * t96 + t92 * t248) * t155 + (t112 * t248 + t113 * t137 + t114 * t138) * t165) / 0.2e1 + t165 * ((-t112 * t165 - t92 * t155 - t93 * t156) * t192 + ((-t189 * t95 + t191 * t97) * t156 + (-t189 * t94 + t191 * t96) * t155 + (-t113 * t189 + t114 * t191) * t165) * t190) / 0.2e1 + ((-t205 * t173 + t176 * t207 + Icges(1,4)) * V_base(5) + (-t205 * t174 + t177 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t173 * t207 + t205 * t176 + Icges(1,2)) * V_base(5) + (t174 * t207 + t205 * t177 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t103 * t248 + t105 * t151 + t107 * t152) * t186 + (t102 * t248 + t104 * t151 + t106 * t152) * t185 + (t119 * t248 + t120 * t151 + t121 * t152) * t193 - t261 * t207 + t262 * t205) * t185 / 0.2e1 + ((t103 * t247 + t153 * t105 + t154 * t107) * t186 + (t102 * t247 + t153 * t104 + t154 * t106) * t185 + (t119 * t247 + t153 * t120 + t154 * t121) * t193 + t262 * t207 + t261 * t205) * t186 / 0.2e1 + ((t144 * t206 + t146 * t204 + (t131 - t103) * t192 + (-t105 * t200 + t107 * t201 + t133) * t190) * t186 + (t143 * t206 + t145 * t204 + (t130 - t102) * t192 + (-t104 * t200 + t106 * t201 + t132) * t190) * t185 + (t172 * t206 + t175 * t204 + Icges(2,3) + (t159 - t119) * t192 + (-t120 * t200 + t121 * t201 + t160) * t190) * t193) * t193 / 0.2e1 + t193 * V_base(4) * (Icges(2,5) * t207 - Icges(2,6) * t205) + V_base(5) * t193 * (Icges(2,5) * t205 + Icges(2,6) * t207) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
