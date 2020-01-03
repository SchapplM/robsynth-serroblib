% Calculate kinetic energy for
% S5RRRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:39
% EndTime: 2019-12-31 22:27:42
% DurationCPUTime: 2.96s
% Computational Cost: add. (1526->331), mult. (2034->512), div. (0->0), fcn. (1984->10), ass. (0->152)
t215 = cos(qJ(3));
t256 = t215 * pkin(3);
t214 = sin(qJ(1));
t254 = Icges(2,4) * t214;
t213 = sin(qJ(2));
t253 = Icges(3,4) * t213;
t216 = cos(qJ(2));
t252 = Icges(3,4) * t216;
t212 = sin(qJ(3));
t251 = t212 * t214;
t217 = cos(qJ(1));
t250 = t212 * t217;
t249 = t213 * t214;
t248 = t213 * t217;
t247 = t214 * t216;
t246 = t216 * t217;
t211 = qJ(3) + qJ(4);
t205 = cos(t211);
t245 = pkin(4) * t205;
t243 = qJD(3) * t213;
t242 = qJD(4) * t213;
t241 = qJD(5) * t213;
t240 = -qJD(3) - qJD(4);
t239 = V_base(5) * pkin(5) + V_base(1);
t194 = qJD(2) * t214 + V_base(4);
t202 = V_base(6) + qJD(1);
t204 = sin(t211);
t236 = pkin(4) * t204;
t162 = t217 * t243 + t194;
t235 = pkin(2) * t216 + pkin(7) * t213;
t193 = -qJD(2) * t217 + V_base(5);
t234 = rSges(3,1) * t216 - rSges(3,2) * t213;
t137 = t217 * t242 + t162;
t233 = Icges(3,1) * t216 - t253;
t232 = -Icges(3,2) * t213 + t252;
t231 = Icges(3,5) * t216 - Icges(3,6) * t213;
t161 = t214 * t243 + t193;
t192 = pkin(1) * t217 + pkin(6) * t214;
t230 = -V_base(4) * pkin(5) + t192 * t202 + V_base(2);
t191 = pkin(1) * t214 - pkin(6) * t217;
t229 = t191 * V_base(4) - t192 * V_base(5) + V_base(3);
t136 = t214 * t242 + t161;
t228 = pkin(8) * t213 + t216 * t256;
t168 = t235 * t214;
t190 = pkin(2) * t213 - pkin(7) * t216;
t227 = t193 * t190 + (-t168 - t191) * t202 + t239;
t226 = (-Icges(3,3) * t217 + t214 * t231) * t193 + (Icges(3,3) * t214 + t217 * t231) * t194 + (Icges(3,5) * t213 + Icges(3,6) * t216) * t202;
t225 = pkin(9) * t213 + t216 * t245;
t169 = t235 * t217;
t224 = t169 * t202 - t190 * t194 + t230;
t223 = t168 * t194 - t169 * t193 + t229;
t117 = -pkin(3) * t250 + t214 * t228;
t130 = -pkin(8) * t216 + t213 * t256;
t186 = -qJD(3) * t216 + t202;
t222 = -t117 * t186 + t130 * t161 + t227;
t118 = pkin(3) * t251 + t217 * t228;
t221 = t118 * t186 - t130 * t162 + t224;
t220 = t117 * t162 - t118 * t161 + t223;
t147 = -Icges(3,6) * t217 + t214 * t232;
t148 = Icges(3,6) * t214 + t217 * t232;
t150 = -Icges(3,5) * t217 + t214 * t233;
t151 = Icges(3,5) * t214 + t217 * t233;
t180 = Icges(3,2) * t216 + t253;
t183 = Icges(3,1) * t213 + t252;
t219 = (-t148 * t213 + t151 * t216) * t194 + (-t147 * t213 + t150 * t216) * t193 + (-t180 * t213 + t183 * t216) * t202;
t207 = qJ(5) + t211;
t206 = Icges(2,4) * t217;
t200 = cos(t207);
t199 = sin(t207);
t189 = rSges(2,1) * t217 - rSges(2,2) * t214;
t188 = rSges(2,1) * t214 + rSges(2,2) * t217;
t187 = rSges(3,1) * t213 + rSges(3,2) * t216;
t185 = Icges(2,1) * t217 - t254;
t184 = Icges(2,1) * t214 + t206;
t182 = -Icges(2,2) * t214 + t206;
t181 = Icges(2,2) * t217 + t254;
t175 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t174 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t173 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t170 = t216 * t240 + t202;
t166 = t215 * t246 + t251;
t165 = -t212 * t246 + t214 * t215;
t164 = t215 * t247 - t250;
t163 = -t212 * t247 - t215 * t217;
t159 = (-qJD(5) + t240) * t216 + t202;
t158 = t204 * t214 + t205 * t246;
t157 = -t204 * t246 + t205 * t214;
t156 = -t204 * t217 + t205 * t247;
t155 = -t204 * t247 - t205 * t217;
t154 = rSges(3,3) * t214 + t217 * t234;
t153 = -rSges(3,3) * t217 + t214 * t234;
t152 = -rSges(4,3) * t216 + (rSges(4,1) * t215 - rSges(4,2) * t212) * t213;
t149 = -Icges(4,5) * t216 + (Icges(4,1) * t215 - Icges(4,4) * t212) * t213;
t146 = -Icges(4,6) * t216 + (Icges(4,4) * t215 - Icges(4,2) * t212) * t213;
t143 = -Icges(4,3) * t216 + (Icges(4,5) * t215 - Icges(4,6) * t212) * t213;
t141 = t199 * t214 + t200 * t246;
t140 = -t199 * t246 + t200 * t214;
t139 = -t199 * t217 + t200 * t247;
t138 = -t199 * t247 - t200 * t217;
t135 = -rSges(5,3) * t216 + (rSges(5,1) * t205 - rSges(5,2) * t204) * t213;
t133 = -Icges(5,5) * t216 + (Icges(5,1) * t205 - Icges(5,4) * t204) * t213;
t132 = -Icges(5,6) * t216 + (Icges(5,4) * t205 - Icges(5,2) * t204) * t213;
t131 = -Icges(5,3) * t216 + (Icges(5,5) * t205 - Icges(5,6) * t204) * t213;
t129 = -rSges(6,3) * t216 + (rSges(6,1) * t200 - rSges(6,2) * t199) * t213;
t128 = -Icges(6,5) * t216 + (Icges(6,1) * t200 - Icges(6,4) * t199) * t213;
t127 = -Icges(6,6) * t216 + (Icges(6,4) * t200 - Icges(6,2) * t199) * t213;
t126 = -Icges(6,3) * t216 + (Icges(6,5) * t200 - Icges(6,6) * t199) * t213;
t125 = V_base(5) * rSges(2,3) - t188 * t202 + t239;
t124 = t189 * t202 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t123 = t188 * V_base(4) - t189 * V_base(5) + V_base(3);
t122 = t217 * t241 + t137;
t121 = t214 * t241 + t136;
t120 = -pkin(9) * t216 + t213 * t245;
t116 = rSges(4,1) * t166 + rSges(4,2) * t165 + rSges(4,3) * t248;
t115 = rSges(4,1) * t164 + rSges(4,2) * t163 + rSges(4,3) * t249;
t114 = Icges(4,1) * t166 + Icges(4,4) * t165 + Icges(4,5) * t248;
t113 = Icges(4,1) * t164 + Icges(4,4) * t163 + Icges(4,5) * t249;
t112 = Icges(4,4) * t166 + Icges(4,2) * t165 + Icges(4,6) * t248;
t111 = Icges(4,4) * t164 + Icges(4,2) * t163 + Icges(4,6) * t249;
t110 = Icges(4,5) * t166 + Icges(4,6) * t165 + Icges(4,3) * t248;
t109 = Icges(4,5) * t164 + Icges(4,6) * t163 + Icges(4,3) * t249;
t108 = rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t248;
t107 = rSges(5,1) * t156 + rSges(5,2) * t155 + rSges(5,3) * t249;
t106 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t248;
t105 = Icges(5,1) * t156 + Icges(5,4) * t155 + Icges(5,5) * t249;
t104 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t248;
t103 = Icges(5,4) * t156 + Icges(5,2) * t155 + Icges(5,6) * t249;
t102 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t248;
t101 = Icges(5,5) * t156 + Icges(5,6) * t155 + Icges(5,3) * t249;
t99 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t248;
t98 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t249;
t97 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t248;
t96 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t249;
t95 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t248;
t94 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t249;
t93 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t248;
t92 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t249;
t90 = t187 * t193 + (-t153 - t191) * t202 + t239;
t89 = t154 * t202 - t187 * t194 + t230;
t88 = t214 * t236 + t217 * t225;
t87 = t214 * t225 - t217 * t236;
t86 = t153 * t194 - t154 * t193 + t229;
t85 = -t115 * t186 + t152 * t161 + t227;
t84 = t116 * t186 - t152 * t162 + t224;
t83 = t115 * t162 - t116 * t161 + t223;
t82 = -t107 * t170 + t135 * t136 + t222;
t81 = t108 * t170 - t135 * t137 + t221;
t80 = t107 * t137 - t108 * t136 + t220;
t79 = t120 * t136 + t121 * t129 - t159 * t98 - t170 * t87 + t222;
t78 = -t120 * t137 - t122 * t129 + t159 * t99 + t170 * t88 + t221;
t77 = -t121 * t99 + t122 * t98 - t136 * t88 + t137 * t87 + t220;
t1 = m(1) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(3) * (t86 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + t194 * (t214 * t226 + t217 * t219) / 0.2e1 + t193 * (t214 * t219 - t217 * t226) / 0.2e1 + m(4) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t162 * ((t110 * t248 + t112 * t165 + t114 * t166) * t162 + (t109 * t248 + t111 * t165 + t113 * t166) * t161 + (t143 * t248 + t146 * t165 + t149 * t166) * t186) / 0.2e1 + t161 * ((t110 * t249 + t112 * t163 + t114 * t164) * t162 + (t109 * t249 + t111 * t163 + t113 * t164) * t161 + (t143 * t249 + t146 * t163 + t149 * t164) * t186) / 0.2e1 + t186 * ((-t109 * t161 - t110 * t162 - t143 * t186) * t216 + ((-t112 * t212 + t114 * t215) * t162 + (-t111 * t212 + t113 * t215) * t161 + (-t146 * t212 + t149 * t215) * t186) * t213) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t137 * ((t102 * t248 + t157 * t104 + t158 * t106) * t137 + (t101 * t248 + t103 * t157 + t105 * t158) * t136 + (t131 * t248 + t132 * t157 + t133 * t158) * t170) / 0.2e1 + t136 * ((t102 * t249 + t104 * t155 + t106 * t156) * t137 + (t101 * t249 + t155 * t103 + t156 * t105) * t136 + (t131 * t249 + t132 * t155 + t133 * t156) * t170) / 0.2e1 + t170 * ((-t101 * t136 - t102 * t137 - t131 * t170) * t216 + ((-t104 * t204 + t106 * t205) * t137 + (-t103 * t204 + t105 * t205) * t136 + (-t132 * t204 + t133 * t205) * t170) * t213) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t122 * ((t140 * t95 + t141 * t97 + t248 * t93) * t122 + (t140 * t94 + t141 * t96 + t248 * t92) * t121 + (t126 * t248 + t127 * t140 + t128 * t141) * t159) / 0.2e1 + t121 * ((t138 * t95 + t139 * t97 + t249 * t93) * t122 + (t138 * t94 + t139 * t96 + t249 * t92) * t121 + (t126 * t249 + t127 * t138 + t128 * t139) * t159) / 0.2e1 + t159 * ((-t121 * t92 - t122 * t93 - t126 * t159) * t216 + ((-t199 * t95 + t200 * t97) * t122 + (-t199 * t94 + t200 * t96) * t121 + (-t127 * t199 + t128 * t200) * t159) * t213) / 0.2e1 + ((t148 * t216 + t151 * t213) * t194 + (t147 * t216 + t150 * t213) * t193 + (t180 * t216 + t183 * t213 + Icges(2,3)) * t202) * t202 / 0.2e1 + ((-t181 * t214 + t184 * t217 + Icges(1,4)) * V_base(5) + (-t182 * t214 + t185 * t217 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t181 * t217 + t184 * t214 + Icges(1,2)) * V_base(5) + (t182 * t217 + t185 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t202 * (Icges(2,5) * t217 - Icges(2,6) * t214) + V_base(5) * t202 * (Icges(2,5) * t214 + Icges(2,6) * t217) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
