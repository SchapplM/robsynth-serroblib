% Calculate kinetic energy for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:53
% EndTime: 2019-12-31 19:00:55
% DurationCPUTime: 1.68s
% Computational Cost: add. (1409->267), mult. (1172->387), div. (0->0), fcn. (1008->10), ass. (0->134)
t181 = sin(qJ(1));
t228 = pkin(1) * t181;
t184 = cos(qJ(1));
t227 = pkin(1) * t184;
t180 = sin(qJ(3));
t226 = pkin(3) * t180;
t183 = cos(qJ(3));
t225 = pkin(3) * t183;
t224 = -pkin(5) - qJ(2);
t222 = Icges(2,4) * t181;
t177 = qJ(1) + pkin(9);
t169 = sin(t177);
t221 = Icges(3,4) * t169;
t220 = Icges(4,4) * t180;
t219 = Icges(4,4) * t183;
t178 = qJ(3) + qJ(4);
t172 = sin(t178);
t218 = Icges(5,4) * t172;
t173 = cos(t178);
t217 = Icges(5,4) * t173;
t216 = t169 * t172;
t170 = cos(t177);
t215 = t170 * t172;
t179 = sin(qJ(5));
t214 = t173 * t179;
t182 = cos(qJ(5));
t213 = t173 * t182;
t212 = qJD(5) * t172;
t171 = V_base(6) + qJD(1);
t211 = t171 * t227 + V_base(2);
t210 = V_base(5) * pkin(5) + V_base(1);
t150 = qJD(3) * t169 + V_base(4);
t138 = t169 * pkin(2) - t170 * pkin(6);
t207 = -t138 - t228;
t206 = V_base(5) * qJ(2) + t210;
t205 = V_base(4) * t228 + qJD(2) + V_base(3);
t127 = qJD(4) * t169 + t150;
t92 = -pkin(7) * t170 + t169 * t225;
t204 = t207 - t92;
t149 = -qJD(3) * t170 + V_base(5);
t203 = t149 * t226 + t206;
t202 = pkin(4) * t173 + pkin(8) * t172;
t201 = rSges(4,1) * t183 - rSges(4,2) * t180;
t200 = rSges(5,1) * t173 - rSges(5,2) * t172;
t199 = Icges(4,1) * t183 - t220;
t198 = Icges(5,1) * t173 - t218;
t197 = -Icges(4,2) * t180 + t219;
t196 = -Icges(5,2) * t172 + t217;
t195 = Icges(4,5) * t183 - Icges(4,6) * t180;
t194 = Icges(5,5) * t173 - Icges(5,6) * t172;
t126 = V_base(5) + (-qJD(3) - qJD(4)) * t170;
t193 = (-Icges(5,3) * t170 + t169 * t194) * t126 + (Icges(5,3) * t169 + t170 * t194) * t127 + (Icges(5,5) * t172 + Icges(5,6) * t173) * t171;
t192 = (-Icges(4,3) * t170 + t169 * t195) * t149 + (Icges(4,3) * t169 + t170 * t195) * t150 + (Icges(4,5) * t180 + Icges(4,6) * t183) * t171;
t139 = t170 * pkin(2) + t169 * pkin(6);
t191 = t171 * t139 + t224 * V_base(4) + t211;
t190 = V_base(4) * t138 + (-t139 - t227) * V_base(5) + t205;
t93 = pkin(7) * t169 + t170 * t225;
t189 = -t150 * t226 + t171 * t93 + t191;
t188 = -t149 * t93 + t150 * t92 + t190;
t141 = Icges(5,2) * t173 + t218;
t142 = Icges(5,1) * t172 + t217;
t96 = -Icges(5,6) * t170 + t169 * t196;
t97 = Icges(5,6) * t169 + t170 * t196;
t98 = -Icges(5,5) * t170 + t169 * t198;
t99 = Icges(5,5) * t169 + t170 * t198;
t187 = (-t172 * t97 + t173 * t99) * t127 + (-t172 * t96 + t173 * t98) * t126 + (-t141 * t172 + t142 * t173) * t171;
t107 = -Icges(4,6) * t170 + t169 * t197;
t108 = Icges(4,6) * t169 + t170 * t197;
t109 = -Icges(4,5) * t170 + t169 * t199;
t110 = Icges(4,5) * t169 + t170 * t199;
t154 = Icges(4,2) * t183 + t220;
t157 = Icges(4,1) * t180 + t219;
t186 = (-t108 * t180 + t110 * t183) * t150 + (-t107 * t180 + t109 * t183) * t149 + (-t154 * t180 + t157 * t183) * t171;
t175 = Icges(2,4) * t184;
t167 = Icges(3,4) * t170;
t162 = rSges(2,1) * t184 - rSges(2,2) * t181;
t161 = rSges(2,1) * t181 + rSges(2,2) * t184;
t160 = rSges(4,1) * t180 + rSges(4,2) * t183;
t159 = Icges(2,1) * t184 - t222;
t158 = Icges(2,1) * t181 + t175;
t156 = -Icges(2,2) * t181 + t175;
t155 = Icges(2,2) * t184 + t222;
t148 = -qJD(5) * t173 + t171;
t147 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t146 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t145 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t144 = pkin(4) * t172 - pkin(8) * t173;
t143 = rSges(5,1) * t172 + rSges(5,2) * t173;
t137 = rSges(3,1) * t170 - rSges(3,2) * t169;
t136 = rSges(3,1) * t169 + rSges(3,2) * t170;
t135 = Icges(3,1) * t170 - t221;
t134 = Icges(3,1) * t169 + t167;
t133 = -Icges(3,2) * t169 + t167;
t132 = Icges(3,2) * t170 + t221;
t124 = t169 * t179 + t170 * t213;
t123 = t169 * t182 - t170 * t214;
t122 = t169 * t213 - t170 * t179;
t121 = -t169 * t214 - t170 * t182;
t120 = t202 * t170;
t119 = t202 * t169;
t118 = -rSges(6,3) * t173 + (rSges(6,1) * t182 - rSges(6,2) * t179) * t172;
t117 = -Icges(6,5) * t173 + (Icges(6,1) * t182 - Icges(6,4) * t179) * t172;
t116 = -Icges(6,6) * t173 + (Icges(6,4) * t182 - Icges(6,2) * t179) * t172;
t115 = -Icges(6,3) * t173 + (Icges(6,5) * t182 - Icges(6,6) * t179) * t172;
t114 = rSges(4,3) * t169 + t170 * t201;
t113 = -rSges(4,3) * t170 + t169 * t201;
t112 = V_base(5) * rSges(2,3) - t161 * t171 + t210;
t111 = t162 * t171 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t104 = t161 * V_base(4) - t162 * V_base(5) + V_base(3);
t103 = t170 * t212 + t127;
t102 = t169 * t212 + t126;
t101 = rSges(5,3) * t169 + t170 * t200;
t100 = -rSges(5,3) * t170 + t169 * t200;
t89 = V_base(5) * rSges(3,3) + (-t136 - t228) * t171 + t206;
t88 = t137 * t171 + (-rSges(3,3) + t224) * V_base(4) + t211;
t87 = t136 * V_base(4) + (-t137 - t227) * V_base(5) + t205;
t86 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t215;
t85 = rSges(6,1) * t122 + rSges(6,2) * t121 + rSges(6,3) * t216;
t84 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t215;
t83 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t216;
t82 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t215;
t81 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t216;
t80 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t215;
t79 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t216;
t78 = t149 * t160 + (-t113 + t207) * t171 + t206;
t77 = t114 * t171 - t150 * t160 + t191;
t76 = t113 * t150 - t114 * t149 + t190;
t75 = t126 * t143 + (-t100 + t204) * t171 + t203;
t74 = t101 * t171 - t127 * t143 + t189;
t73 = t100 * t127 - t101 * t126 + t188;
t72 = t102 * t118 + t126 * t144 - t148 * t85 + (-t119 + t204) * t171 + t203;
t71 = -t103 * t118 + t120 * t171 - t127 * t144 + t148 * t86 + t189;
t70 = -t102 * t86 + t103 * t85 + t119 * t127 - t120 * t126 + t188;
t1 = m(1) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + t150 * (t192 * t169 + t186 * t170) / 0.2e1 + t149 * (t186 * t169 - t192 * t170) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t127 * (t193 * t169 + t187 * t170) / 0.2e1 + t126 * (t187 * t169 - t193 * t170) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t103 * ((t123 * t82 + t124 * t84 + t80 * t215) * t103 + (t123 * t81 + t124 * t83 + t215 * t79) * t102 + (t115 * t215 + t116 * t123 + t117 * t124) * t148) / 0.2e1 + t102 * ((t121 * t82 + t122 * t84 + t216 * t80) * t103 + (t121 * t81 + t122 * t83 + t79 * t216) * t102 + (t115 * t216 + t116 * t121 + t117 * t122) * t148) / 0.2e1 + t148 * ((-t79 * t102 - t80 * t103 - t115 * t148) * t173 + ((-t179 * t82 + t182 * t84) * t103 + (-t179 * t81 + t182 * t83) * t102 + (-t116 * t179 + t117 * t182) * t148) * t172) / 0.2e1 + ((-t132 * t169 + t134 * t170 - t155 * t181 + t158 * t184 + Icges(1,4)) * V_base(5) + (-t133 * t169 + t135 * t170 - t156 * t181 + t159 * t184 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t132 * t170 + t134 * t169 + t155 * t184 + t158 * t181 + Icges(1,2)) * V_base(5) + (t133 * t170 + t135 * t169 + t156 * t184 + t159 * t181 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t108 * t183 + t110 * t180) * t150 + (t107 * t183 + t109 * t180) * t149 + (t172 * t99 + t173 * t97) * t127 + (t172 * t98 + t173 * t96) * t126 + (t141 * t173 + t142 * t172 + t154 * t183 + t157 * t180 + Icges(2,3) + Icges(3,3)) * t171) * t171 / 0.2e1 + t171 * V_base(5) * (Icges(2,5) * t181 + Icges(3,5) * t169 + Icges(2,6) * t184 + Icges(3,6) * t170) + t171 * V_base(4) * (Icges(2,5) * t184 + Icges(3,5) * t170 - Icges(2,6) * t181 - Icges(3,6) * t169) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
