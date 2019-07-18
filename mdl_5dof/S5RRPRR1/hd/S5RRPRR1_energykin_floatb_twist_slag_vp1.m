% Calculate kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:38
% EndTime: 2019-07-18 17:20:41
% DurationCPUTime: 2.91s
% Computational Cost: add. (952->248), mult. (1294->367), div. (0->0), fcn. (1158->8), ass. (0->126)
t242 = Icges(3,4) + Icges(4,4);
t241 = Icges(3,1) + Icges(4,1);
t240 = Icges(3,2) + Icges(4,2);
t172 = cos(qJ(2));
t239 = t242 * t172;
t169 = sin(qJ(2));
t238 = t242 * t169;
t237 = Icges(3,5) + Icges(4,5);
t236 = Icges(3,6) + Icges(4,6);
t235 = -t240 * t169 + t239;
t234 = t241 * t172 - t238;
t233 = Icges(3,3) + Icges(4,3);
t170 = sin(qJ(1));
t173 = cos(qJ(1));
t232 = t235 * t170 - t236 * t173;
t231 = t236 * t170 + t235 * t173;
t230 = t234 * t170 - t237 * t173;
t229 = t237 * t170 + t234 * t173;
t228 = t240 * t172 + t238;
t227 = t241 * t169 + t239;
t226 = -t236 * t169 + t237 * t172;
t157 = -qJD(2) * t173 + V_base(5);
t158 = qJD(2) * t170 + V_base(4);
t159 = V_base(6) + qJD(1);
t225 = (-t228 * t169 + t227 * t172) * t159 + (-t231 * t169 + t229 * t172) * t158 + (-t232 * t169 + t230 * t172) * t157;
t224 = (t237 * t169 + t236 * t172) * t159 + (t233 * t170 + t226 * t173) * t158 + (t226 * t170 - t233 * t173) * t157;
t220 = pkin(1) * t172;
t219 = t169 * pkin(1);
t134 = -t173 * qJ(3) + t170 * t220;
t198 = pkin(2) * t172;
t86 = -pkin(3) * t173 + t170 * t198;
t217 = -t134 - t86;
t216 = Icges(2,4) * t170;
t166 = qJ(2) + qJ(4);
t163 = sin(t166);
t211 = Icges(5,4) * t163;
t164 = cos(t166);
t210 = Icges(5,4) * t164;
t209 = t163 * t170;
t208 = t163 * t173;
t168 = sin(qJ(5));
t207 = t168 * t170;
t206 = t168 * t173;
t171 = cos(qJ(5));
t205 = t170 * t171;
t204 = t171 * t173;
t202 = qJD(5) * t163;
t201 = t158 * t134 + V_base(3);
t197 = qJD(3) * t170 + t157 * t219 + V_base(1);
t132 = qJD(4) * t170 + t158;
t140 = pkin(2) * t169;
t196 = t157 * t140 + t197;
t195 = rSges(3,1) * t172 - rSges(3,2) * t169;
t194 = rSges(4,1) * t172 - rSges(4,2) * t169;
t193 = rSges(5,1) * t164 - rSges(5,2) * t163;
t190 = Icges(5,1) * t164 - t211;
t187 = -Icges(5,2) * t163 + t210;
t184 = Icges(5,5) * t164 - Icges(5,6) * t163;
t135 = t170 * qJ(3) + t173 * t220;
t183 = -qJD(3) * t173 + t159 * t135 + V_base(2);
t131 = V_base(5) + (-qJD(2) - qJD(4)) * t173;
t87 = pkin(3) * t170 + t173 * t198;
t182 = t158 * t86 + (-t135 - t87) * t157 + t201;
t181 = (Icges(5,5) * t163 + Icges(5,6) * t164) * t159 + t131 * (-Icges(5,3) * t173 + t170 * t184) + t132 * (Icges(5,3) * t170 + t173 * t184);
t178 = t159 * t87 + (-t140 - t219) * t158 + t183;
t128 = Icges(5,2) * t164 + t211;
t129 = Icges(5,1) * t163 + t210;
t94 = -Icges(5,6) * t173 + t170 * t187;
t95 = Icges(5,6) * t170 + t173 * t187;
t96 = -Icges(5,5) * t173 + t170 * t190;
t97 = Icges(5,5) * t170 + t173 * t190;
t177 = (-t163 * t95 + t164 * t97) * t132 + (-t163 * t94 + t164 * t96) * t131 + (-t128 * t163 + t129 * t164) * t159;
t165 = Icges(2,4) * t173;
t156 = rSges(2,1) * t173 - rSges(2,2) * t170;
t155 = rSges(2,1) * t170 + rSges(2,2) * t173;
t154 = rSges(3,1) * t169 + rSges(3,2) * t172;
t153 = rSges(4,1) * t169 + rSges(4,2) * t172;
t152 = Icges(2,1) * t173 - t216;
t151 = Icges(2,1) * t170 + t165;
t148 = -Icges(2,2) * t170 + t165;
t147 = Icges(2,2) * t173 + t216;
t139 = -qJD(5) * t164 + t159;
t138 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t137 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t136 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(5,1) * t163 + rSges(5,2) * t164;
t126 = t164 * t204 + t207;
t125 = -t164 * t206 + t205;
t124 = t164 * t205 - t206;
t123 = -t164 * t207 - t204;
t121 = rSges(3,3) * t170 + t173 * t195;
t120 = rSges(4,3) * t170 + t173 * t194;
t119 = -rSges(3,3) * t173 + t170 * t195;
t118 = -rSges(4,3) * t173 + t170 * t194;
t103 = V_base(5) * rSges(2,3) - t155 * t159 + V_base(1);
t102 = -V_base(4) * rSges(2,3) + t156 * t159 + V_base(2);
t101 = t173 * t202 + t132;
t100 = t170 * t202 + t131;
t99 = rSges(5,3) * t170 + t173 * t193;
t98 = -rSges(5,3) * t173 + t170 * t193;
t91 = -rSges(6,3) * t164 + (rSges(6,1) * t171 - rSges(6,2) * t168) * t163;
t90 = -Icges(6,5) * t164 + (Icges(6,1) * t171 - Icges(6,4) * t168) * t163;
t89 = -Icges(6,6) * t164 + (Icges(6,4) * t171 - Icges(6,2) * t168) * t163;
t88 = -Icges(6,3) * t164 + (Icges(6,5) * t171 - Icges(6,6) * t168) * t163;
t85 = t155 * V_base(4) - t156 * V_base(5) + V_base(3);
t82 = -t119 * t159 + t154 * t157 + V_base(1);
t81 = t121 * t159 - t154 * t158 + V_base(2);
t80 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t208;
t79 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t209;
t78 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t208;
t77 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t209;
t76 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t208;
t75 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t209;
t74 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t208;
t73 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t209;
t72 = t119 * t158 - t121 * t157 + V_base(3);
t71 = t153 * t157 + (-t118 - t134) * t159 + t197;
t70 = t120 * t159 + (-t153 - t219) * t158 + t183;
t69 = t118 * t158 + (-t120 - t135) * t157 + t201;
t68 = t130 * t131 + (-t98 + t217) * t159 + t196;
t67 = -t130 * t132 + t159 * t99 + t178;
t66 = -t131 * t99 + t132 * t98 + t182;
t65 = -pkin(4) * t131 * t164 + t100 * t91 - t139 * t79 + (-pkin(4) * t209 + t217) * t159 + t196;
t64 = -t101 * t91 + t139 * t80 + (t132 * t164 + t159 * t208) * pkin(4) + t178;
t63 = -t100 * t80 + t101 * t79 + (-t131 * t173 + t132 * t170) * t163 * pkin(4) + t182;
t1 = m(1) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t132 * (t181 * t170 + t177 * t173) / 0.2e1 + t131 * (t177 * t170 - t181 * t173) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t101 * ((t125 * t76 + t126 * t78 + t74 * t208) * t101 + (t125 * t75 + t126 * t77 + t208 * t73) * t100 + (t125 * t89 + t126 * t90 + t208 * t88) * t139) / 0.2e1 + t100 * ((t123 * t76 + t124 * t78 + t209 * t74) * t101 + (t123 * t75 + t124 * t77 + t73 * t209) * t100 + (t123 * t89 + t124 * t90 + t209 * t88) * t139) / 0.2e1 + t139 * ((-t100 * t73 - t101 * t74 - t139 * t88) * t164 + ((-t168 * t76 + t171 * t78) * t101 + (-t168 * t75 + t171 * t77) * t100 + (-t168 * t89 + t171 * t90) * t139) * t163) / 0.2e1 + (t225 * t170 - t224 * t173) * t157 / 0.2e1 + (t224 * t170 + t225 * t173) * t158 / 0.2e1 + ((-t147 * t170 + t151 * t173 + Icges(1,4)) * V_base(5) + (-t148 * t170 + t152 * t173 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t147 * t173 + t151 * t170 + Icges(1,2)) * V_base(5) + (t148 * t173 + t152 * t170 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t163 * t97 + t164 * t95) * t132 + (t163 * t96 + t164 * t94) * t131 + (t229 * t169 + t231 * t172) * t158 + (t230 * t169 + t232 * t172) * t157 + (t128 * t164 + t129 * t163 + t227 * t169 + t228 * t172 + Icges(2,3)) * t159) * t159 / 0.2e1 + t159 * V_base(4) * (Icges(2,5) * t173 - Icges(2,6) * t170) + V_base(5) * t159 * (Icges(2,5) * t170 + Icges(2,6) * t173) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
