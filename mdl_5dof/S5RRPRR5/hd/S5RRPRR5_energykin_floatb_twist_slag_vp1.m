% Calculate kinetic energy for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:34
% EndTime: 2019-12-05 18:33:36
% DurationCPUTime: 1.75s
% Computational Cost: add. (1258->246), mult. (928->341), div. (0->0), fcn. (708->10), ass. (0->131)
t170 = qJ(1) + qJ(2);
t163 = sin(t170);
t164 = cos(t170);
t171 = sin(pkin(9));
t172 = cos(pkin(9));
t219 = Icges(4,4) * t172;
t195 = -Icges(4,2) * t171 + t219;
t102 = Icges(4,6) * t164 - t195 * t163;
t220 = Icges(4,4) * t171;
t198 = Icges(4,1) * t172 - t220;
t104 = Icges(4,5) * t164 - t198 * t163;
t221 = Icges(3,4) * t164;
t236 = -Icges(3,1) * t163 - t102 * t171 + t104 * t172 - t221;
t103 = Icges(4,6) * t163 + t195 * t164;
t105 = Icges(4,5) * t163 + t198 * t164;
t222 = Icges(3,4) * t163;
t235 = Icges(3,1) * t164 - t103 * t171 + t105 * t172 - t222;
t169 = pkin(9) + qJ(4);
t161 = qJ(5) + t169;
t156 = cos(t161);
t155 = sin(t161);
t216 = Icges(6,4) * t155;
t110 = Icges(6,2) * t156 + t216;
t215 = Icges(6,4) * t156;
t111 = Icges(6,1) * t155 + t215;
t138 = qJD(4) * t163 + V_base(6);
t112 = qJD(5) * t163 + t138;
t139 = qJD(4) * t164 + V_base(5);
t113 = qJD(5) * t164 + t139;
t162 = V_base(4) + qJD(1);
t157 = qJD(2) + t162;
t193 = -Icges(6,2) * t155 + t215;
t83 = Icges(6,6) * t164 - t193 * t163;
t84 = Icges(6,6) * t163 + t193 * t164;
t196 = Icges(6,1) * t156 - t216;
t85 = Icges(6,5) * t164 - t196 * t163;
t86 = Icges(6,5) * t163 + t196 * t164;
t234 = (t110 * t155 - t111 * t156) * t157 + (t155 * t83 - t156 * t85) * t113 + (t155 * t84 - t156 * t86) * t112;
t160 = cos(t169);
t159 = sin(t169);
t218 = Icges(5,4) * t159;
t117 = Icges(5,2) * t160 + t218;
t217 = Icges(5,4) * t160;
t118 = Icges(5,1) * t159 + t217;
t194 = -Icges(5,2) * t159 + t217;
t91 = Icges(5,6) * t164 - t194 * t163;
t92 = Icges(5,6) * t163 + t194 * t164;
t197 = Icges(5,1) * t160 - t218;
t93 = Icges(5,5) * t164 - t197 * t163;
t94 = Icges(5,5) * t163 + t197 * t164;
t233 = (t117 * t159 - t118 * t160) * t157 + (t159 * t91 - t160 * t93) * t139 + (t159 * t92 - t160 * t94) * t138;
t232 = -pkin(5) - pkin(6);
t174 = sin(qJ(1));
t230 = pkin(1) * t174;
t175 = cos(qJ(1));
t229 = pkin(1) * t175;
t228 = pkin(3) * t171;
t227 = pkin(4) * t159;
t226 = t172 * pkin(3);
t128 = pkin(2) * t164 + qJ(3) * t163;
t80 = pkin(7) * t163 + t226 * t164;
t225 = -t128 - t80;
t224 = Icges(2,4) * t174;
t223 = Icges(2,4) * t175;
t213 = pkin(4) * t160;
t211 = V_base(6) * pkin(5) + V_base(2);
t208 = V_base(5) * t229 + V_base(6) * t230 + V_base(1);
t207 = V_base(5) * t128 + t208;
t206 = rSges(4,1) * t172 - rSges(4,2) * t171;
t205 = rSges(5,1) * t160 - rSges(5,2) * t159;
t204 = rSges(6,1) * t156 - rSges(6,2) * t155;
t203 = -t162 * t230 + V_base(3);
t192 = Icges(4,5) * t172 - Icges(4,6) * t171;
t191 = Icges(5,5) * t160 - Icges(5,6) * t159;
t190 = Icges(6,5) * t156 - Icges(6,6) * t155;
t135 = Icges(4,2) * t172 + t220;
t136 = Icges(4,1) * t171 + t219;
t185 = t135 * t171 - t136 * t172;
t126 = -pkin(2) * t163 + qJ(3) * t164;
t184 = qJD(3) * t163 + t157 * t126 + t203;
t183 = V_base(6) * pkin(6) - t162 * t229 + t211;
t182 = (Icges(6,5) * t155 + Icges(6,6) * t156) * t157 + (Icges(6,3) * t163 + t190 * t164) * t112 + (Icges(6,3) * t164 - t190 * t163) * t113;
t181 = (Icges(5,5) * t159 + Icges(5,6) * t160) * t157 + (Icges(5,3) * t163 + t191 * t164) * t138 + (Icges(5,3) * t164 - t191 * t163) * t139;
t180 = qJD(3) * t164 + t183;
t179 = V_base(6) * t228 + t180;
t178 = (Icges(4,3) * t164 - t192 * t163) * V_base(5) + (Icges(4,3) * t163 + t192 * t164) * V_base(6) + (Icges(4,5) * t171 + Icges(4,6) * t172) * t157;
t79 = pkin(7) * t164 - t226 * t163;
t177 = V_base(5) * t80 + (-t126 - t79) * V_base(6) + t207;
t176 = t157 * t79 + (-t228 + t232) * V_base(5) + t184;
t147 = rSges(2,1) * t175 - t174 * rSges(2,2);
t146 = -t174 * rSges(2,1) - rSges(2,2) * t175;
t145 = Icges(2,1) * t175 - t224;
t144 = -Icges(2,1) * t174 - t223;
t143 = -Icges(2,2) * t174 + t223;
t142 = -Icges(2,2) * t175 - t224;
t137 = rSges(4,1) * t171 + rSges(4,2) * t172;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t164 - rSges(3,2) * t163;
t127 = -rSges(3,1) * t163 - rSges(3,2) * t164;
t123 = -Icges(3,2) * t163 + t221;
t122 = -Icges(3,2) * t164 - t222;
t121 = Icges(3,5) * t164 - Icges(3,6) * t163;
t120 = -Icges(3,5) * t163 - Icges(3,6) * t164;
t119 = rSges(5,1) * t159 + rSges(5,2) * t160;
t115 = rSges(6,1) * t155 + rSges(6,2) * t156;
t107 = rSges(4,3) * t163 + t206 * t164;
t106 = rSges(4,3) * t164 - t206 * t163;
t99 = V_base(6) * rSges(2,3) - t147 * t162 + t211;
t98 = t146 * t162 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t146 * V_base(6) + t147 * V_base(5) + V_base(1);
t96 = rSges(5,3) * t163 + t205 * t164;
t95 = rSges(5,3) * t164 - t205 * t163;
t88 = rSges(6,3) * t163 + t204 * t164;
t87 = rSges(6,3) * t164 - t204 * t163;
t76 = V_base(6) * rSges(3,3) - t157 * t129 + t183;
t75 = t127 * t157 + (-rSges(3,3) + t232) * V_base(5) + t203;
t74 = pkin(8) * t163 + t213 * t164;
t73 = pkin(8) * t164 - t213 * t163;
t72 = -t127 * V_base(6) + t129 * V_base(5) + t208;
t71 = V_base(6) * t137 + (-t107 - t128) * t157 + t180;
t70 = t106 * t157 + (-t137 + t232) * V_base(5) + t184;
t69 = t107 * V_base(5) + (-t106 - t126) * V_base(6) + t207;
t68 = t138 * t119 + (-t96 + t225) * t157 + t179;
t67 = -t119 * t139 + t157 * t95 + t176;
t66 = -t138 * t95 + t139 * t96 + t177;
t65 = t138 * t227 + t112 * t115 + (-t74 - t88 + t225) * t157 + t179;
t64 = -t139 * t227 - t113 * t115 + (t73 + t87) * t157 + t176;
t63 = -t112 * t87 + t113 * t88 - t138 * t73 + t139 * t74 + t177;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t139 * (t233 * t163 + t181 * t164) / 0.2e1 + t138 * (t181 * t163 - t233 * t164) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t113 * (t234 * t163 + t182 * t164) / 0.2e1 + t112 * (t182 * t163 - t234 * t164) / 0.2e1 + ((t159 * t93 + t160 * t91) * t139 + (t159 * t94 + t160 * t92) * t138 + (t155 * t85 + t156 * t83) * t113 + (t155 * t86 + t156 * t84) * t112 + (t103 * t172 + t105 * t171 + t121) * V_base(6) + (t102 * t172 + t104 * t171 + t120) * V_base(5) + (t156 * t110 + t155 * t111 + t160 * t117 + t159 * t118 + t172 * t135 + t171 * t136 + Icges(3,3)) * t157) * t157 / 0.2e1 + (t178 * t164 + (t185 * t163 + t120) * t157 + (-t123 * t164 - t143 * t175 - t174 * t145 - t235 * t163 + Icges(1,6)) * V_base(6) + (-t164 * t122 - t175 * t142 - t174 * t144 - t236 * t163 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t178 * t163 + (-t185 * t164 + t121) * t157 + (-t163 * t123 - t174 * t143 + t175 * t145 + t235 * t164 + Icges(1,3)) * V_base(6) + (-t122 * t163 - t174 * t142 + t144 * t175 + t236 * t164 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t174 - Icges(2,6) * t175) * V_base(5) + (Icges(2,5) * t175 - Icges(2,6) * t174) * V_base(6) + Icges(2,3) * t162 / 0.2e1) * t162;
T = t1;
