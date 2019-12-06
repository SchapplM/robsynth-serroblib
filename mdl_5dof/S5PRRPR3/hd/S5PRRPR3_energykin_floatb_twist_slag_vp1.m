% Calculate kinetic energy for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:18:59
% EndTime: 2019-12-05 16:19:01
% DurationCPUTime: 1.67s
% Computational Cost: add. (1249->242), mult. (946->326), div. (0->0), fcn. (726->10), ass. (0->126)
t233 = Icges(4,3) + Icges(5,3);
t170 = qJ(3) + pkin(9);
t159 = sin(t170);
t161 = cos(t170);
t174 = sin(qJ(3));
t175 = cos(qJ(3));
t232 = Icges(4,5) * t175 + Icges(5,5) * t161 - Icges(4,6) * t174 - Icges(5,6) * t159;
t169 = pkin(8) + qJ(2);
t158 = sin(t169);
t160 = cos(t169);
t214 = Icges(4,4) * t175;
t193 = -Icges(4,2) * t174 + t214;
t101 = -Icges(4,6) * t160 + t158 * t193;
t102 = Icges(4,6) * t158 + t160 * t193;
t215 = Icges(4,4) * t174;
t196 = Icges(4,1) * t175 - t215;
t103 = -Icges(4,5) * t160 + t158 * t196;
t104 = Icges(4,5) * t158 + t160 * t196;
t213 = Icges(5,4) * t159;
t121 = Icges(5,2) * t161 + t213;
t212 = Icges(5,4) * t161;
t124 = Icges(5,1) * t159 + t212;
t142 = -qJD(3) * t160 + V_base(5);
t143 = qJD(3) * t158 + V_base(4);
t147 = Icges(4,2) * t175 + t215;
t148 = Icges(4,1) * t174 + t214;
t164 = V_base(6) + qJD(2);
t192 = -Icges(5,2) * t159 + t212;
t93 = -Icges(5,6) * t160 + t158 * t192;
t94 = Icges(5,6) * t158 + t160 * t192;
t195 = Icges(5,1) * t161 - t213;
t95 = -Icges(5,5) * t160 + t158 * t195;
t96 = Icges(5,5) * t158 + t160 * t195;
t229 = (-t121 * t159 + t124 * t161 - t147 * t174 + t148 * t175) * t164 + (-t102 * t174 + t104 * t175 - t159 * t94 + t161 * t96) * t143 + (-t101 * t174 + t103 * t175 - t159 * t93 + t161 * t95) * t142;
t228 = (Icges(4,5) * t174 + Icges(5,5) * t159 + Icges(4,6) * t175 + Icges(5,6) * t161) * t164 + (t233 * t158 + t232 * t160) * t143 + (t232 * t158 - t233 * t160) * t142;
t172 = cos(pkin(8));
t224 = pkin(1) * t172;
t223 = pkin(3) * t174;
t222 = pkin(4) * t159;
t221 = t175 * pkin(3);
t220 = -pkin(5) - qJ(1);
t130 = pkin(2) * t158 - pkin(6) * t160;
t80 = -qJ(4) * t160 + t158 * t221;
t218 = -t130 - t80;
t171 = sin(pkin(8));
t217 = Icges(2,4) * t171;
t216 = Icges(3,4) * t158;
t163 = qJ(5) + t170;
t155 = sin(t163);
t211 = Icges(6,4) * t155;
t156 = cos(t163);
t210 = Icges(6,4) * t156;
t209 = pkin(4) * t161;
t201 = pkin(1) * V_base(6);
t207 = t172 * t201 + V_base(2);
t206 = V_base(5) * qJ(1) + V_base(1);
t202 = qJD(1) + V_base(3);
t200 = V_base(4) * t171 * pkin(1) + t202;
t199 = rSges(4,1) * t175 - rSges(4,2) * t174;
t198 = rSges(5,1) * t161 - rSges(5,2) * t159;
t197 = rSges(6,1) * t156 - rSges(6,2) * t155;
t194 = Icges(6,1) * t156 - t211;
t191 = -Icges(6,2) * t155 + t210;
t188 = Icges(6,5) * t156 - Icges(6,6) * t155;
t110 = V_base(5) + (-qJD(3) - qJD(5)) * t160;
t111 = qJD(5) * t158 + t143;
t187 = t110 * (-Icges(6,3) * t160 + t158 * t188) + t111 * (Icges(6,3) * t158 + t160 * t188) + (Icges(6,5) * t155 + Icges(6,6) * t156) * t164;
t184 = V_base(5) * pkin(5) - t171 * t201 + t206;
t131 = pkin(2) * t160 + pkin(6) * t158;
t183 = t164 * t131 + t220 * V_base(4) + t207;
t182 = qJD(4) * t158 + t142 * t223 + t184;
t181 = V_base(4) * t130 + (-t131 - t224) * V_base(5) + t200;
t180 = t143 * t80 + t181;
t81 = qJ(4) * t158 + t160 * t221;
t179 = -qJD(4) * t160 + t164 * t81 + t183;
t114 = Icges(6,2) * t156 + t211;
t115 = Icges(6,1) * t155 + t210;
t84 = -Icges(6,6) * t160 + t158 * t191;
t85 = Icges(6,6) * t158 + t160 * t191;
t86 = -Icges(6,5) * t160 + t158 * t194;
t87 = Icges(6,5) * t158 + t160 * t194;
t178 = (-t155 * t85 + t156 * t87) * t111 + (-t155 * t84 + t156 * t86) * t110 + (-t114 * t155 + t115 * t156) * t164;
t162 = Icges(2,4) * t172;
t154 = Icges(3,4) * t160;
t149 = t174 * rSges(4,1) + rSges(4,2) * t175;
t145 = rSges(2,1) * t172 - rSges(2,2) * t171;
t144 = rSges(2,1) * t171 + rSges(2,2) * t172;
t141 = Icges(2,1) * t172 - t217;
t140 = Icges(2,1) * t171 + t162;
t139 = -Icges(2,2) * t171 + t162;
t138 = Icges(2,2) * t172 + t217;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t160 - rSges(3,2) * t158;
t128 = rSges(5,1) * t159 + rSges(5,2) * t161;
t127 = rSges(3,1) * t158 + rSges(3,2) * t160;
t126 = Icges(3,1) * t160 - t216;
t125 = Icges(3,1) * t158 + t154;
t123 = -Icges(3,2) * t158 + t154;
t122 = Icges(3,2) * t160 + t216;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = V_base(5) * rSges(2,3) - t144 * V_base(6) + t206;
t107 = t145 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t106 = t158 * rSges(4,3) + t160 * t199;
t105 = -t160 * rSges(4,3) + t158 * t199;
t98 = rSges(5,3) * t158 + t160 * t198;
t97 = -rSges(5,3) * t160 + t158 * t198;
t90 = t144 * V_base(4) - t145 * V_base(5) + t202;
t89 = rSges(6,3) * t158 + t160 * t197;
t88 = -rSges(6,3) * t160 + t158 * t197;
t78 = V_base(5) * rSges(3,3) - t127 * t164 + t184;
t77 = t129 * t164 + (-rSges(3,3) + t220) * V_base(4) + t207;
t75 = pkin(7) * t158 + t160 * t209;
t74 = -pkin(7) * t160 + t158 * t209;
t73 = t127 * V_base(4) + (-t129 - t224) * V_base(5) + t200;
t72 = t142 * t149 + (-t105 - t130) * t164 + t184;
t71 = t106 * t164 - t143 * t149 + t183;
t70 = t105 * t143 - t106 * t142 + t181;
t69 = t128 * t142 + (-t97 + t218) * t164 + t182;
t68 = t164 * t98 + (-t128 - t223) * t143 + t179;
t67 = t143 * t97 + (-t81 - t98) * t142 + t180;
t66 = t142 * t222 + t110 * t116 + (-t74 - t88 + t218) * t164 + t182;
t65 = -t111 * t116 + (t75 + t89) * t164 + (-t222 - t223) * t143 + t179;
t64 = -t110 * t89 + t111 * t88 + t143 * t74 + (-t75 - t81) * t142 + t180;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t108 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t111 * (t187 * t158 + t178 * t160) / 0.2e1 + t110 * (t178 * t158 - t187 * t160) / 0.2e1 + (t158 * t229 - t160 * t228) * t142 / 0.2e1 + (t158 * t228 + t229 * t160) * t143 / 0.2e1 + ((-t122 * t158 + t125 * t160 - t138 * t171 + t140 * t172 + Icges(1,4)) * V_base(5) + (-t123 * t158 + t126 * t160 - t139 * t171 + t141 * t172 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t122 * t160 + t125 * t158 + t138 * t172 + t140 * t171 + Icges(1,2)) * V_base(5) + (t123 * t160 + t126 * t158 + t139 * t172 + t141 * t171 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t155 * t87 + t156 * t85) * t111 + (t155 * t86 + t156 * t84) * t110 + (t102 * t175 + t174 * t104 + t159 * t96 + t161 * t94) * t143 + (t101 * t175 + t174 * t103 + t159 * t95 + t161 * t93) * t142 + (t114 * t156 + t115 * t155 + t121 * t161 + t124 * t159 + t147 * t175 + t174 * t148 + Icges(3,3)) * t164) * t164 / 0.2e1 + V_base(4) * t164 * (Icges(3,5) * t160 - Icges(3,6) * t158) + V_base(5) * t164 * (Icges(3,5) * t158 + Icges(3,6) * t160) + ((Icges(2,5) * t171 + Icges(2,6) * t172 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t172 - Icges(2,6) * t171 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
