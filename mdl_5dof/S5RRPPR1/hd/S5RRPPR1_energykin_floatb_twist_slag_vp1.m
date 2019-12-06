% Calculate kinetic energy for
% S5RRPPR1
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:07
% EndTime: 2019-12-05 18:18:09
% DurationCPUTime: 1.43s
% Computational Cost: add. (1114->232), mult. (736->297), div. (0->0), fcn. (514->10), ass. (0->118)
t151 = qJ(1) + qJ(2);
t143 = pkin(8) + t151;
t137 = sin(t143);
t138 = cos(t143);
t145 = sin(t151);
t146 = cos(t151);
t209 = -Icges(3,5) * t145 - Icges(4,5) * t137 - Icges(3,6) * t146 - Icges(4,6) * t138;
t208 = Icges(3,5) * t146 + Icges(4,5) * t138 - Icges(3,6) * t145 - Icges(4,6) * t137;
t152 = sin(pkin(9));
t153 = cos(pkin(9));
t192 = Icges(4,4) * t137;
t189 = Icges(5,4) * t153;
t170 = -Icges(5,2) * t152 + t189;
t82 = Icges(5,6) * t137 + t138 * t170;
t190 = Icges(5,4) * t152;
t172 = Icges(5,1) * t153 - t190;
t84 = Icges(5,5) * t137 + t138 * t172;
t207 = -Icges(4,1) * t138 + t152 * t82 - t153 * t84 + t192;
t191 = Icges(4,4) * t138;
t81 = Icges(5,6) * t138 - t137 * t170;
t83 = Icges(5,5) * t138 - t137 * t172;
t206 = Icges(4,1) * t137 + t152 * t81 - t153 * t83 + t191;
t150 = pkin(9) + qJ(5);
t142 = cos(t150);
t141 = sin(t150);
t188 = Icges(6,4) * t141;
t103 = Icges(6,2) * t142 + t188;
t187 = Icges(6,4) * t142;
t104 = Icges(6,1) * t141 + t187;
t114 = qJD(5) * t137 + V_base(6);
t115 = qJD(5) * t138 + V_base(5);
t144 = V_base(4) + qJD(1);
t139 = qJD(2) + t144;
t169 = -Icges(6,2) * t141 + t187;
t73 = Icges(6,6) * t138 - t137 * t169;
t74 = Icges(6,6) * t137 + t138 * t169;
t171 = Icges(6,1) * t142 - t188;
t75 = Icges(6,5) * t138 - t137 * t171;
t76 = Icges(6,5) * t137 + t138 * t171;
t205 = (t103 * t141 - t104 * t142) * t139 + (t141 * t73 - t142 * t75) * t115 + (t141 * t74 - t142 * t76) * t114;
t204 = -pkin(5) - pkin(6);
t155 = sin(qJ(1));
t202 = pkin(1) * t155;
t156 = cos(qJ(1));
t201 = pkin(1) * t156;
t200 = pkin(2) * t145;
t199 = pkin(2) * t146;
t198 = pkin(4) * t152;
t197 = pkin(4) * t153;
t196 = Icges(2,4) * t155;
t195 = Icges(2,4) * t156;
t194 = Icges(3,4) * t145;
t193 = Icges(3,4) * t146;
t185 = -qJ(3) + t204;
t184 = V_base(6) * pkin(5) + V_base(2);
t100 = pkin(3) * t138 + qJ(4) * t137;
t181 = -t100 - t199;
t180 = V_base(5) * t201 + V_base(6) * t202 + V_base(1);
t179 = rSges(5,1) * t153 - rSges(5,2) * t152;
t178 = rSges(6,1) * t142 - rSges(6,2) * t141;
t177 = -t144 * t202 + V_base(3);
t168 = Icges(5,5) * t153 - Icges(5,6) * t152;
t167 = Icges(6,5) * t142 - Icges(6,6) * t141;
t120 = Icges(5,2) * t153 + t190;
t121 = Icges(5,1) * t152 + t189;
t165 = t120 * t152 - t121 * t153;
t164 = V_base(5) * t199 + V_base(6) * t200 + qJD(3) + t180;
t98 = -pkin(3) * t137 + qJ(4) * t138;
t163 = qJD(4) * t137 + t139 * t98 + t177;
t162 = V_base(6) * pkin(6) - t144 * t201 + t184;
t161 = (Icges(6,5) * t141 + Icges(6,6) * t142) * t139 + t114 * (Icges(6,3) * t137 + t138 * t167) + t115 * (Icges(6,3) * t138 - t137 * t167);
t160 = V_base(5) * t100 + t164;
t159 = V_base(6) * qJ(3) + t162;
t158 = qJD(4) * t138 + t159;
t157 = (Icges(5,5) * t152 + Icges(5,6) * t153) * t139 + (Icges(5,3) * t138 - t137 * t168) * V_base(5) + (Icges(5,3) * t137 + t138 * t168) * V_base(6);
t130 = rSges(2,1) * t156 - t155 * rSges(2,2);
t129 = -t155 * rSges(2,1) - rSges(2,2) * t156;
t128 = Icges(2,1) * t156 - t196;
t127 = -Icges(2,1) * t155 - t195;
t126 = -Icges(2,2) * t155 + t195;
t125 = -Icges(2,2) * t156 - t196;
t122 = rSges(5,1) * t152 + rSges(5,2) * t153;
t118 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t117 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t116 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = rSges(3,1) * t146 - rSges(3,2) * t145;
t112 = -rSges(3,1) * t145 - rSges(3,2) * t146;
t111 = Icges(3,1) * t146 - t194;
t110 = -Icges(3,1) * t145 - t193;
t109 = -Icges(3,2) * t145 + t193;
t108 = -Icges(3,2) * t146 - t194;
t105 = rSges(6,1) * t141 + rSges(6,2) * t142;
t101 = rSges(4,1) * t138 - rSges(4,2) * t137;
t99 = -rSges(4,1) * t137 - rSges(4,2) * t138;
t95 = -Icges(4,2) * t137 + t191;
t94 = -Icges(4,2) * t138 - t192;
t89 = V_base(6) * rSges(2,3) - t130 * t144 + t184;
t88 = t129 * t144 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t129 * V_base(6) + t130 * V_base(5) + V_base(1);
t86 = rSges(5,3) * t137 + t138 * t179;
t85 = rSges(5,3) * t138 - t137 * t179;
t78 = rSges(6,3) * t137 + t138 * t178;
t77 = rSges(6,3) * t138 - t137 * t178;
t70 = pkin(7) * t137 + t138 * t197;
t69 = pkin(7) * t138 - t137 * t197;
t68 = V_base(6) * rSges(3,3) - t139 * t113 + t162;
t67 = t112 * t139 + (-rSges(3,3) + t204) * V_base(5) + t177;
t66 = -t112 * V_base(6) + t113 * V_base(5) + t180;
t65 = V_base(6) * rSges(4,3) + (-t101 - t199) * t139 + t159;
t64 = (t99 - t200) * t139 + (-rSges(4,3) + t185) * V_base(5) + t177;
t63 = t101 * V_base(5) - t99 * V_base(6) + t164;
t62 = V_base(6) * t122 + (t181 - t86) * t139 + t158;
t61 = (t85 - t200) * t139 + (-t122 + t185) * V_base(5) + t163;
t60 = t86 * V_base(5) + (-t85 - t98) * V_base(6) + t160;
t59 = V_base(6) * t198 + t114 * t105 + (t181 - t70 - t78) * t139 + t158;
t58 = -t105 * t115 + (t69 + t77 - t200) * t139 + (t185 - t198) * V_base(5) + t163;
t57 = -t114 * t77 + t115 * t78 + t70 * V_base(5) + (-t69 - t98) * V_base(6) + t160;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t115 * (t205 * t137 + t161 * t138) / 0.2e1 + t114 * (t161 * t137 - t205 * t138) / 0.2e1 + ((t141 * t75 + t142 * t73) * t115 + (t141 * t76 + t142 * t74) * t114 + (t152 * t84 + t153 * t82 + t208) * V_base(6) + (t152 * t83 + t153 * t81 + t209) * V_base(5) + (t103 * t142 + t104 * t141 + t120 * t153 + t121 * t152 + Icges(3,3) + Icges(4,3)) * t139) * t139 / 0.2e1 + (t157 * t138 + (t165 * t137 + t209) * t139 + (-t109 * t146 - t111 * t145 - t126 * t156 - t155 * t128 + t207 * t137 - t138 * t95 + Icges(1,6)) * V_base(6) + (-t108 * t146 - t110 * t145 - t125 * t156 - t155 * t127 + t206 * t137 - t138 * t94 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t157 * t137 + (-t165 * t138 + t208) * t139 + (-t109 * t145 + t111 * t146 - t155 * t126 + t128 * t156 - t137 * t95 - t207 * t138 + Icges(1,3)) * V_base(6) + (-t108 * t145 + t110 * t146 - t155 * t125 + t127 * t156 - t137 * t94 - t206 * t138 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t155 - Icges(2,6) * t156) * V_base(5) + (Icges(2,5) * t156 - Icges(2,6) * t155) * V_base(6) + Icges(2,3) * t144 / 0.2e1) * t144;
T = t1;
