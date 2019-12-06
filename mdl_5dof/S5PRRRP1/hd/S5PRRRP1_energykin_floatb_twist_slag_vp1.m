% Calculate kinetic energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:43
% EndTime: 2019-12-05 16:39:45
% DurationCPUTime: 2.15s
% Computational Cost: add. (1054->197), mult. (750->251), div. (0->0), fcn. (528->8), ass. (0->104)
t223 = Icges(5,4) + Icges(6,4);
t222 = Icges(5,1) + Icges(6,1);
t221 = Icges(5,2) + Icges(6,2);
t153 = cos(qJ(4));
t220 = t223 * t153;
t152 = sin(qJ(4));
t219 = t223 * t152;
t218 = Icges(5,5) + Icges(6,5);
t217 = Icges(5,6) + Icges(6,6);
t216 = -t221 * t152 + t220;
t215 = t222 * t153 - t219;
t214 = rSges(6,1) + pkin(4);
t148 = pkin(8) + qJ(2);
t143 = qJ(3) + t148;
t136 = sin(t143);
t137 = cos(t143);
t213 = t216 * t136 - t217 * t137;
t212 = t217 * t136 + t216 * t137;
t211 = t215 * t136 - t218 * t137;
t210 = t218 * t136 + t215 * t137;
t209 = Icges(5,3) + Icges(6,3);
t208 = t221 * t153 + t219;
t207 = t222 * t152 + t220;
t206 = -t217 * t152 + t218 * t153;
t205 = rSges(6,3) + qJ(5);
t204 = -rSges(6,2) * t152 + t214 * t153;
t109 = -qJD(4) * t137 + V_base(5);
t110 = qJD(4) * t136 + V_base(4);
t144 = V_base(6) + qJD(2);
t138 = qJD(3) + t144;
t201 = (-t208 * t152 + t207 * t153) * t138 + (-t212 * t152 + t210 * t153) * t110 + (-t213 * t152 + t211 * t153) * t109;
t200 = (t218 * t152 + t217 * t153) * t138 + (t209 * t136 + t206 * t137) * t110 + (t206 * t136 - t209 * t137) * t109;
t150 = cos(pkin(8));
t196 = pkin(1) * t150;
t195 = pkin(2) * t144;
t193 = -pkin(5) - qJ(1);
t191 = t204 * t136 - t205 * t137;
t190 = t205 * t136 + t204 * t137;
t149 = sin(pkin(8));
t189 = Icges(2,4) * t149;
t140 = sin(t148);
t188 = Icges(3,4) * t140;
t187 = Icges(4,4) * t136;
t182 = -pkin(6) + t193;
t175 = pkin(1) * V_base(6);
t181 = t150 * t175 + V_base(2);
t180 = V_base(5) * qJ(1) + V_base(1);
t176 = qJD(1) + V_base(3);
t174 = rSges(6,2) * t153 + t214 * t152;
t141 = cos(t148);
t173 = t141 * t195 + t181;
t172 = V_base(4) * t149 * pkin(1) + t176;
t171 = -pkin(2) * t141 - t196;
t170 = V_base(4) * pkin(2) * t140 + t172;
t169 = rSges(5,1) * t153 - rSges(5,2) * t152;
t159 = V_base(5) * pkin(5) - t149 * t175 + t180;
t100 = pkin(3) * t137 + pkin(7) * t136;
t158 = t138 * t100 + t182 * V_base(4) + t173;
t157 = V_base(5) * pkin(6) - t140 * t195 + t159;
t99 = pkin(3) * t136 - pkin(7) * t137;
t156 = V_base(4) * t99 + (-t100 + t171) * V_base(5) + t170;
t142 = Icges(2,4) * t150;
t135 = Icges(3,4) * t141;
t132 = Icges(4,4) * t137;
t130 = t152 * rSges(5,1) + rSges(5,2) * t153;
t122 = rSges(2,1) * t150 - rSges(2,2) * t149;
t121 = rSges(2,1) * t149 + rSges(2,2) * t150;
t120 = Icges(2,1) * t150 - t189;
t119 = Icges(2,1) * t149 + t142;
t118 = -Icges(2,2) * t149 + t142;
t117 = Icges(2,2) * t150 + t189;
t113 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t112 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t111 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t108 = rSges(3,1) * t141 - rSges(3,2) * t140;
t107 = rSges(3,1) * t140 + rSges(3,2) * t141;
t106 = Icges(3,1) * t141 - t188;
t105 = Icges(3,1) * t140 + t135;
t104 = -Icges(3,2) * t140 + t135;
t103 = Icges(3,2) * t141 + t188;
t98 = rSges(4,1) * t137 - rSges(4,2) * t136;
t97 = rSges(4,1) * t136 + rSges(4,2) * t137;
t96 = Icges(4,1) * t137 - t187;
t95 = Icges(4,1) * t136 + t132;
t94 = -Icges(4,2) * t136 + t132;
t93 = Icges(4,2) * t137 + t187;
t88 = V_base(5) * rSges(2,3) - t121 * V_base(6) + t180;
t87 = t122 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t86 = t136 * rSges(5,3) + t137 * t169;
t84 = -t137 * rSges(5,3) + t136 * t169;
t70 = t121 * V_base(4) - t122 * V_base(5) + t176;
t67 = V_base(5) * rSges(3,3) - t107 * t144 + t159;
t66 = t108 * t144 + (-rSges(3,3) + t193) * V_base(4) + t181;
t65 = t107 * V_base(4) + (-t108 - t196) * V_base(5) + t172;
t64 = V_base(5) * rSges(4,3) - t138 * t97 + t157;
t63 = t138 * t98 + (-rSges(4,3) + t182) * V_base(4) + t173;
t62 = t97 * V_base(4) + (t171 - t98) * V_base(5) + t170;
t61 = t109 * t130 + (-t84 - t99) * t138 + t157;
t60 = -t110 * t130 + t138 * t86 + t158;
t59 = -t109 * t86 + t110 * t84 + t156;
t58 = qJD(5) * t136 + t174 * t109 + (-t99 - t191) * t138 + t157;
t57 = -qJD(5) * t137 - t110 * t174 + t138 * t190 + t158;
t56 = -t109 * t190 + t110 * t191 + t156;
t1 = m(1) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t201 * t136 - t200 * t137) * t109 / 0.2e1 + (t200 * t136 + t201 * t137) * t110 / 0.2e1 + ((t210 * t152 + t212 * t153) * t110 + (t211 * t152 + t213 * t153) * t109 + (t207 * t152 + t208 * t153 + Icges(4,3)) * t138) * t138 / 0.2e1 + ((-t103 * t140 + t105 * t141 - t117 * t149 + t119 * t150 - t136 * t93 + t137 * t95 + Icges(1,4)) * V_base(5) + (-t140 * t104 + t141 * t106 - t149 * t118 + t150 * t120 - t136 * t94 + t137 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t141 * t103 + t140 * t105 + t150 * t117 + t149 * t119 + t136 * t95 + t137 * t93 + Icges(1,2)) * V_base(5) + (t104 * t141 + t106 * t140 + t118 * t150 + t120 * t149 + t136 * t96 + t137 * t94 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t138 * (Icges(4,5) * t137 - Icges(4,6) * t136) + V_base(5) * t138 * (Icges(4,5) * t136 + Icges(4,6) * t137) + ((Icges(2,5) * t149 + Icges(2,6) * t150 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t150 - Icges(2,6) * t149 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t140 + Icges(3,6) * t141) * V_base(5) + (Icges(3,5) * t141 - Icges(3,6) * t140) * V_base(4) + Icges(3,3) * t144 / 0.2e1) * t144;
T = t1;
