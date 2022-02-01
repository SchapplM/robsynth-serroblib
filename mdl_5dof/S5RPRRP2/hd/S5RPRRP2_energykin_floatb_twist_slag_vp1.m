% Calculate kinetic energy for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:48
% EndTime: 2022-01-23 09:27:49
% DurationCPUTime: 1.92s
% Computational Cost: add. (1067->198), mult. (750->248), div. (0->0), fcn. (528->8), ass. (0->104)
t223 = Icges(5,4) + Icges(6,4);
t222 = Icges(5,1) + Icges(6,1);
t221 = Icges(5,2) + Icges(6,2);
t152 = cos(qJ(4));
t220 = t223 * t152;
t150 = sin(qJ(4));
t219 = t223 * t150;
t218 = Icges(5,5) + Icges(6,5);
t217 = Icges(5,6) + Icges(6,6);
t216 = -t221 * t150 + t220;
t215 = t222 * t152 - t219;
t214 = rSges(6,1) + pkin(4);
t148 = qJ(1) + pkin(8);
t142 = qJ(3) + t148;
t136 = sin(t142);
t137 = cos(t142);
t213 = t216 * t136 - t217 * t137;
t212 = t217 * t136 + t216 * t137;
t211 = t215 * t136 - t218 * t137;
t210 = t218 * t136 + t215 * t137;
t209 = Icges(5,3) + Icges(6,3);
t208 = t221 * t152 + t219;
t207 = t222 * t150 + t220;
t206 = -t217 * t150 + t218 * t152;
t205 = rSges(6,3) + qJ(5);
t204 = -rSges(6,2) * t150 + t214 * t152;
t109 = -qJD(4) * t137 + V_base(5);
t110 = qJD(4) * t136 + V_base(4);
t143 = V_base(6) + qJD(1);
t138 = qJD(3) + t143;
t201 = (-t208 * t150 + t207 * t152) * t138 + (-t212 * t150 + t210 * t152) * t110 + (-t213 * t150 + t211 * t152) * t109;
t200 = (t218 * t150 + t217 * t152) * t138 + (t209 * t136 + t206 * t137) * t110 + (t206 * t136 - t209 * t137) * t109;
t151 = sin(qJ(1));
t195 = pkin(1) * t151;
t153 = cos(qJ(1));
t194 = pkin(1) * t153;
t140 = sin(t148);
t193 = pkin(2) * t140;
t141 = cos(t148);
t192 = pkin(2) * t141;
t190 = -pkin(5) - qJ(2);
t188 = t204 * t136 - t205 * t137;
t187 = t205 * t136 + t204 * t137;
t186 = Icges(2,4) * t151;
t185 = Icges(3,4) * t140;
t184 = Icges(4,4) * t136;
t179 = -pkin(6) + t190;
t178 = t143 * t194 + V_base(2);
t177 = V_base(5) * pkin(5) + V_base(1);
t174 = rSges(6,2) * t152 + t214 * t150;
t173 = t143 * t192 + t178;
t172 = V_base(5) * qJ(2) + t177;
t171 = V_base(4) * t195 + qJD(2) + V_base(3);
t170 = -t192 - t194;
t169 = V_base(4) * t193 + t171;
t168 = rSges(5,1) * t152 - rSges(5,2) * t150;
t100 = pkin(3) * t137 + pkin(7) * t136;
t158 = t138 * t100 + t179 * V_base(4) + t173;
t157 = V_base(5) * pkin(6) + (-t193 - t195) * t143 + t172;
t99 = pkin(3) * t136 - pkin(7) * t137;
t156 = V_base(4) * t99 + (-t100 + t170) * V_base(5) + t169;
t145 = Icges(2,4) * t153;
t135 = Icges(3,4) * t141;
t133 = Icges(4,4) * t137;
t130 = rSges(2,1) * t153 - t151 * rSges(2,2);
t129 = t151 * rSges(2,1) + rSges(2,2) * t153;
t128 = rSges(5,1) * t150 + rSges(5,2) * t152;
t126 = Icges(2,1) * t153 - t186;
t125 = Icges(2,1) * t151 + t145;
t122 = -Icges(2,2) * t151 + t145;
t121 = Icges(2,2) * t153 + t186;
t113 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t112 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t111 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t108 = rSges(3,1) * t141 - rSges(3,2) * t140;
t107 = rSges(3,1) * t140 + rSges(3,2) * t141;
t106 = Icges(3,1) * t141 - t185;
t105 = Icges(3,1) * t140 + t135;
t104 = -Icges(3,2) * t140 + t135;
t103 = Icges(3,2) * t141 + t185;
t98 = rSges(4,1) * t137 - rSges(4,2) * t136;
t97 = rSges(4,1) * t136 + rSges(4,2) * t137;
t96 = Icges(4,1) * t137 - t184;
t95 = Icges(4,1) * t136 + t133;
t94 = -Icges(4,2) * t136 + t133;
t93 = Icges(4,2) * t137 + t184;
t88 = V_base(5) * rSges(2,3) - t129 * t143 + t177;
t87 = t130 * t143 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t86 = t129 * V_base(4) - t130 * V_base(5) + V_base(3);
t85 = rSges(5,3) * t136 + t137 * t168;
t83 = -rSges(5,3) * t137 + t136 * t168;
t67 = V_base(5) * rSges(3,3) + (-t107 - t195) * t143 + t172;
t66 = t108 * t143 + (-rSges(3,3) + t190) * V_base(4) + t178;
t65 = V_base(4) * t107 + (-t108 - t194) * V_base(5) + t171;
t64 = V_base(5) * rSges(4,3) - t138 * t97 + t157;
t63 = t138 * t98 + (-rSges(4,3) + t179) * V_base(4) + t173;
t62 = V_base(4) * t97 + (t170 - t98) * V_base(5) + t169;
t61 = t109 * t128 + (-t83 - t99) * t138 + t157;
t60 = -t110 * t128 + t138 * t85 + t158;
t59 = -t109 * t85 + t110 * t83 + t156;
t58 = qJD(5) * t136 + t174 * t109 + (-t99 - t188) * t138 + t157;
t57 = -qJD(5) * t137 - t110 * t174 + t138 * t187 + t158;
t56 = -t109 * t187 + t110 * t188 + t156;
t1 = m(1) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + (t201 * t136 - t200 * t137) * t109 / 0.2e1 + (t200 * t136 + t201 * t137) * t110 / 0.2e1 + ((t210 * t150 + t212 * t152) * t110 + (t211 * t150 + t213 * t152) * t109 + (t207 * t150 + t208 * t152 + Icges(4,3)) * t138) * t138 / 0.2e1 + ((-t103 * t140 + t105 * t141 - t151 * t121 + t125 * t153 - t136 * t93 + t137 * t95 + Icges(1,4)) * V_base(5) + (-t140 * t104 + t141 * t106 - t151 * t122 + t153 * t126 - t136 * t94 + t137 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t141 * t103 + t140 * t105 + t153 * t121 + t151 * t125 + t136 * t95 + t137 * t93 + Icges(1,2)) * V_base(5) + (t104 * t141 + t106 * t140 + t122 * t153 + t151 * t126 + t136 * t96 + t137 * t94 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t138 * (Icges(4,5) * t137 - Icges(4,6) * t136) + V_base(5) * t138 * (Icges(4,5) * t136 + Icges(4,6) * t137) + ((Icges(2,5) * t153 + Icges(3,5) * t141 - Icges(2,6) * t151 - Icges(3,6) * t140) * V_base(4) + (Icges(2,5) * t151 + Icges(3,5) * t140 + Icges(2,6) * t153 + Icges(3,6) * t141) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t143) * t143 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
