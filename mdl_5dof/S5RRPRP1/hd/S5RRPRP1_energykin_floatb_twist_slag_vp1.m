% Calculate kinetic energy for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:51
% EndTime: 2020-01-03 11:58:53
% DurationCPUTime: 1.76s
% Computational Cost: add. (1078->199), mult. (750->249), div. (0->0), fcn. (528->8), ass. (0->104)
t226 = Icges(5,4) + Icges(6,4);
t225 = Icges(5,1) + Icges(6,1);
t224 = Icges(5,2) + Icges(6,2);
t149 = cos(qJ(4));
t223 = t226 * t149;
t147 = sin(qJ(4));
t222 = t226 * t147;
t221 = Icges(5,5) + Icges(6,5);
t220 = -Icges(5,6) - Icges(6,6);
t219 = -t224 * t147 + t223;
t218 = t225 * t149 - t222;
t217 = rSges(6,1) + pkin(4);
t145 = qJ(1) + qJ(2);
t137 = pkin(8) + t145;
t133 = sin(t137);
t134 = cos(t137);
t216 = t219 * t133 + t220 * t134;
t215 = t220 * t133 - t219 * t134;
t214 = t218 * t133 - t221 * t134;
t213 = -t221 * t133 - t218 * t134;
t212 = Icges(5,3) + Icges(6,3);
t211 = t224 * t149 + t222;
t210 = t225 * t147 + t223;
t209 = t220 * t147 + t221 * t149;
t208 = -rSges(6,3) - qJ(5);
t207 = -rSges(6,2) * t147 + t217 * t149;
t108 = -qJD(4) * t133 + V_base(6);
t109 = -qJD(4) * t134 + V_base(5);
t138 = V_base(4) + qJD(1);
t135 = qJD(2) + t138;
t206 = t108 * (t215 * t147 - t213 * t149) + t109 * (t216 * t147 - t214 * t149) + t135 * (t211 * t147 - t210 * t149);
t203 = (-t221 * t147 + t220 * t149) * t135 + (-t209 * t133 + t212 * t134) * t109 + (t212 * t133 + t209 * t134) * t108;
t199 = -pkin(5) - pkin(6);
t148 = sin(qJ(1));
t195 = pkin(1) * t148;
t150 = cos(qJ(1));
t194 = pkin(1) * t150;
t139 = sin(t145);
t193 = pkin(2) * t139;
t140 = cos(t145);
t192 = pkin(2) * t140;
t189 = t207 * t133 + t208 * t134;
t188 = t208 * t133 - t207 * t134;
t187 = Icges(2,4) * t150;
t186 = Icges(3,4) * t140;
t185 = Icges(4,4) * t134;
t180 = -qJ(3) + t199;
t179 = t138 * t195 + V_base(3);
t178 = V_base(6) * pkin(5) + V_base(2);
t175 = qJD(3) + V_base(1);
t174 = rSges(6,2) * t149 + t217 * t147;
t173 = t135 * t193 + t179;
t172 = V_base(6) * pkin(6) + t138 * t194 + t178;
t171 = -t193 - t195;
t170 = -t192 - t194;
t169 = rSges(5,1) * t149 - rSges(5,2) * t147;
t155 = V_base(6) * qJ(3) + t135 * t192 + t172;
t96 = pkin(3) * t133 - pkin(7) * t134;
t152 = t135 * t96 + t180 * V_base(5) + t173;
t97 = -pkin(3) * t134 - pkin(7) * t133;
t151 = (t171 - t96) * V_base(6) + t175 + (t97 + t170) * V_base(5);
t142 = Icges(2,4) * t148;
t132 = Icges(3,4) * t139;
t131 = Icges(4,4) * t133;
t128 = -rSges(2,1) * t150 + t148 * rSges(2,2);
t127 = t148 * rSges(2,1) + rSges(2,2) * t150;
t126 = rSges(5,1) * t147 + rSges(5,2) * t149;
t124 = -Icges(2,1) * t150 + t142;
t123 = Icges(2,1) * t148 + t187;
t120 = Icges(2,2) * t148 - t187;
t119 = Icges(2,2) * t150 + t142;
t112 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t111 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t110 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t105 = -rSges(3,1) * t140 + rSges(3,2) * t139;
t104 = rSges(3,1) * t139 + rSges(3,2) * t140;
t103 = -Icges(3,1) * t140 + t132;
t102 = Icges(3,1) * t139 + t186;
t101 = Icges(3,2) * t139 - t186;
t100 = Icges(3,2) * t140 + t132;
t95 = -rSges(4,1) * t134 + rSges(4,2) * t133;
t94 = rSges(4,1) * t133 + rSges(4,2) * t134;
t93 = -Icges(4,1) * t134 + t131;
t92 = Icges(4,1) * t133 + t185;
t91 = Icges(4,2) * t133 - t185;
t90 = Icges(4,2) * t134 + t131;
t85 = V_base(6) * rSges(2,3) - t128 * t138 + t178;
t84 = t127 * t138 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t83 = -t127 * V_base(6) + t128 * V_base(5) + V_base(1);
t82 = -rSges(5,3) * t133 - t134 * t169;
t80 = -rSges(5,3) * t134 + t133 * t169;
t64 = V_base(6) * rSges(3,3) - t105 * t135 + t172;
t63 = t104 * t135 + (-rSges(3,3) + t199) * V_base(5) + t179;
t62 = -V_base(6) * t104 + V_base(5) * t105 + V_base(1) + (-V_base(6) * t148 - t150 * V_base(5)) * pkin(1);
t61 = V_base(6) * rSges(4,3) - t135 * t95 + t155;
t60 = t135 * t94 + (-rSges(4,3) + t180) * V_base(5) + t173;
t59 = (t171 - t94) * V_base(6) + (t170 + t95) * V_base(5) + t175;
t58 = t108 * t126 + (-t82 - t97) * t135 + t155;
t57 = -t109 * t126 + t135 * t80 + t152;
t56 = -t108 * t80 + t109 * t82 + t151;
t55 = -qJD(5) * t134 + t174 * t108 + (-t97 - t188) * t135 + t155;
t54 = -qJD(5) * t133 - t109 * t174 + t135 * t189 + t152;
t53 = -t108 * t189 + t109 * t188 + t151;
t1 = m(1) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(2) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(6) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + (t203 * t133 + t206 * t134) * t108 / 0.2e1 + (-t206 * t133 + t203 * t134) * t109 / 0.2e1 + ((t214 * t147 + t216 * t149) * t109 + (t213 * t147 + t215 * t149) * t108 + (t210 * t147 + t211 * t149 + Icges(3,3) + Icges(4,3)) * t135) * t135 / 0.2e1 + ((t101 * t140 + t103 * t139 + t120 * t150 + t148 * t124 + t133 * t93 + t134 * t91 + Icges(1,6)) * V_base(6) + (t140 * t100 + t139 * t102 + t150 * t119 + t148 * t123 + t133 * t92 + t134 * t90 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t139 * t101 - t140 * t103 + t148 * t120 - t150 * t124 + t133 * t91 - t134 * t93 + Icges(1,3)) * V_base(6) + (t100 * t139 - t102 * t140 + t148 * t119 - t123 * t150 + t133 * t90 - t134 * t92 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + V_base(6) * t135 * (-Icges(3,5) * t140 - Icges(4,5) * t134 + Icges(3,6) * t139 + Icges(4,6) * t133) + V_base(5) * t135 * (Icges(3,5) * t139 + Icges(4,5) * t133 + Icges(3,6) * t140 + Icges(4,6) * t134) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t148 + Icges(2,6) * t150) * V_base(5) + (-Icges(2,5) * t150 + Icges(2,6) * t148) * V_base(6) + Icges(2,3) * t138 / 0.2e1) * t138;
T = t1;
