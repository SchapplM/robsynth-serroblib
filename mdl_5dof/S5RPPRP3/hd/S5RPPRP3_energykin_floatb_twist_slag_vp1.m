% Calculate kinetic energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:41
% EndTime: 2019-12-31 17:50:43
% DurationCPUTime: 1.68s
% Computational Cost: add. (791->188), mult. (754->225), div. (0->0), fcn. (534->6), ass. (0->101)
t232 = Icges(5,4) + Icges(6,4);
t231 = Icges(5,1) + Icges(6,1);
t230 = Icges(5,2) + Icges(6,2);
t141 = cos(qJ(4));
t229 = t232 * t141;
t139 = sin(qJ(4));
t228 = t232 * t139;
t227 = Icges(5,5) + Icges(6,5);
t226 = Icges(5,6) + Icges(6,6);
t225 = t230 * t141 + t228;
t224 = t231 * t139 + t229;
t223 = rSges(6,1) + pkin(4);
t222 = Icges(3,4) + Icges(4,6);
t137 = qJ(1) + pkin(7);
t130 = sin(t137);
t131 = cos(t137);
t221 = t225 * t130 + t226 * t131;
t220 = t226 * t130 - t225 * t131;
t219 = t224 * t130 + t227 * t131;
t218 = t227 * t130 - t224 * t131;
t217 = Icges(3,1) + Icges(4,2);
t216 = -Icges(4,4) + Icges(3,5);
t215 = Icges(4,5) - Icges(3,6);
t214 = Icges(3,2) + Icges(4,3);
t213 = Icges(5,3) + Icges(6,3);
t212 = -t230 * t139 + t229;
t211 = t231 * t141 - t228;
t210 = t227 * t139 + t226 * t141;
t209 = rSges(6,3) + qJ(5);
t208 = t222 * t131;
t207 = rSges(6,2) * t141 + t223 * t139;
t206 = t222 * t130;
t107 = qJD(4) * t130 + V_base(5);
t108 = qJD(4) * t131 + V_base(4);
t132 = V_base(6) + qJD(1);
t205 = (t218 * t139 + t220 * t141) * t107 + (t219 * t139 + t221 * t141) * t108 + (t211 * t139 + t212 * t141) * t132;
t204 = -t214 * t131 - t206;
t203 = t214 * t130 - t208;
t202 = t217 * t130 + t208;
t201 = t217 * t131 - t206;
t200 = (-t226 * t139 + t227 * t141) * t132 + (t210 * t130 + t213 * t131) * t108 + (t213 * t130 - t210 * t131) * t107;
t140 = sin(qJ(1));
t191 = pkin(1) * t140;
t142 = cos(qJ(1));
t190 = pkin(1) * t142;
t188 = t130 * pkin(6);
t187 = t131 * pkin(6);
t186 = -pkin(5) - qJ(2);
t184 = t207 * t130 + t209 * t131;
t183 = t209 * t130 - t207 * t131;
t182 = Icges(2,4) * t140;
t174 = t132 * t190 + V_base(2);
t173 = V_base(5) * pkin(5) + V_base(1);
t97 = t130 * pkin(2) - t131 * qJ(3);
t170 = -t97 - t191;
t100 = t131 * pkin(2) + t130 * qJ(3);
t169 = -t100 - t190;
t168 = -t139 * rSges(6,2) + t223 * t141;
t167 = V_base(5) * qJ(2) + t173;
t166 = V_base(4) * t191 + qJD(2) + V_base(3);
t165 = V_base(4) * t97 + t166;
t164 = qJD(3) * t130 + t167;
t163 = rSges(5,1) * t139 + rSges(5,2) * t141;
t149 = t170 - t188;
t148 = V_base(5) * pkin(3) + t164;
t147 = -qJD(3) * t131 + t132 * t100 + t174;
t144 = t132 * t187 + (-pkin(3) + t186) * V_base(4) + t147;
t143 = V_base(4) * t188 + (t169 - t187) * V_base(5) + t165;
t134 = Icges(2,4) * t142;
t124 = t142 * rSges(2,1) - t140 * rSges(2,2);
t123 = t141 * rSges(5,1) - t139 * rSges(5,2);
t121 = t140 * rSges(2,1) + t142 * rSges(2,2);
t120 = Icges(2,1) * t142 - t182;
t119 = Icges(2,1) * t140 + t134;
t116 = -Icges(2,2) * t140 + t134;
t115 = Icges(2,2) * t142 + t182;
t105 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t104 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t103 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t102 = t131 * rSges(3,1) - t130 * rSges(3,2);
t101 = -t131 * rSges(4,2) + t130 * rSges(4,3);
t99 = t130 * rSges(3,1) + t131 * rSges(3,2);
t98 = -t130 * rSges(4,2) - t131 * rSges(4,3);
t80 = t130 * rSges(5,3) - t163 * t131;
t78 = t131 * rSges(5,3) + t163 * t130;
t76 = V_base(5) * rSges(2,3) - t132 * t121 + t173;
t75 = t132 * t124 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t62 = V_base(4) * t121 - V_base(5) * t124 + V_base(3);
t61 = V_base(5) * rSges(3,3) + (-t99 - t191) * t132 + t167;
t60 = t132 * t102 + (-rSges(3,3) + t186) * V_base(4) + t174;
t59 = V_base(4) * t99 + (-t102 - t190) * V_base(5) + t166;
t58 = V_base(5) * rSges(4,1) + (t170 - t98) * t132 + t164;
t57 = t132 * t101 + (-rSges(4,1) + t186) * V_base(4) + t147;
t56 = V_base(4) * t98 + (-t101 + t169) * V_base(5) + t165;
t55 = t107 * t123 + (t149 - t80) * t132 + t148;
t54 = -t108 * t123 + t132 * t78 + t144;
t53 = -t107 * t78 + t108 * t80 + t143;
t52 = qJD(5) * t131 + t168 * t107 + (t149 - t183) * t132 + t148;
t51 = qJD(5) * t130 - t168 * t108 + t184 * t132 + t144;
t50 = -t184 * t107 + t183 * t108 + t143;
t1 = m(1) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(2) * (t62 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(6) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + (t200 * t130 - t205 * t131) * t107 / 0.2e1 + (t205 * t130 + t200 * t131) * t108 / 0.2e1 + ((-t140 * t115 + t142 * t119 + t204 * t130 + t202 * t131 + Icges(1,4)) * V_base(5) + (-t140 * t116 + t142 * t120 + t203 * t130 + t201 * t131 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t142 * t115 + t140 * t119 + t202 * t130 - t204 * t131 + Icges(1,2)) * V_base(5) + (t142 * t116 + t140 * t120 + t201 * t130 - t203 * t131 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t221 * t139 + t219 * t141) * t108 + (-t220 * t139 + t218 * t141) * t107 + (-t212 * t139 + t211 * t141 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t132) * t132 / 0.2e1 + t132 * V_base(5) * (Icges(2,5) * t140 + Icges(2,6) * t142 + t216 * t130 - t215 * t131) + t132 * V_base(4) * (Icges(2,5) * t142 - Icges(2,6) * t140 + t215 * t130 + t216 * t131) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
