% Calculate kinetic energy for
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:29
% EndTime: 2019-12-31 18:05:31
% DurationCPUTime: 1.84s
% Computational Cost: add. (566->225), mult. (988->308), div. (0->0), fcn. (826->6), ass. (0->112)
t208 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t207 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t206 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t205 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t204 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t151 = sin(qJ(1));
t203 = t208 * t151;
t154 = cos(qJ(1));
t202 = t208 * t154;
t201 = t154 * t207 - t203;
t200 = t151 * t207 + t202;
t199 = -t154 * t204 - t203;
t198 = t151 * t204 - t202;
t195 = -pkin(2) - pkin(5);
t193 = pkin(6) * t151;
t192 = pkin(6) * t154;
t150 = sin(qJ(4));
t190 = Icges(5,4) * t150;
t153 = cos(qJ(4));
t189 = Icges(5,4) * t153;
t185 = qJ(3) * t151;
t184 = qJ(3) * t154;
t149 = sin(qJ(5));
t183 = t149 * t154;
t182 = t151 * t149;
t152 = cos(qJ(5));
t181 = t151 * t152;
t180 = t151 * t153;
t179 = t152 * t154;
t178 = t153 * t154;
t177 = qJD(5) * t153;
t123 = pkin(1) * t151 - qJ(2) * t154;
t176 = t123 * V_base(4) + V_base(3);
t175 = V_base(5) * pkin(5) + V_base(1);
t134 = qJD(4) * t154 + V_base(5);
t139 = V_base(6) + qJD(1);
t172 = t185 * V_base(4) + t176;
t171 = qJD(2) * t151 + t175;
t170 = -t123 - t185;
t128 = pkin(1) * t154 + qJ(2) * t151;
t169 = -t128 - t184;
t168 = pkin(4) * t150 - pkin(7) * t153;
t135 = -qJD(4) * t151 + V_base(4);
t167 = rSges(5,1) * t150 + rSges(5,2) * t153;
t166 = Icges(5,1) * t150 + t189;
t165 = Icges(5,2) * t153 + t190;
t164 = Icges(5,5) * t150 + Icges(5,6) * t153;
t163 = -qJD(2) * t154 + t128 * t139 + V_base(2);
t162 = V_base(5) * pkin(2) + qJD(3) * t154 + t171;
t161 = t170 - t192;
t160 = V_base(5) * pkin(3) + t162;
t159 = qJD(3) * t151 + t139 * t184 + t163;
t158 = (Icges(5,5) * t153 - Icges(5,6) * t150) * t139 + (Icges(5,3) * t154 + t151 * t164) * t134 + (-Icges(5,3) * t151 + t154 * t164) * t135;
t157 = V_base(4) * t192 + t172 + (t169 + t193) * V_base(5);
t156 = (-pkin(3) + t195) * V_base(4) + t159;
t112 = -Icges(5,2) * t150 + t189;
t119 = Icges(5,1) * t153 - t190;
t80 = Icges(5,6) * t154 + t151 * t165;
t81 = -Icges(5,6) * t151 + t154 * t165;
t83 = Icges(5,5) * t154 + t151 * t166;
t84 = -Icges(5,5) * t151 + t154 * t166;
t155 = (t150 * t84 + t153 * t81) * t135 + (t150 * t83 + t153 * t80) * t134 + (t112 * t153 + t119 * t150) * t139;
t133 = pkin(4) * t153 + pkin(7) * t150;
t131 = rSges(2,1) * t154 - rSges(2,2) * t151;
t130 = -rSges(3,2) * t154 + rSges(3,3) * t151;
t129 = -rSges(4,2) * t154 + rSges(4,3) * t151;
t127 = rSges(5,1) * t153 - rSges(5,2) * t150;
t126 = rSges(2,1) * t151 + rSges(2,2) * t154;
t125 = -rSges(3,2) * t151 - rSges(3,3) * t154;
t124 = rSges(4,2) * t151 + rSges(4,3) * t154;
t122 = qJD(5) * t150 + t139;
t100 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t99 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t98 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t96 = t168 * t154;
t95 = t168 * t151;
t93 = t150 * t179 - t182;
t92 = -t150 * t183 - t181;
t91 = t150 * t181 + t183;
t90 = -t150 * t182 + t179;
t89 = -t154 * t177 + t135;
t88 = -t151 * t177 + t134;
t87 = -rSges(5,3) * t151 + t154 * t167;
t86 = rSges(6,3) * t150 + (rSges(6,1) * t152 - rSges(6,2) * t149) * t153;
t85 = rSges(5,3) * t154 + t151 * t167;
t82 = Icges(6,5) * t150 + (Icges(6,1) * t152 - Icges(6,4) * t149) * t153;
t79 = Icges(6,6) * t150 + (Icges(6,4) * t152 - Icges(6,2) * t149) * t153;
t76 = Icges(6,3) * t150 + (Icges(6,5) * t152 - Icges(6,6) * t149) * t153;
t75 = V_base(5) * rSges(2,3) - t126 * t139 + t175;
t74 = t131 * t139 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t73 = t126 * V_base(4) - t131 * V_base(5) + V_base(3);
t72 = rSges(6,1) * t93 + rSges(6,2) * t92 - rSges(6,3) * t178;
t71 = rSges(6,1) * t91 + rSges(6,2) * t90 - rSges(6,3) * t180;
t70 = Icges(6,1) * t93 + Icges(6,4) * t92 - Icges(6,5) * t178;
t69 = Icges(6,1) * t91 + Icges(6,4) * t90 - Icges(6,5) * t180;
t68 = Icges(6,4) * t93 + Icges(6,2) * t92 - Icges(6,6) * t178;
t67 = Icges(6,4) * t91 + Icges(6,2) * t90 - Icges(6,6) * t180;
t66 = Icges(6,5) * t93 + Icges(6,6) * t92 - Icges(6,3) * t178;
t65 = Icges(6,5) * t91 + Icges(6,6) * t90 - Icges(6,3) * t180;
t64 = V_base(5) * rSges(3,1) + (-t123 - t125) * t139 + t171;
t63 = t139 * t130 + (-rSges(3,1) - pkin(5)) * V_base(4) + t163;
t62 = t125 * V_base(4) + (-t128 - t130) * V_base(5) + t176;
t61 = V_base(5) * rSges(4,1) + (-t129 + t170) * t139 + t162;
t60 = t139 * t124 + (-rSges(4,1) + t195) * V_base(4) + t159;
t59 = V_base(4) * t129 + (-t124 + t169) * V_base(5) + t172;
t58 = t134 * t127 + (t161 - t85) * t139 + t160;
t57 = -t135 * t127 + (t87 - t193) * t139 + t156;
t56 = -t134 * t87 + t135 * t85 + t157;
t55 = -t122 * t71 + t134 * t133 + t88 * t86 + (t161 - t95) * t139 + t160;
t54 = t122 * t72 - t135 * t133 - t89 * t86 + (t96 - t193) * t139 + t156;
t53 = -t134 * t96 + t135 * t95 + t89 * t71 - t88 * t72 + t157;
t1 = m(1) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(2) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + t135 * (-t158 * t151 + t155 * t154) / 0.2e1 + t134 * (t155 * t151 + t158 * t154) / 0.2e1 + m(6) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + t89 * ((-t66 * t178 + t92 * t68 + t93 * t70) * t89 + (-t178 * t65 + t67 * t92 + t69 * t93) * t88 + (-t178 * t76 + t79 * t92 + t82 * t93) * t122) / 0.2e1 + t88 * ((-t180 * t66 + t68 * t90 + t70 * t91) * t89 + (-t65 * t180 + t90 * t67 + t91 * t69) * t88 + (-t180 * t76 + t79 * t90 + t82 * t91) * t122) / 0.2e1 + t122 * ((t122 * t76 + t65 * t88 + t66 * t89) * t150 + ((-t149 * t68 + t152 * t70) * t89 + (-t149 * t67 + t152 * t69) * t88 + (-t149 * t79 + t152 * t82) * t122) * t153) / 0.2e1 + ((-t150 * t81 + t153 * t84) * t135 + (-t150 * t80 + t153 * t83) * t134 + (-t150 * t112 + t153 * t119 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t139) * t139 / 0.2e1 + ((t151 * t199 + t154 * t200 + Icges(1,4)) * V_base(5) + (t198 * t151 + t154 * t201 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t151 * t200 - t154 * t199 + Icges(1,2)) * V_base(5) + (t151 * t201 - t154 * t198 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t139 * (t151 * t206 - t154 * t205) + V_base(4) * t139 * (t151 * t205 + t154 * t206) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
