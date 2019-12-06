% Calculate kinetic energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:52
% EndTime: 2019-12-05 18:52:54
% DurationCPUTime: 1.74s
% Computational Cost: add. (1247->243), mult. (1060->369), div. (0->0), fcn. (932->10), ass. (0->126)
t158 = sin(qJ(1));
t201 = pkin(1) * t158;
t161 = cos(qJ(1));
t200 = pkin(1) * t161;
t157 = sin(qJ(3));
t199 = pkin(2) * t157;
t160 = cos(qJ(3));
t198 = pkin(2) * t160;
t197 = Icges(2,4) * t158;
t155 = qJ(1) + qJ(2);
t150 = sin(t155);
t196 = Icges(3,4) * t150;
t195 = Icges(4,4) * t157;
t194 = Icges(4,4) * t160;
t154 = qJ(3) + qJ(4);
t149 = sin(t154);
t193 = Icges(5,4) * t149;
t151 = cos(t154);
t192 = Icges(5,4) * t151;
t191 = t149 * t150;
t152 = cos(t155);
t190 = t149 * t152;
t156 = sin(qJ(5));
t189 = t150 * t156;
t159 = cos(qJ(5));
t188 = t150 * t159;
t187 = t152 * t156;
t186 = t152 * t159;
t185 = qJD(5) * t149;
t148 = V_base(6) + qJD(1);
t184 = t148 * t200 + V_base(2);
t183 = V_base(4) * t201 + V_base(3);
t182 = t150 * t198;
t181 = t152 * t198;
t129 = qJD(3) * t150 + V_base(4);
t109 = qJD(4) * t150 + t129;
t147 = qJD(2) + t148;
t178 = rSges(4,1) * t160 - rSges(4,2) * t157;
t177 = rSges(5,1) * t151 - rSges(5,2) * t149;
t176 = -t148 * t201 + V_base(1);
t175 = Icges(4,1) * t160 - t195;
t174 = Icges(5,1) * t151 - t193;
t173 = -Icges(4,2) * t157 + t194;
t172 = -Icges(5,2) * t149 + t192;
t171 = Icges(4,5) * t160 - Icges(4,6) * t157;
t170 = Icges(5,5) * t151 - Icges(5,6) * t149;
t128 = -qJD(3) * t152 + V_base(5);
t169 = t128 * t199 + t176;
t108 = V_base(5) + (-qJD(3) - qJD(4)) * t152;
t168 = -t129 * t199 + t147 * t181 + t184;
t167 = -t200 * V_base(5) + t183;
t166 = (-Icges(5,3) * t152 + t150 * t170) * t108 + (Icges(5,3) * t150 + t152 * t170) * t109 + (Icges(5,5) * t149 + Icges(5,6) * t151) * t147;
t165 = (-Icges(4,3) * t152 + t150 * t171) * t128 + (Icges(4,3) * t150 + t152 * t171) * t129 + (Icges(4,5) * t157 + Icges(4,6) * t160) * t147;
t164 = -t128 * t181 + t129 * t182 + t167;
t115 = Icges(5,2) * t151 + t193;
t118 = Icges(5,1) * t149 + t192;
t80 = -Icges(5,6) * t152 + t150 * t172;
t81 = Icges(5,6) * t150 + t152 * t172;
t82 = -Icges(5,5) * t152 + t150 * t174;
t83 = Icges(5,5) * t150 + t152 * t174;
t163 = (-t149 * t81 + t151 * t83) * t109 + (-t149 * t80 + t151 * t82) * t108 + (-t115 * t149 + t118 * t151) * t147;
t133 = Icges(4,2) * t160 + t195;
t136 = Icges(4,1) * t157 + t194;
t93 = -Icges(4,6) * t152 + t150 * t173;
t94 = Icges(4,6) * t150 + t152 * t173;
t96 = -Icges(4,5) * t152 + t150 * t175;
t97 = Icges(4,5) * t150 + t152 * t175;
t162 = (-t157 * t94 + t160 * t97) * t129 + (-t157 * t93 + t160 * t96) * t128 + (-t133 * t157 + t136 * t160) * t147;
t153 = Icges(2,4) * t161;
t146 = Icges(3,4) * t152;
t141 = rSges(2,1) * t161 - t158 * rSges(2,2);
t140 = t158 * rSges(2,1) + rSges(2,2) * t161;
t139 = rSges(4,1) * t157 + rSges(4,2) * t160;
t138 = Icges(2,1) * t161 - t197;
t137 = Icges(2,1) * t158 + t153;
t135 = -Icges(2,2) * t158 + t153;
t134 = Icges(2,2) * t161 + t197;
t127 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t126 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t125 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t124 = -qJD(5) * t151 + t147;
t123 = rSges(3,1) * t152 - rSges(3,2) * t150;
t122 = rSges(3,1) * t150 + rSges(3,2) * t152;
t121 = rSges(5,1) * t149 + rSges(5,2) * t151;
t120 = Icges(3,1) * t152 - t196;
t119 = Icges(3,1) * t150 + t146;
t117 = -Icges(3,2) * t150 + t146;
t116 = Icges(3,2) * t152 + t196;
t106 = t151 * t186 + t189;
t105 = -t151 * t187 + t188;
t104 = t151 * t188 - t187;
t103 = -t151 * t189 - t186;
t102 = V_base(5) * rSges(2,3) - t140 * t148 + V_base(1);
t101 = -V_base(4) * rSges(2,3) + t141 * t148 + V_base(2);
t100 = rSges(4,3) * t150 + t152 * t178;
t99 = -rSges(4,3) * t152 + t150 * t178;
t98 = -rSges(6,3) * t151 + (rSges(6,1) * t159 - rSges(6,2) * t156) * t149;
t95 = -Icges(6,5) * t151 + (Icges(6,1) * t159 - Icges(6,4) * t156) * t149;
t92 = -Icges(6,6) * t151 + (Icges(6,4) * t159 - Icges(6,2) * t156) * t149;
t89 = -Icges(6,3) * t151 + (Icges(6,5) * t159 - Icges(6,6) * t156) * t149;
t88 = t152 * t185 + t109;
t87 = t150 * t185 + t108;
t86 = t140 * V_base(4) - t141 * V_base(5) + V_base(3);
t85 = rSges(5,3) * t150 + t152 * t177;
t84 = -rSges(5,3) * t152 + t150 * t177;
t77 = V_base(5) * rSges(3,3) - t122 * t147 + t176;
t76 = -V_base(4) * rSges(3,3) + t123 * t147 + t184;
t75 = V_base(4) * t122 + (-t123 - t200) * V_base(5) + t183;
t74 = rSges(6,1) * t106 + rSges(6,2) * t105 + rSges(6,3) * t190;
t73 = rSges(6,1) * t104 + rSges(6,2) * t103 + rSges(6,3) * t191;
t72 = Icges(6,1) * t106 + Icges(6,4) * t105 + Icges(6,5) * t190;
t71 = Icges(6,1) * t104 + Icges(6,4) * t103 + Icges(6,5) * t191;
t70 = Icges(6,4) * t106 + Icges(6,2) * t105 + Icges(6,6) * t190;
t69 = Icges(6,4) * t104 + Icges(6,2) * t103 + Icges(6,6) * t191;
t68 = Icges(6,5) * t106 + Icges(6,6) * t105 + Icges(6,3) * t190;
t67 = Icges(6,5) * t104 + Icges(6,6) * t103 + Icges(6,3) * t191;
t66 = t128 * t139 - t147 * t99 + t176;
t65 = t100 * t147 - t129 * t139 + t184;
t64 = -t128 * t100 + t129 * t99 + t167;
t63 = t108 * t121 + (-t84 - t182) * t147 + t169;
t62 = -t109 * t121 + t147 * t85 + t168;
t61 = -t108 * t85 + t109 * t84 + t164;
t60 = -t124 * t73 - t147 * t182 + t87 * t98 + t169;
t59 = t124 * t74 - t88 * t98 + t168;
t58 = t88 * t73 - t87 * t74 + t164;
t1 = m(1) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t86 ^ 2) / 0.2e1 + m(3) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t129 * (t165 * t150 + t162 * t152) / 0.2e1 + t128 * (t162 * t150 - t165 * t152) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + t109 * (t166 * t150 + t163 * t152) / 0.2e1 + t108 * (t163 * t150 - t166 * t152) / 0.2e1 + m(6) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + t88 * ((t105 * t70 + t106 * t72 + t68 * t190) * t88 + (t105 * t69 + t106 * t71 + t190 * t67) * t87 + (t105 * t92 + t106 * t95 + t190 * t89) * t124) / 0.2e1 + t87 * ((t103 * t70 + t104 * t72 + t191 * t68) * t88 + (t103 * t69 + t104 * t71 + t67 * t191) * t87 + (t103 * t92 + t104 * t95 + t191 * t89) * t124) / 0.2e1 + t124 * ((-t89 * t124 - t67 * t87 - t68 * t88) * t151 + ((-t156 * t70 + t159 * t72) * t88 + (-t156 * t69 + t159 * t71) * t87 + (-t156 * t92 + t159 * t95) * t124) * t149) / 0.2e1 + ((t157 * t97 + t160 * t94) * t129 + (t157 * t96 + t160 * t93) * t128 + (t149 * t83 + t151 * t81) * t109 + (t149 * t82 + t151 * t80) * t108 + (t151 * t115 + t149 * t118 + t160 * t133 + t157 * t136 + Icges(3,3)) * t147) * t147 / 0.2e1 + ((-t116 * t150 + t119 * t152 - t158 * t134 + t137 * t161 + Icges(1,4)) * V_base(5) + (-t117 * t150 + t120 * t152 - t158 * t135 + t138 * t161 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t116 * t152 + t119 * t150 + t134 * t161 + t158 * t137 + Icges(1,2)) * V_base(5) + (t117 * t152 + t120 * t150 + t135 * t161 + t158 * t138 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t147 * (Icges(3,5) * t152 - Icges(3,6) * t150) + V_base(5) * t147 * (Icges(3,5) * t150 + Icges(3,6) * t152) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t158 + Icges(2,6) * t161) * V_base(5) + (Icges(2,5) * t161 - Icges(2,6) * t158) * V_base(4) + Icges(2,3) * t148 / 0.2e1) * t148;
T = t1;
