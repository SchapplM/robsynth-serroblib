% Calculate kinetic energy for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:25
% EndTime: 2019-12-31 20:15:27
% DurationCPUTime: 1.47s
% Computational Cost: add. (921->217), mult. (764->291), div. (0->0), fcn. (544->8), ass. (0->115)
t216 = Icges(3,4) + Icges(4,6);
t215 = Icges(3,1) + Icges(4,2);
t214 = -Icges(4,4) + Icges(3,5);
t213 = Icges(4,5) - Icges(3,6);
t212 = Icges(3,2) + Icges(4,3);
t149 = qJ(1) + qJ(2);
t143 = cos(t149);
t211 = t216 * t143;
t141 = sin(t149);
t210 = t216 * t141;
t208 = t212 * t143 + t210;
t207 = -t212 * t141 + t211;
t206 = t215 * t141 + t211;
t205 = t215 * t143 - t210;
t148 = qJ(4) + qJ(5);
t142 = cos(t148);
t140 = sin(t148);
t188 = Icges(6,4) * t140;
t103 = Icges(6,1) * t142 - t188;
t139 = V_base(6) + qJD(1);
t138 = qJD(2) + t139;
t166 = Icges(6,2) * t142 + t188;
t68 = Icges(6,6) * t143 + t166 * t141;
t69 = Icges(6,6) * t141 - t166 * t143;
t187 = Icges(6,4) * t142;
t168 = Icges(6,1) * t140 + t187;
t70 = Icges(6,5) * t143 + t168 * t141;
t71 = Icges(6,5) * t141 - t168 * t143;
t117 = qJD(4) * t141 + V_base(5);
t88 = qJD(5) * t141 + t117;
t118 = qJD(4) * t143 + V_base(4);
t89 = qJD(5) * t143 + t118;
t98 = -Icges(6,2) * t140 + t187;
t203 = (t140 * t70 + t142 * t68) * t89 + (t140 * t71 + t142 * t69) * t88 + (t103 * t140 + t142 * t98) * t138;
t150 = sin(qJ(4));
t152 = cos(qJ(4));
t189 = Icges(5,4) * t152;
t122 = -Icges(5,2) * t150 + t189;
t190 = Icges(5,4) * t150;
t125 = Icges(5,1) * t152 - t190;
t167 = Icges(5,2) * t152 + t190;
t79 = Icges(5,6) * t143 + t167 * t141;
t80 = Icges(5,6) * t141 - t167 * t143;
t169 = Icges(5,1) * t150 + t189;
t81 = Icges(5,5) * t143 + t169 * t141;
t82 = Icges(5,5) * t141 - t169 * t143;
t202 = (t150 * t81 + t152 * t79) * t118 + (t150 * t82 + t152 * t80) * t117 + (t122 * t152 + t125 * t150) * t138;
t201 = -pkin(5) - pkin(6);
t151 = sin(qJ(1));
t199 = pkin(1) * t151;
t153 = cos(qJ(1));
t198 = pkin(1) * t153;
t197 = pkin(4) * t150;
t196 = pkin(4) * t152;
t195 = t141 * pkin(7);
t194 = t143 * pkin(7);
t192 = Icges(2,4) * t151;
t184 = t139 * t198 + V_base(2);
t183 = V_base(4) * t199 + V_base(3);
t182 = V_base(5) * pkin(5) + V_base(1);
t110 = pkin(2) * t143 + qJ(3) * t141;
t179 = -t110 - t198;
t106 = pkin(2) * t141 - qJ(3) * t143;
t178 = -t106 - t195;
t177 = V_base(4) * t106 + t183;
t176 = rSges(5,1) * t150 + rSges(5,2) * t152;
t175 = rSges(6,1) * t140 + rSges(6,2) * t142;
t165 = Icges(5,5) * t150 + Icges(5,6) * t152;
t164 = Icges(6,5) * t140 + Icges(6,6) * t142;
t162 = -qJD(3) * t143 + t138 * t110 + t184;
t161 = V_base(5) * pkin(6) - t139 * t199 + t182;
t160 = (Icges(6,5) * t142 - Icges(6,6) * t140) * t138 + (Icges(6,3) * t143 + t164 * t141) * t89 + (Icges(6,3) * t141 - t164 * t143) * t88;
t159 = (Icges(5,3) * t141 - t165 * t143) * t117 + (Icges(5,3) * t143 + t165 * t141) * t118 + (Icges(5,5) * t152 - Icges(5,6) * t150) * t138;
t158 = qJD(3) * t141 + t161;
t157 = V_base(5) * pkin(3) + t158;
t156 = t138 * t194 + (-pkin(3) + t201) * V_base(4) + t162;
t155 = V_base(4) * t195 + (t179 - t194) * V_base(5) + t177;
t144 = Icges(2,4) * t153;
t130 = rSges(2,1) * t153 - rSges(2,2) * t151;
t129 = rSges(5,1) * t152 - rSges(5,2) * t150;
t128 = rSges(2,1) * t151 + rSges(2,2) * t153;
t127 = Icges(2,1) * t153 - t192;
t126 = Icges(2,1) * t151 + t144;
t124 = -Icges(2,2) * t151 + t144;
t123 = Icges(2,2) * t153 + t192;
t116 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t115 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t114 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t112 = rSges(3,1) * t143 - rSges(3,2) * t141;
t111 = -rSges(4,2) * t143 + rSges(4,3) * t141;
t109 = rSges(6,1) * t142 - rSges(6,2) * t140;
t108 = rSges(3,1) * t141 + rSges(3,2) * t143;
t107 = -rSges(4,2) * t141 - rSges(4,3) * t143;
t86 = pkin(8) * t143 + t141 * t197;
t85 = pkin(8) * t141 - t143 * t197;
t84 = rSges(5,3) * t141 - t176 * t143;
t83 = rSges(5,3) * t143 + t176 * t141;
t76 = V_base(5) * rSges(2,3) - t128 * t139 + t182;
t75 = t130 * t139 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t74 = t128 * V_base(4) - t130 * V_base(5) + V_base(3);
t73 = rSges(6,3) * t141 - t175 * t143;
t72 = rSges(6,3) * t143 + t175 * t141;
t65 = V_base(5) * rSges(3,3) - t108 * t138 + t161;
t64 = t112 * t138 + (-rSges(3,3) + t201) * V_base(4) + t184;
t63 = t108 * V_base(4) + (-t112 - t198) * V_base(5) + t183;
t62 = V_base(5) * rSges(4,1) + (-t106 - t107) * t138 + t158;
t61 = t111 * t138 + (-rSges(4,1) + t201) * V_base(4) + t162;
t60 = t107 * V_base(4) + (-t111 + t179) * V_base(5) + t177;
t59 = t117 * t129 + (t178 - t84) * t138 + t157;
t58 = -t118 * t129 + t138 * t83 + t156;
t57 = -t117 * t83 + t118 * t84 + t155;
t56 = t117 * t196 + t109 * t88 + (t178 - t73 - t85) * t138 + t157;
t55 = -t118 * t196 - t109 * t89 + (t72 + t86) * t138 + t156;
t54 = -t117 * t86 + t118 * t85 - t72 * t88 + t73 * t89 + t155;
t1 = m(1) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t118 * (t202 * t141 + t159 * t143) / 0.2e1 + t117 * (t159 * t141 - t202 * t143) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + t89 * (t203 * t141 + t160 * t143) / 0.2e1 + t88 * (t160 * t141 - t203 * t143) / 0.2e1 + ((-t150 * t79 + t152 * t81) * t118 + (-t150 * t80 + t152 * t82) * t117 + (-t140 * t68 + t142 * t70) * t89 + (-t140 * t69 + t142 * t71) * t88 + (t103 * t142 - t122 * t150 + t125 * t152 - t140 * t98 + Icges(4,1) + Icges(3,3)) * t138) * t138 / 0.2e1 + ((-t123 * t151 + t126 * t153 - t141 * t208 + t143 * t206 + Icges(1,4)) * V_base(5) + (-t124 * t151 + t127 * t153 - t207 * t141 + t205 * t143 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t123 * t153 + t126 * t151 + t206 * t141 + t208 * t143 + Icges(1,2)) * V_base(5) + (t124 * t153 + t127 * t151 + t141 * t205 + t143 * t207 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t138 * (t214 * t141 - t213 * t143) + V_base(4) * t138 * (t213 * t141 + t214 * t143) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t151 + Icges(2,6) * t153) * V_base(5) + (Icges(2,5) * t153 - Icges(2,6) * t151) * V_base(4) + Icges(2,3) * t139 / 0.2e1) * t139;
T = t1;
