% Calculate joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:02
% EndTime: 2019-12-31 19:54:07
% DurationCPUTime: 2.26s
% Computational Cost: add. (2362->214), mult. (2184->309), div. (0->0), fcn. (1965->8), ass. (0->114)
t220 = Icges(5,4) - Icges(6,5);
t219 = Icges(5,1) + Icges(6,1);
t218 = -Icges(5,2) - Icges(6,3);
t115 = qJ(2) + pkin(8);
t107 = qJ(4) + t115;
t103 = cos(t107);
t217 = t220 * t103;
t102 = sin(t107);
t216 = t220 * t102;
t215 = Icges(6,4) + Icges(5,5);
t214 = Icges(5,6) - Icges(6,6);
t213 = t218 * t102 + t217;
t212 = t219 * t103 - t216;
t120 = sin(qJ(1));
t122 = cos(qJ(1));
t211 = t213 * t120 - t214 * t122;
t210 = t214 * t120 + t213 * t122;
t209 = t212 * t120 - t215 * t122;
t208 = t215 * t120 + t212 * t122;
t207 = Icges(6,2) + Icges(5,3);
t206 = -t214 * t102 + t215 * t103;
t190 = rSges(6,3) + qJ(5);
t201 = rSges(6,1) + pkin(4);
t205 = t190 * t102 + t201 * t103;
t204 = t218 * t103 - t216;
t203 = t219 * t102 + t217;
t200 = -t120 * t206 + t122 * t207;
t199 = t120 * t207 + t122 * t206;
t198 = Icges(3,3) + Icges(4,3);
t105 = sin(t115);
t106 = cos(t115);
t119 = sin(qJ(2));
t121 = cos(qJ(2));
t197 = Icges(3,5) * t121 + Icges(4,5) * t106 - Icges(3,6) * t119 - Icges(4,6) * t105;
t196 = t210 * t102 - t208 * t103;
t195 = t211 * t102 - t209 * t103;
t116 = t120 ^ 2;
t193 = t120 * pkin(6);
t192 = t215 * t102 + t214 * t103;
t191 = t120 * rSges(6,2);
t189 = t204 * t102 + t203 * t103;
t163 = t103 * t122;
t164 = t102 * t122;
t188 = t201 * t163 + t190 * t164 + t191;
t152 = rSges(4,1) * t106 - rSges(4,2) * t105;
t117 = t122 ^ 2;
t187 = t200 * t117 + (t196 * t120 + (-t195 + t199) * t122) * t120;
t186 = t198 * t120 + t197 * t122;
t185 = -t197 * t120 + t198 * t122;
t184 = (t199 * t116 + ((-t196 + t200) * t120 + t195 * t122) * t122) * t120;
t183 = t120 / 0.2e1;
t182 = -t122 / 0.2e1;
t181 = pkin(2) * t119;
t118 = -qJ(3) - pkin(6);
t104 = t121 * pkin(2) + pkin(1);
t126 = rSges(5,1) * t163 - rSges(5,2) * t164 + t120 * rSges(5,3);
t151 = rSges(5,1) * t103 - rSges(5,2) * t102;
t19 = t120 * (-t122 * rSges(5,3) + t151 * t120) + t122 * t126;
t113 = t122 * pkin(6);
t98 = t122 * t104;
t180 = t120 * (t113 + (-pkin(1) + t104) * t120) + t122 * (-t122 * pkin(1) - t193 + t98);
t179 = -t102 * t201 + t190 * t103;
t177 = rSges(3,1) * t121;
t175 = rSges(3,2) * t119;
t173 = t122 * rSges(3,3);
t172 = Icges(3,4) * t119;
t171 = Icges(3,4) * t121;
t170 = Icges(4,4) * t105;
t169 = Icges(4,4) * t106;
t162 = t120 * rSges(3,3) + t122 * t177;
t160 = t116 + t117;
t158 = -t105 * rSges(4,1) - t106 * rSges(4,2) - t181;
t84 = pkin(3) * t106 + t104;
t79 = t122 * t84;
t157 = t122 * (t79 - t98) + t180 + (-t104 + t84) * t116;
t8 = (t188 - t191) * t122 + t205 * t116;
t114 = -pkin(7) + t118;
t156 = -t120 * t114 + t79;
t154 = -pkin(3) * t105 - t181;
t153 = -t175 + t177;
t140 = t187 * t122 + t184;
t139 = Icges(3,1) * t121 - t172;
t138 = Icges(4,1) * t106 - t170;
t135 = -Icges(3,2) * t119 + t171;
t134 = -Icges(4,2) * t105 + t169;
t127 = t120 * rSges(4,3) + t152 * t122;
t78 = t102 * rSges(5,1) + t103 * rSges(5,2);
t125 = t154 - t78;
t124 = t154 + t179;
t123 = (t208 * t102 + t210 * t103 + t192 * t120 + t189 * t122) * t183 + (t209 * t102 + t211 * t103 + t189 * t120 - t192 * t122) * t182;
t96 = t122 * rSges(2,1) - t120 * rSges(2,2);
t95 = -t120 * rSges(2,1) - t122 * rSges(2,2);
t94 = t119 * rSges(3,1) + t121 * rSges(3,2);
t56 = t158 * t122;
t55 = t158 * t120;
t40 = t193 + (pkin(1) - t175) * t122 + t162;
t39 = t173 + t113 + (-pkin(1) - t153) * t120;
t32 = t125 * t122;
t31 = t125 * t120;
t30 = -t120 * t118 + t127 + t98;
t29 = (rSges(4,3) - t118) * t122 + (-t104 - t152) * t120;
t28 = t179 * t122;
t27 = t179 * t120;
t24 = t122 * (-t122 * t175 + t162) + (t153 * t120 - t173) * t120;
t23 = t126 + t156;
t22 = (rSges(5,3) - t114) * t122 + (-t151 - t84) * t120;
t21 = t124 * t122;
t20 = t124 * t120;
t18 = t156 + t188;
t17 = (rSges(6,2) - t114) * t122 + (-t205 - t84) * t120;
t7 = t122 * t127 + (-t122 * rSges(4,3) + t152 * t120) * t120 + t180;
t6 = t157 + t19;
t1 = t8 + t157;
t2 = [t105 * (Icges(4,1) * t105 + t169) + t106 * (Icges(4,2) * t106 + t170) + t119 * (Icges(3,1) * t119 + t171) + t121 * (Icges(3,2) * t121 + t172) + Icges(2,3) - t204 * t103 + t203 * t102 + m(5) * (t22 ^ 2 + t23 ^ 2) + m(6) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t39 ^ 2 + t40 ^ 2) + m(2) * (t95 ^ 2 + t96 ^ 2); m(3) * (-t120 * t40 - t122 * t39) * t94 + t123 + m(5) * (t32 * t22 + t31 * t23) + m(6) * (t21 * t17 + t20 * t18) + m(4) * (t56 * t29 + t55 * t30) + (t105 * (Icges(4,5) * t120 + t138 * t122) + t106 * (Icges(4,6) * t120 + t134 * t122) + t119 * (Icges(3,5) * t120 + t139 * t122) + t121 * (Icges(3,6) * t120 + t135 * t122)) * t183 + (t105 * (-Icges(4,5) * t122 + t138 * t120) + t106 * (-Icges(4,6) * t122 + t134 * t120) + t119 * (-Icges(3,5) * t122 + t139 * t120) + t121 * (-Icges(3,6) * t122 + t135 * t120)) * t182 + (Icges(3,5) * t119 + Icges(4,5) * t105 + Icges(3,6) * t121 + Icges(4,6) * t106) * (t117 / 0.2e1 + t116 / 0.2e1); m(6) * (t1 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2 + t6 ^ 2) + m(4) * (t55 ^ 2 + t56 ^ 2 + t7 ^ 2) + m(3) * (t160 * t94 ^ 2 + t24 ^ 2) + t184 + t186 * t120 * t116 + (t185 * t117 + (t185 * t120 + t186 * t122) * t120 + t187) * t122; m(5) * (t120 * t22 - t122 * t23) + m(6) * (t120 * t17 - t122 * t18) + m(4) * (t120 * t29 - t122 * t30); m(6) * (t120 * t21 - t122 * t20) + m(5) * (t120 * t32 - t122 * t31) + m(4) * (t120 * t56 - t122 * t55); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t160; m(6) * (t28 * t17 + t27 * t18) + m(5) * (-t120 * t23 - t122 * t22) * t78 + t123; m(6) * (t8 * t1 + t27 * t20 + t28 * t21) + m(5) * (t19 * t6 + (-t120 * t31 - t122 * t32) * t78) + t140; m(6) * (t28 * t120 - t27 * t122); m(5) * (t160 * t78 ^ 2 + t19 ^ 2) + m(6) * (t27 ^ 2 + t28 ^ 2 + t8 ^ 2) + t140; m(6) * (t120 * t18 + t122 * t17) * t102; m(6) * (-t103 * t1 + (t120 * t20 + t122 * t21) * t102); 0; m(6) * (-t103 * t8 + (t120 * t27 + t122 * t28) * t102); m(6) * (t160 * t102 ^ 2 + t103 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
