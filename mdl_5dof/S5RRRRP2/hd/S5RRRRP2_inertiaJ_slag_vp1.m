% Calculate joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:43
% EndTime: 2022-01-20 11:48:44
% DurationCPUTime: 1.70s
% Computational Cost: add. (2961->202), mult. (2286->282), div. (0->0), fcn. (2028->8), ass. (0->115)
t208 = Icges(5,4) + Icges(6,4);
t207 = Icges(5,1) + Icges(6,1);
t206 = Icges(5,2) + Icges(6,2);
t122 = qJ(3) + qJ(4);
t117 = cos(t122);
t205 = t208 * t117;
t115 = sin(t122);
t204 = t208 * t115;
t203 = Icges(5,5) + Icges(6,5);
t202 = Icges(5,6) + Icges(6,6);
t201 = -t206 * t115 + t205;
t200 = t207 * t117 - t204;
t123 = qJ(1) + qJ(2);
t116 = sin(t123);
t118 = cos(t123);
t199 = t201 * t116 - t202 * t118;
t198 = t202 * t116 + t201 * t118;
t197 = t200 * t116 - t203 * t118;
t196 = t203 * t116 + t200 * t118;
t195 = Icges(5,3) + Icges(6,3);
t194 = -t202 * t115 + t203 * t117;
t193 = t206 * t117 + t204;
t192 = t207 * t115 + t205;
t191 = t195 * t116 + t194 * t118;
t190 = -t194 * t116 + t195 * t118;
t189 = t198 * t115 - t196 * t117;
t188 = t199 * t115 - t197 * t117;
t187 = t203 * t115 + t202 * t117;
t128 = -pkin(8) - pkin(7);
t121 = qJ(5) - t128;
t161 = t117 * t118;
t163 = t115 * t118;
t126 = cos(qJ(3));
t112 = t126 * pkin(3) + pkin(2);
t86 = pkin(4) * t117 + t112;
t25 = rSges(6,1) * t161 - rSges(6,2) * t163 + t118 * t86 + (rSges(6,3) + t121) * t116;
t186 = -t193 * t115 + t192 * t117;
t114 = t118 ^ 2;
t185 = t190 * t114 + (t189 * t116 + (-t188 + t191) * t118) * t116;
t113 = t116 ^ 2;
t184 = (t191 * t113 + ((-t189 + t190) * t116 + t188 * t118) * t118) * t116;
t124 = sin(qJ(3));
t96 = t124 * rSges(4,1) + t126 * rSges(4,2);
t183 = m(4) * t96;
t81 = t115 * rSges(5,1) + t117 * rSges(5,2);
t182 = m(5) * t81;
t181 = t116 / 0.2e1;
t180 = -t118 / 0.2e1;
t179 = pkin(3) * t124;
t125 = sin(qJ(1));
t178 = t125 * pkin(1);
t110 = t118 * pkin(7);
t155 = -t118 * t112 + t116 * t128;
t159 = -t118 * pkin(2) - t116 * pkin(7);
t160 = t118 * t128;
t177 = t116 * (t160 + t110 + (-pkin(2) + t112) * t116) + t118 * (-t155 + t159);
t135 = rSges(5,1) * t161 - rSges(5,2) * t163 + t116 * rSges(5,3);
t162 = t116 * t117;
t164 = t115 * t116;
t172 = rSges(5,2) * t164 + t118 * rSges(5,3);
t16 = t116 * (rSges(5,1) * t162 - t172) + t118 * t135;
t176 = rSges(4,1) * t126;
t175 = rSges(4,2) * t124;
t173 = rSges(6,2) * t164 + t118 * rSges(6,3);
t171 = t118 * rSges(4,3) + t116 * t175;
t170 = Icges(4,4) * t124;
t169 = Icges(4,4) * t126;
t158 = t113 + t114;
t157 = -t81 - t179;
t156 = -t117 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t115;
t6 = (t155 + t25) * t118 + ((-t121 - t128) * t118 + rSges(6,1) * t162 - t173 + (-t112 + t86) * t116) * t116;
t83 = t118 * rSges(3,1) - t116 * rSges(3,2);
t82 = -t116 * rSges(3,1) - t118 * rSges(3,2);
t94 = Icges(4,2) * t126 + t170;
t95 = Icges(4,1) * t124 + t169;
t146 = -t124 * t94 + t126 * t95;
t145 = t185 * t118 + t184;
t144 = Icges(4,1) * t126 - t170;
t141 = -Icges(4,2) * t124 + t169;
t138 = Icges(4,5) * t126 - Icges(4,6) * t124;
t133 = t116 * rSges(4,3) + (-t175 + t176) * t118;
t132 = t156 - t179;
t131 = t192 * t115 + t193 * t117 + t124 * t95 + t126 * t94 + Icges(3,3);
t130 = (t196 * t115 + t187 * t116 + t198 * t117 + t186 * t118) * t181 + (t197 * t115 + t186 * t116 + t199 * t117 - t187 * t118) * t180;
t37 = t133 - t159;
t36 = t110 + (-pkin(2) - t176) * t116 + t171;
t31 = t135 - t155;
t24 = t118 * t121 + (-rSges(6,1) * t117 - t86) * t116 + t173;
t30 = -t160 + (-rSges(5,1) * t117 - t112) * t116 + t172;
t93 = Icges(4,5) * t124 + Icges(4,6) * t126;
t129 = t130 + (t116 * t93 + t146 * t118 + t124 * (Icges(4,5) * t116 + t144 * t118) + t126 * (Icges(4,6) * t116 + t141 * t118)) * t181 + (t146 * t116 - t118 * t93 + t124 * (-Icges(4,5) * t118 + t144 * t116) + t126 * (-Icges(4,6) * t118 + t141 * t116)) * t180;
t127 = cos(qJ(1));
t120 = t127 * pkin(1);
t98 = t127 * rSges(2,1) - t125 * rSges(2,2);
t97 = -t125 * rSges(2,1) - t127 * rSges(2,2);
t71 = t120 + t83;
t70 = t82 - t178;
t61 = Icges(4,3) * t116 + t138 * t118;
t60 = -Icges(4,3) * t118 + t138 * t116;
t59 = t157 * t118;
t58 = t157 * t116;
t45 = t156 * t118;
t44 = t156 * t116;
t35 = t132 * t118;
t34 = t132 * t116;
t33 = t120 + t37;
t32 = t36 - t178;
t27 = t120 + t31;
t26 = t30 - t178;
t21 = t120 + t25;
t20 = t24 - t178;
t19 = t116 * (t116 * t176 - t171) + t118 * t133;
t7 = t16 + t177;
t1 = t6 + t177;
t2 = [Icges(2,3) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(3) * (t70 ^ 2 + t71 ^ 2) + m(2) * (t97 ^ 2 + t98 ^ 2) + t131; m(6) * (t24 * t20 + t25 * t21) + m(5) * (t30 * t26 + t31 * t27) + m(4) * (t36 * t32 + t37 * t33) + m(3) * (t82 * t70 + t83 * t71) + t131; m(6) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t82 ^ 2 + t83 ^ 2) + t131; t129 + m(6) * (t35 * t20 + t34 * t21) + m(5) * (t59 * t26 + t58 * t27) + (-t116 * t33 - t118 * t32) * t183; t129 + m(6) * (t35 * t24 + t34 * t25) + m(5) * (t59 * t30 + t58 * t31) + (-t116 * t37 - t118 * t36) * t183; m(6) * (t1 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2 + t7 ^ 2) + m(4) * (t158 * t96 ^ 2 + t19 ^ 2) + t116 * t113 * t61 + t184 + (-t114 * t60 + (-t116 * t60 + t118 * t61) * t116 + t185) * t118; m(6) * (t45 * t20 + t44 * t21) + (-t116 * t27 - t118 * t26) * t182 + t130; m(6) * (t45 * t24 + t44 * t25) + (-t116 * t31 - t118 * t30) * t182 + t130; m(6) * (t6 * t1 + t44 * t34 + t45 * t35) + m(5) * (t16 * t7 + (-t116 * t58 - t118 * t59) * t81) + t145; m(5) * (t158 * t81 ^ 2 + t16 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2 + t6 ^ 2) + t145; m(6) * (t116 * t20 - t118 * t21); m(6) * (t116 * t24 - t118 * t25); m(6) * (t116 * t35 - t118 * t34); m(6) * (t116 * t45 - t118 * t44); m(6) * t158;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
