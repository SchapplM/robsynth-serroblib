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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:10
% DurationCPUTime: 2.03s
% Computational Cost: add. (2961->204), mult. (2286->283), div. (0->0), fcn. (2028->8), ass. (0->120)
t211 = Icges(5,4) + Icges(6,4);
t210 = Icges(5,1) + Icges(6,1);
t209 = Icges(5,2) + Icges(6,2);
t128 = qJ(3) + qJ(4);
t122 = cos(t128);
t208 = t211 * t122;
t120 = sin(t128);
t207 = t211 * t120;
t206 = Icges(5,5) + Icges(6,5);
t205 = -Icges(5,6) - Icges(6,6);
t204 = -t209 * t120 + t208;
t203 = t210 * t122 - t207;
t129 = qJ(1) + qJ(2);
t121 = sin(t129);
t123 = cos(t129);
t202 = t204 * t121 + t205 * t123;
t201 = t205 * t121 - t204 * t123;
t200 = t203 * t121 - t206 * t123;
t199 = -t206 * t121 - t203 * t123;
t198 = Icges(5,3) + Icges(6,3);
t197 = t205 * t120 + t206 * t122;
t196 = t209 * t122 + t207;
t195 = t210 * t120 + t208;
t194 = -t197 * t121 + t198 * t123;
t193 = t198 * t121 + t197 * t123;
t192 = -t201 * t120 + t199 * t122;
t191 = t202 * t120 - t200 * t122;
t190 = -t206 * t120 + t205 * t122;
t134 = -pkin(8) - pkin(7);
t127 = -qJ(5) + t134;
t167 = t121 * t122;
t169 = t120 * t121;
t132 = cos(qJ(3));
t117 = t132 * pkin(3) + pkin(2);
t88 = pkin(4) * t122 + t117;
t24 = rSges(6,1) * t167 - rSges(6,2) * t169 + t121 * t88 + (-rSges(6,3) + t127) * t123;
t189 = t196 * t120 - t195 * t122;
t119 = t123 ^ 2;
t188 = t194 * t119 + (t192 * t121 + (-t191 + t193) * t123) * t121;
t118 = t121 ^ 2;
t187 = t193 * t118 + (t191 * t123 + (-t192 + t194) * t121) * t123;
t166 = t122 * t123;
t168 = t120 * t123;
t186 = rSges(6,1) * t166 - rSges(6,2) * t168 + t121 * rSges(6,3) + t123 * t88;
t83 = t120 * rSges(5,1) + t122 * rSges(5,2);
t185 = m(5) * t83;
t184 = -t121 / 0.2e1;
t183 = -t123 / 0.2e1;
t130 = sin(qJ(3));
t100 = t130 * rSges(4,1) + t132 * rSges(4,2);
t182 = m(4) * t100;
t181 = pkin(3) * t130;
t176 = t121 * t117 + t123 * t134;
t180 = (-t176 + t24) * t121;
t90 = t123 * t117;
t179 = -t90 - (t127 - t134) * t121 + t186;
t82 = t120 * rSges(6,1) + t122 * rSges(6,2);
t43 = pkin(4) * t168 + t123 * t82;
t178 = rSges(4,1) * t132;
t175 = Icges(4,4) * t130;
t174 = Icges(4,4) * t132;
t165 = t123 * t130;
t84 = t121 * rSges(3,1) + t123 * rSges(3,2);
t164 = t123 * pkin(2) + t121 * pkin(7);
t163 = t119 + t118;
t162 = -pkin(4) * t120 - t82;
t85 = t123 * rSges(3,1) - t121 * rSges(3,2);
t161 = t121 * t134 - t90;
t160 = (-rSges(4,2) * t130 + t178) * t121;
t98 = Icges(4,2) * t132 + t175;
t99 = Icges(4,1) * t130 + t174;
t151 = t130 * t98 - t132 * t99;
t150 = Icges(4,1) * t132 - t175;
t147 = -Icges(4,2) * t130 + t174;
t144 = Icges(4,5) * t132 - Icges(4,6) * t130;
t57 = -rSges(5,1) * t166 + rSges(5,2) * t168 - t121 * rSges(5,3);
t141 = rSges(4,2) * t165 - t121 * rSges(4,3) - t123 * t178;
t140 = t195 * t120 + t196 * t122 + t130 * t99 + t132 * t98 + Icges(3,3);
t139 = (t199 * t120 + t190 * t121 + t201 * t122 + t189 * t123) * t184 + (t200 * t120 - t189 * t121 + t202 * t122 + t190 * t123) * t183;
t138 = rSges(5,1) * t167 - rSges(5,2) * t169 - t123 * rSges(5,3);
t136 = t187 * t121 + t188 * t123;
t37 = -t141 + t164;
t29 = t138 + t176;
t30 = -t161 - t57;
t25 = -t121 * t127 + t186;
t115 = t121 * pkin(2);
t36 = t115 + (-rSges(4,3) - pkin(7)) * t123 + t160;
t97 = Icges(4,5) * t130 + Icges(4,6) * t132;
t135 = t139 + (-t121 * t97 + t151 * t123 + t130 * (-Icges(4,5) * t121 - t150 * t123) + t132 * (-Icges(4,6) * t121 - t147 * t123)) * t184 + (-t151 * t121 - t123 * t97 + t130 * (-Icges(4,5) * t123 + t150 * t121) + t132 * (-Icges(4,6) * t123 + t147 * t121)) * t183;
t133 = cos(qJ(1));
t131 = sin(qJ(1));
t126 = t133 * pkin(1);
t124 = t131 * pkin(1);
t106 = pkin(3) * t165;
t102 = t133 * rSges(2,1) - t131 * rSges(2,2);
t101 = t131 * rSges(2,1) + t133 * rSges(2,2);
t72 = t126 + t85;
t71 = t124 + t84;
t61 = -Icges(4,3) * t121 - t144 * t123;
t60 = -Icges(4,3) * t123 + t144 * t121;
t59 = t123 * t83 + t106;
t58 = (-t83 - t181) * t121;
t42 = t162 * t121;
t41 = t161 + t164;
t40 = t121 * t138;
t38 = t121 * (t123 * pkin(7) - t115 + t176);
t35 = t106 + t43;
t34 = (t162 - t181) * t121;
t33 = t126 + t37;
t32 = t124 + t36;
t27 = t126 + t30;
t26 = t124 + t29;
t21 = t126 + t25;
t20 = t124 + t24;
t19 = -t123 * t141 + t121 * (-t123 * rSges(4,3) + t160);
t16 = -t123 * t57 + t40;
t7 = t38 + t40 + (-t41 - t57) * t123;
t6 = t179 * t123 + t180;
t1 = t38 + (-t41 + t179) * t123 + t180;
t2 = [Icges(2,3) + m(2) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t71 ^ 2 + t72 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2) + t140; m(3) * (t84 * t71 + t85 * t72) + m(4) * (t36 * t32 + t37 * t33) + m(5) * (t29 * t26 + t30 * t27) + m(6) * (t24 * t20 + t25 * t21) + t140; m(6) * (t24 ^ 2 + t25 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t84 ^ 2 + t85 ^ 2) + t140; t135 + m(5) * (t59 * t26 + t58 * t27) + m(6) * (t35 * t20 + t34 * t21) + (-t121 * t33 + t123 * t32) * t182; t135 + m(6) * (t35 * t24 + t34 * t25) + m(5) * (t59 * t29 + t58 * t30) + (-t121 * t37 + t123 * t36) * t182; m(6) * (t1 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2 + t7 ^ 2) + m(4) * (t163 * t100 ^ 2 + t19 ^ 2) + (-t118 * t61 + t187) * t121 + (-t119 * t60 + (-t121 * t60 - t123 * t61) * t121 + t188) * t123; m(6) * (t43 * t20 + t42 * t21) + (-t121 * t27 + t123 * t26) * t185 + t139; m(6) * (t43 * t24 + t42 * t25) + (-t121 * t30 + t123 * t29) * t185 + t139; m(6) * (t6 * t1 + t42 * t34 + t43 * t35) + m(5) * (t16 * t7 + (-t121 * t58 + t123 * t59) * t83) + t136; m(5) * (t163 * t83 ^ 2 + t16 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2 + t6 ^ 2) + t136; m(6) * (-t121 * t20 - t123 * t21); m(6) * (-t121 * t24 - t123 * t25); m(6) * (-t121 * t35 - t123 * t34); m(6) * (-t121 * t43 - t123 * t42); m(6) * t163;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
