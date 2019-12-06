% Calculate joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:43
% DurationCPUTime: 1.93s
% Computational Cost: add. (4510->226), mult. (6185->352), div. (0->0), fcn. (6784->8), ass. (0->121)
t124 = qJ(2) + qJ(3);
t121 = cos(t124);
t126 = cos(pkin(8));
t129 = cos(qJ(4));
t157 = t126 * t129;
t125 = sin(pkin(8));
t127 = sin(qJ(4));
t160 = t125 * t127;
t109 = t121 * t160 + t157;
t158 = t126 * t127;
t159 = t125 * t129;
t110 = t121 * t159 - t158;
t120 = sin(t124);
t163 = t120 * t125;
t57 = Icges(6,5) * t110 + Icges(6,6) * t163 + Icges(6,3) * t109;
t63 = Icges(5,4) * t110 - Icges(5,2) * t109 + Icges(5,6) * t163;
t200 = t57 - t63;
t111 = t121 * t158 - t159;
t112 = t121 * t157 + t160;
t162 = t120 * t126;
t58 = Icges(6,5) * t112 + Icges(6,6) * t162 + Icges(6,3) * t111;
t64 = Icges(5,4) * t112 - Icges(5,2) * t111 + Icges(5,6) * t162;
t199 = t58 - t64;
t59 = Icges(5,5) * t110 - Icges(5,6) * t109 + Icges(5,3) * t163;
t61 = Icges(6,4) * t110 + Icges(6,2) * t163 + Icges(6,6) * t109;
t198 = t59 + t61;
t60 = Icges(5,5) * t112 - Icges(5,6) * t111 + Icges(5,3) * t162;
t62 = Icges(6,4) * t112 + Icges(6,2) * t162 + Icges(6,6) * t111;
t197 = t60 + t62;
t65 = Icges(6,1) * t110 + Icges(6,4) * t163 + Icges(6,5) * t109;
t67 = Icges(5,1) * t110 - Icges(5,4) * t109 + Icges(5,5) * t163;
t196 = t65 + t67;
t66 = Icges(6,1) * t112 + Icges(6,4) * t162 + Icges(6,5) * t111;
t68 = Icges(5,1) * t112 - Icges(5,4) * t111 + Icges(5,5) * t162;
t195 = t66 + t68;
t194 = rSges(6,1) + pkin(4);
t193 = -t200 * t109 - t196 * t110 - t198 * t163;
t192 = t199 * t109 + t195 * t110 + t197 * t163;
t191 = t200 * t111 + t196 * t112 + t198 * t162;
t190 = t199 * t111 + t195 * t112 + t197 * t162;
t84 = -Icges(6,6) * t121 + (Icges(6,5) * t129 + Icges(6,3) * t127) * t120;
t87 = -Icges(5,6) * t121 + (Icges(5,4) * t129 - Icges(5,2) * t127) * t120;
t189 = -t84 + t87;
t88 = -Icges(6,4) * t121 + (Icges(6,1) * t129 + Icges(6,5) * t127) * t120;
t89 = -Icges(5,5) * t121 + (Icges(5,1) * t129 - Icges(5,4) * t127) * t120;
t188 = -t88 - t89;
t187 = rSges(6,3) + qJ(5);
t85 = -Icges(5,3) * t121 + (Icges(5,5) * t129 - Icges(5,6) * t127) * t120;
t86 = -Icges(6,2) * t121 + (Icges(6,4) * t129 + Icges(6,6) * t127) * t120;
t186 = (-t86 - t85) * t121;
t122 = t125 ^ 2;
t123 = t126 ^ 2;
t185 = t122 + t123;
t184 = (t189 * t109 + t188 * t110) * t121 + (t192 * t126 + (t186 - t193) * t125) * t120;
t183 = (t189 * t111 + t188 * t112) * t121 + ((t186 + t190) * t126 + t191 * t125) * t120;
t182 = t192 * t125 + t193 * t126;
t181 = t190 * t125 - t191 * t126;
t170 = rSges(6,2) * t163 + t187 * t109 + t194 * t110;
t180 = rSges(6,2) * t162 + t187 * t111 + t194 * t112;
t139 = Icges(4,4) * t121 - Icges(4,2) * t120;
t141 = Icges(4,1) * t121 - Icges(4,4) * t120;
t143 = -t120 * (Icges(4,6) * t125 + t139 * t126) + t121 * (Icges(4,5) * t125 + t141 * t126);
t144 = t120 * (-Icges(4,6) * t126 + t139 * t125) - t121 * (-Icges(4,5) * t126 + t141 * t125);
t137 = Icges(4,5) * t121 - Icges(4,6) * t120;
t94 = -Icges(4,3) * t126 + t137 * t125;
t95 = Icges(4,3) * t125 + t137 * t126;
t178 = -t123 * t94 - (t143 * t125 + (t144 - t95) * t126) * t125 - t182;
t177 = t121 ^ 2;
t128 = sin(qJ(2));
t173 = pkin(2) * t128;
t130 = cos(qJ(2));
t168 = t185 * pkin(2) * t130;
t50 = t185 * (rSges(4,1) * t121 - rSges(4,2) * t120);
t165 = -rSges(6,2) * t121 + (t187 * t127 + t194 * t129) * t120;
t115 = pkin(3) * t120 - pkin(7) * t121;
t91 = -rSges(5,3) * t121 + (rSges(5,1) * t129 - rSges(5,2) * t127) * t120;
t164 = -t115 - t91;
t161 = t120 * t127;
t156 = t185 * (pkin(3) * t121 + pkin(7) * t120);
t154 = (t122 * t95 + (t144 * t126 + (t143 - t94) * t125) * t126 + t181) * t125;
t153 = -t115 - t165;
t114 = rSges(4,1) * t120 + rSges(4,2) * t121;
t150 = -t114 - t173;
t149 = -t115 - t173;
t70 = rSges(5,1) * t110 - rSges(5,2) * t109 + rSges(5,3) * t163;
t72 = rSges(5,1) * t112 - rSges(5,2) * t111 + rSges(5,3) * t162;
t41 = t125 * t70 + t126 * t72 + t156;
t148 = t149 - t91;
t138 = Icges(3,5) * t130 - Icges(3,6) * t128;
t134 = t149 - t165;
t22 = t170 * t125 + t180 * t126 + t156;
t133 = t178 * t126 + t154;
t35 = -t121 * t61 + (t127 * t57 + t129 * t65) * t120;
t36 = -t121 * t62 + (t127 * t58 + t129 * t66) * t120;
t37 = -t121 * t59 + (-t127 * t63 + t129 * t67) * t120;
t38 = -t121 * t60 + (-t127 * t64 + t129 * t68) * t120;
t132 = -((-t35 - t37) * t126 + (t36 + t38) * t125) * t121 / 0.2e1 + t183 * t125 / 0.2e1 - t184 * t126 / 0.2e1 + t182 * t163 / 0.2e1 + t181 * t162 / 0.2e1;
t117 = rSges(3,1) * t128 + rSges(3,2) * t130;
t103 = Icges(3,3) * t125 + t138 * t126;
t102 = -Icges(3,3) * t126 + t138 * t125;
t93 = t150 * t126;
t92 = t150 * t125;
t75 = t164 * t126;
t74 = t164 * t125;
t73 = t185 * (rSges(3,1) * t130 - rSges(3,2) * t128);
t52 = t148 * t126;
t51 = t148 * t125;
t49 = t153 * t126;
t48 = t153 * t125;
t47 = t134 * t126;
t46 = t134 * t125;
t45 = -t121 * t72 - t91 * t162;
t44 = t121 * t70 + t91 * t163;
t43 = t50 + t168;
t42 = (-t125 * t72 + t126 * t70) * t120;
t40 = -t121 * t180 - t165 * t162;
t39 = t170 * t121 + t165 * t163;
t24 = t41 + t168;
t23 = (-t125 * t180 + t170 * t126) * t120;
t21 = t22 + t168;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t73 + m(4) * t43 + m(5) * t24 + m(6) * t21; m(6) * (t21 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t24 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(4) * (t43 ^ 2 + t92 ^ 2 + t93 ^ 2) + t125 * t122 * t103 + m(3) * (t117 ^ 2 * t185 + t73 ^ 2) + t154 + (-t123 * t102 + (-t125 * t102 + t126 * t103) * t125 + t178) * t126; m(4) * t50 + m(5) * t41 + m(6) * t22; m(6) * (t21 * t22 + t46 * t48 + t47 * t49) + m(5) * (t24 * t41 + t51 * t74 + t52 * t75) + m(4) * (t43 * t50 + (-t125 * t92 - t126 * t93) * t114) + t133; m(6) * (t22 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t41 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t114 ^ 2 * t185 + t50 ^ 2) + t133; m(5) * t42 + m(6) * t23; m(6) * (t21 * t23 + t39 * t47 + t40 * t46) + m(5) * (t24 * t42 + t44 * t52 + t45 * t51) + t132; m(6) * (t22 * t23 + t39 * t49 + t40 * t48) + m(5) * (t41 * t42 + t44 * t75 + t45 * t74) + t132; m(6) * (t23 ^ 2 + t39 ^ 2 + t40 ^ 2) - t121 * (t177 * t86 + (t36 * t126 + t35 * t125 - (t127 * t84 + t129 * t88) * t121) * t120) - t121 * (t177 * t85 + (t38 * t126 + t37 * t125 - (-t127 * t87 + t129 * t89) * t121) * t120) + m(5) * (t42 ^ 2 + t44 ^ 2 + t45 ^ 2) + t184 * t163 + t183 * t162; m(6) * t161; m(6) * (t109 * t46 + t111 * t47 + t21 * t161); m(6) * (t109 * t48 + t111 * t49 + t22 * t161); m(6) * (t109 * t40 + t111 * t39 + t23 * t161); m(6) * (t120 ^ 2 * t127 ^ 2 + t109 ^ 2 + t111 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
