% Calculate joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:49
% EndTime: 2022-01-20 11:16:50
% DurationCPUTime: 1.06s
% Computational Cost: add. (4216->255), mult. (4188->368), div. (0->0), fcn. (4406->10), ass. (0->119)
t124 = cos(qJ(4));
t112 = t124 * pkin(4) + pkin(3);
t119 = qJ(1) + qJ(2);
t116 = cos(t119);
t120 = sin(pkin(9));
t126 = -pkin(8) - pkin(7);
t141 = t120 * t126;
t121 = cos(pkin(9));
t143 = t116 * t121;
t114 = sin(t119);
t122 = sin(qJ(4));
t145 = t114 * t122;
t144 = t116 * t120;
t118 = qJ(4) + qJ(5);
t113 = sin(t118);
t115 = cos(t118);
t75 = -t113 * t143 + t114 * t115;
t76 = t114 * t113 + t115 * t143;
t48 = t76 * rSges(6,1) + t75 * rSges(6,2) + rSges(6,3) * t144;
t155 = pkin(4) * t145 + t112 * t143 - t116 * t141 + t48;
t123 = sin(qJ(1));
t154 = t123 * pkin(1);
t153 = -pkin(3) + t112;
t138 = pkin(3) * t143 + pkin(7) * t144;
t152 = t138 - t155;
t147 = t114 * t120;
t146 = t114 * t121;
t73 = -t113 * t146 - t116 * t115;
t74 = -t116 * t113 + t115 * t146;
t131 = -t74 * rSges(6,1) - t73 * rSges(6,2);
t47 = rSges(6,3) * t147 - t131;
t72 = -t121 * rSges(6,3) + (rSges(6,1) * t115 - rSges(6,2) * t113) * t120;
t29 = t121 * t47 + t72 * t147;
t70 = -Icges(6,6) * t121 + (Icges(6,4) * t115 - Icges(6,2) * t113) * t120;
t151 = t113 * t70;
t80 = -Icges(5,6) * t121 + (Icges(5,4) * t124 - Icges(5,2) * t122) * t120;
t150 = t122 * t80;
t71 = -Icges(6,5) * t121 + (Icges(6,1) * t115 - Icges(6,4) * t113) * t120;
t63 = t120 * t115 * t71;
t69 = -Icges(6,3) * t121 + (Icges(6,5) * t115 - Icges(6,6) * t113) * t120;
t31 = -t120 * t151 - t121 * t69 + t63;
t149 = t31 * t121;
t142 = t116 * t122;
t148 = -pkin(4) * t142 - t114 * t141;
t140 = t121 * t122;
t139 = t121 * t124;
t137 = t116 * pkin(2) + t114 * qJ(3);
t85 = t114 * t124 - t116 * t140;
t86 = t116 * t139 + t145;
t58 = t86 * rSges(5,1) + t85 * rSges(5,2) + rSges(5,3) * t144;
t136 = t147 / 0.2e1;
t135 = t144 / 0.2e1;
t18 = t69 * t147 + t73 * t70 + t74 * t71;
t19 = t69 * t144 + t75 * t70 + t76 * t71;
t41 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t147;
t43 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t147;
t45 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t147;
t7 = -t121 * t41 + (-t113 * t43 + t115 * t45) * t120;
t42 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t144;
t44 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t144;
t46 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t144;
t8 = -t121 * t42 + (-t113 * t44 + t115 * t46) * t120;
t134 = (t18 + t7) * t136 + (t19 + t8) * t135;
t93 = t116 * rSges(3,1) - t114 * rSges(3,2);
t133 = -t121 * (-t149 + (t114 * t7 + t116 * t8) * t120) + ((t42 * t147 + t73 * t44 + t74 * t46) * t144 + (t41 * t147 + t73 * t43 + t74 * t45) * t147 - t18 * t121) * t147 + ((t42 * t144 + t75 * t44 + t76 * t46) * t144 + (t41 * t144 + t75 * t43 + t76 * t45) * t147 - t19 * t121) * t144;
t83 = -t114 * t140 - t116 * t124;
t84 = t114 * t139 - t142;
t132 = -t84 * rSges(5,1) - t83 * rSges(5,2);
t92 = -t114 * rSges(3,1) - t116 * rSges(3,2);
t129 = t134 - t149;
t37 = t58 + t137 + t138;
t62 = rSges(4,1) * t143 - rSges(4,2) * t144 + t114 * rSges(4,3) + t137;
t107 = t116 * qJ(3);
t61 = t107 + t116 * rSges(4,3) + rSges(4,2) * t147 + (-rSges(4,1) * t121 - pkin(2)) * t114;
t28 = t137 + t155;
t49 = Icges(5,5) * t84 + Icges(5,6) * t83 + Icges(5,3) * t147;
t51 = Icges(5,4) * t84 + Icges(5,2) * t83 + Icges(5,6) * t147;
t53 = Icges(5,1) * t84 + Icges(5,4) * t83 + Icges(5,5) * t147;
t13 = -t121 * t49 + (-t122 * t51 + t124 * t53) * t120;
t50 = Icges(5,5) * t86 + Icges(5,6) * t85 + Icges(5,3) * t144;
t52 = Icges(5,4) * t86 + Icges(5,2) * t85 + Icges(5,6) * t144;
t54 = Icges(5,1) * t86 + Icges(5,4) * t85 + Icges(5,5) * t144;
t14 = -t121 * t50 + (-t122 * t52 + t124 * t54) * t120;
t79 = -Icges(5,3) * t121 + (Icges(5,5) * t124 - Icges(5,6) * t122) * t120;
t81 = -Icges(5,5) * t121 + (Icges(5,1) * t124 - Icges(5,4) * t122) * t120;
t23 = t79 * t147 + t83 * t80 + t84 * t81;
t24 = t79 * t144 + t85 * t80 + t86 * t81;
t65 = t120 * t124 * t81;
t38 = -t120 * t150 - t121 * t79 + t65;
t128 = (-t31 - t38) * t121 + t134 + (t13 + t23) * t136 + (t14 + t24) * t135;
t36 = t107 + (-pkin(3) * t121 - pkin(2) + (-rSges(5,3) - pkin(7)) * t120) * t114 + t132;
t27 = t107 + (-rSges(6,3) * t120 - t112 * t121 - pkin(2)) * t114 + t131 - t148;
t127 = Icges(3,3) + t63 + t65 + (Icges(4,2) * t121 - t69 - t79) * t121 + (Icges(4,1) * t120 + 0.2e1 * Icges(4,4) * t121 - t150 - t151) * t120;
t125 = cos(qJ(1));
t117 = t125 * pkin(1);
t97 = t125 * rSges(2,1) - t123 * rSges(2,2);
t96 = -t123 * rSges(2,1) - t125 * rSges(2,2);
t88 = t117 + t93;
t87 = t92 - t154;
t82 = -t121 * rSges(5,3) + (rSges(5,1) * t124 - rSges(5,2) * t122) * t120;
t66 = (pkin(7) + t126) * t121 + t153 * t120;
t60 = t117 + t62;
t59 = t61 - t154;
t57 = rSges(5,3) * t147 - t132;
t55 = (-pkin(7) * t120 + t153 * t121) * t114 + t148;
t39 = t47 * t144;
t35 = -t121 * t58 - t82 * t144;
t34 = t121 * t57 + t82 * t147;
t33 = t117 + t37;
t32 = t36 - t154;
t30 = -t121 * t48 - t72 * t144;
t26 = t117 + t28;
t25 = t27 - t154;
t20 = (-t114 * t58 + t116 * t57) * t120;
t17 = -t48 * t147 + t39;
t10 = t152 * t121 + (-t66 - t72) * t144;
t9 = t121 * t55 + t66 * t147 + t29;
t4 = t39 + (t152 * t114 + t116 * t55) * t120;
t1 = [Icges(2,3) + m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(4) * (t59 ^ 2 + t60 ^ 2) + m(3) * (t87 ^ 2 + t88 ^ 2) + m(2) * (t96 ^ 2 + t97 ^ 2) + t127; m(6) * (t27 * t25 + t28 * t26) + m(5) * (t36 * t32 + t37 * t33) + m(4) * (t61 * t59 + t62 * t60) + m(3) * (t92 * t87 + t93 * t88) + t127; m(6) * (t27 ^ 2 + t28 ^ 2) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t61 ^ 2 + t62 ^ 2) + m(3) * (t92 ^ 2 + t93 ^ 2) + t127; m(6) * (t114 * t25 - t116 * t26) + m(5) * (t114 * t32 - t116 * t33) + m(4) * (t114 * t59 - t116 * t60); m(6) * (t114 * t27 - t116 * t28) + m(5) * (t114 * t36 - t116 * t37) + m(4) * (t114 * t61 - t116 * t62); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t114 ^ 2 + t116 ^ 2); m(6) * (t10 * t26 + t9 * t25) + m(5) * (t34 * t32 + t35 * t33) + t128; m(6) * (t10 * t28 + t9 * t27) + m(5) * (t34 * t36 + t35 * t37) + t128; m(5) * (t34 * t114 - t35 * t116) + m(6) * (-t10 * t116 + t9 * t114); ((t50 * t144 + t85 * t52 + t86 * t54) * t144 + (t49 * t144 + t85 * t51 + t86 * t53) * t147 - t24 * t121) * t144 + ((t50 * t147 + t83 * t52 + t84 * t54) * t144 + (t49 * t147 + t83 * t51 + t84 * t53) * t147 - t23 * t121) * t147 - t121 * (-t38 * t121 + (t114 * t13 + t116 * t14) * t120) + m(5) * (t20 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t10 ^ 2 + t4 ^ 2 + t9 ^ 2) + t133; m(6) * (t29 * t25 + t30 * t26) + t129; m(6) * (t29 * t27 + t30 * t28) + t129; m(6) * (t29 * t114 - t30 * t116); m(6) * (t30 * t10 + t17 * t4 + t29 * t9) + t133; m(6) * (t17 ^ 2 + t29 ^ 2 + t30 ^ 2) + t133;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
