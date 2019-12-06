% Calculate joint inertia matrix for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:27
% EndTime: 2019-12-05 15:53:33
% DurationCPUTime: 1.53s
% Computational Cost: add. (5194->266), mult. (5729->433), div. (0->0), fcn. (6207->10), ass. (0->133)
t127 = sin(pkin(8));
t122 = t127 ^ 2;
t129 = cos(pkin(8));
t123 = t129 ^ 2;
t149 = t122 + t123;
t132 = cos(qJ(2));
t171 = t132 ^ 2;
t170 = t127 / 0.2e1;
t169 = -t129 / 0.2e1;
t168 = -t132 / 0.2e1;
t128 = cos(pkin(9));
t166 = t128 * pkin(3);
t131 = sin(qJ(2));
t124 = pkin(9) + qJ(4);
t118 = cos(t124);
t151 = pkin(4) * t118;
t133 = pkin(7) * t131 + t132 * t151;
t117 = sin(t124);
t146 = pkin(4) * t117;
t156 = t129 * t131;
t119 = qJ(5) + t124;
t114 = sin(t119);
t115 = cos(t119);
t155 = t129 * t132;
t89 = -t114 * t155 + t115 * t127;
t90 = t114 * t127 + t115 * t155;
t56 = rSges(6,1) * t90 + rSges(6,2) * t89 + rSges(6,3) * t156;
t165 = t127 * t146 + t129 * t133 + t56;
t74 = -pkin(7) * t132 + t131 * t151;
t81 = -rSges(6,3) * t132 + (rSges(6,1) * t115 - rSges(6,2) * t114) * t131;
t164 = -t74 - t81;
t159 = t127 * t131;
t158 = t127 * t132;
t87 = -t114 * t158 - t115 * t129;
t88 = -t114 * t129 + t115 * t158;
t55 = rSges(6,1) * t88 + rSges(6,2) * t87 + rSges(6,3) * t159;
t39 = t132 * t55 + t159 * t81;
t78 = -Icges(6,3) * t132 + (Icges(6,5) * t115 - Icges(6,6) * t114) * t131;
t163 = t132 * t78;
t83 = -Icges(5,3) * t132 + (Icges(5,5) * t118 - Icges(5,6) * t117) * t131;
t162 = t132 * t83;
t112 = t131 * pkin(2) - qJ(3) * t132;
t161 = pkin(6) * t132 - t131 * t166 - t112;
t126 = sin(pkin(9));
t160 = t127 * t126;
t157 = t129 * t126;
t153 = rSges(4,3) * t132 - (rSges(4,1) * t128 - rSges(4,2) * t126) * t131 - t112;
t152 = t149 * (pkin(2) * t132 + qJ(3) * t131);
t148 = -m(4) - m(5) - m(6);
t86 = -rSges(5,3) * t132 + (rSges(5,1) * t118 - rSges(5,2) * t117) * t131;
t147 = -t86 + t161;
t145 = t161 + t164;
t134 = pkin(6) * t131 + t132 * t166;
t144 = t127 * (-pkin(3) * t157 + t127 * t134) + t129 * (pkin(3) * t160 + t129 * t134) + t152;
t49 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t159;
t51 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t159;
t53 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t159;
t30 = -t132 * t49 + (-t114 * t51 + t115 * t53) * t131;
t50 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t156;
t52 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t156;
t54 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t156;
t31 = -t132 * t50 + (-t114 * t52 + t115 * t54) * t131;
t20 = t159 * t49 + t51 * t87 + t53 * t88;
t21 = t159 * t50 + t52 * t87 + t54 * t88;
t79 = -Icges(6,6) * t132 + (Icges(6,4) * t115 - Icges(6,2) * t114) * t131;
t80 = -Icges(6,5) * t132 + (Icges(6,1) * t115 - Icges(6,4) * t114) * t131;
t5 = -(t79 * t87 + t80 * t88) * t132 + (t21 * t129 + (t20 - t163) * t127) * t131;
t22 = t156 * t49 + t51 * t89 + t53 * t90;
t23 = t156 * t50 + t52 * t89 + t54 * t90;
t6 = -(t79 * t89 + t80 * t90) * t132 + (t22 * t127 + (t23 - t163) * t129) * t131;
t143 = -t132 * (t171 * t78 + (t31 * t129 + t30 * t127 - (-t114 * t79 + t115 * t80) * t132) * t131) + t6 * t156 + t5 * t159;
t12 = t127 * t21 - t129 * t20;
t13 = t127 * t23 - t129 * t22;
t142 = t12 * t159 / 0.2e1 + t5 * t169 + t6 * t170 + t13 * t156 / 0.2e1 + (t127 * t31 - t129 * t30) * t168;
t135 = Icges(3,5) * t132 - Icges(3,6) * t131;
t113 = rSges(3,1) * t131 + rSges(3,2) * t132;
t108 = t128 * t155 + t160;
t107 = -t126 * t155 + t127 * t128;
t106 = t128 * t158 - t157;
t105 = -t126 * t158 - t128 * t129;
t100 = t117 * t127 + t118 * t155;
t99 = -t117 * t155 + t118 * t127;
t98 = -t117 * t129 + t118 * t158;
t97 = -t117 * t158 - t118 * t129;
t92 = Icges(3,3) * t127 + t129 * t135;
t91 = -Icges(3,3) * t129 + t127 * t135;
t85 = -Icges(5,5) * t132 + (Icges(5,1) * t118 - Icges(5,4) * t117) * t131;
t84 = -Icges(5,6) * t132 + (Icges(5,4) * t118 - Icges(5,2) * t117) * t131;
t76 = t153 * t129;
t75 = t153 * t127;
t73 = Icges(4,1) * t108 + Icges(4,4) * t107 + Icges(4,5) * t156;
t72 = Icges(4,1) * t106 + Icges(4,4) * t105 + Icges(4,5) * t159;
t71 = Icges(4,4) * t108 + Icges(4,2) * t107 + Icges(4,6) * t156;
t70 = Icges(4,4) * t106 + Icges(4,2) * t105 + Icges(4,6) * t159;
t69 = Icges(4,5) * t108 + Icges(4,6) * t107 + Icges(4,3) * t156;
t68 = Icges(4,5) * t106 + Icges(4,6) * t105 + Icges(4,3) * t159;
t65 = t149 * (rSges(3,1) * t132 - rSges(3,2) * t131);
t64 = rSges(5,1) * t100 + rSges(5,2) * t99 + rSges(5,3) * t156;
t63 = rSges(5,1) * t98 + rSges(5,2) * t97 + rSges(5,3) * t159;
t62 = Icges(5,1) * t100 + Icges(5,4) * t99 + Icges(5,5) * t156;
t61 = Icges(5,1) * t98 + Icges(5,4) * t97 + Icges(5,5) * t159;
t60 = Icges(5,4) * t100 + Icges(5,2) * t99 + Icges(5,6) * t156;
t59 = Icges(5,4) * t98 + Icges(5,2) * t97 + Icges(5,6) * t159;
t58 = Icges(5,5) * t100 + Icges(5,6) * t99 + Icges(5,3) * t156;
t57 = Icges(5,5) * t98 + Icges(5,6) * t97 + Icges(5,3) * t159;
t47 = t55 * t156;
t46 = t147 * t129;
t45 = t147 * t127;
t43 = t127 * t133 - t129 * t146;
t42 = -t132 * t64 - t156 * t86;
t41 = t132 * t63 + t159 * t86;
t40 = -t132 * t56 - t156 * t81;
t38 = t145 * t129;
t37 = t145 * t127;
t36 = (-t127 * t64 + t129 * t63) * t131;
t35 = -t159 * t56 + t47;
t34 = t127 * (rSges(4,1) * t106 + rSges(4,2) * t105 + rSges(4,3) * t159) + t129 * (rSges(4,1) * t108 + rSges(4,2) * t107 + rSges(4,3) * t156) + t152;
t33 = -t132 * t58 + (-t117 * t60 + t118 * t62) * t131;
t32 = -t132 * t57 + (-t117 * t59 + t118 * t61) * t131;
t29 = t100 * t62 + t156 * t58 + t60 * t99;
t28 = t100 * t61 + t156 * t57 + t59 * t99;
t27 = t159 * t58 + t60 * t97 + t62 * t98;
t26 = t159 * t57 + t59 * t97 + t61 * t98;
t25 = -t132 * t165 + t156 * t164;
t24 = t132 * t43 + t159 * t74 + t39;
t19 = t127 * t63 + t129 * t64 + t144;
t18 = t47 + (-t127 * t165 + t129 * t43) * t131;
t16 = t165 * t129 + (t43 + t55) * t127 + t144;
t15 = t127 * t29 - t129 * t28;
t14 = t127 * t27 - t129 * t26;
t8 = -(t100 * t85 + t84 * t99) * t132 + (t28 * t127 + (t29 - t162) * t129) * t131;
t7 = -(t84 * t97 + t85 * t98) * t132 + (t27 * t129 + (t26 - t162) * t127) * t131;
t1 = [m(2) + m(3) - t148; m(3) * t65 + m(4) * t34 + m(5) * t19 + m(6) * t16; m(6) * (t16 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t19 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(4) * (t34 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(3) * (t113 ^ 2 * t149 + t65 ^ 2) + (-t123 * t91 - t12 - t14 + (t105 * t70 + t106 * t72 + t68 * t159) * t129) * t129 + (t13 + t15 + t122 * t92 + (t107 * t71 + t108 * t73 + t69 * t156) * t127 + (-t105 * t71 - t106 * t73 - t107 * t70 - t108 * t72 - t127 * t91 + t129 * t92 - t156 * t68 - t159 * t69) * t129) * t127; t148 * t132; m(6) * (-t132 * t16 + (t127 * t37 + t129 * t38) * t131) + m(5) * (-t132 * t19 + (t127 * t45 + t129 * t46) * t131) + m(4) * (-t132 * t34 + (t127 * t75 + t129 * t76) * t131); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t131 ^ 2 * t149 + t171); m(5) * t36 + m(6) * t18; (t127 * t33 - t129 * t32) * t168 + t8 * t170 + t7 * t169 + (t14 * t170 + t129 * t15 / 0.2e1) * t131 + m(6) * (t16 * t18 + t24 * t38 + t25 * t37) + m(5) * (t19 * t36 + t41 * t46 + t42 * t45) + t142; m(5) * (-t132 * t36 + (t127 * t42 + t129 * t41) * t131) + m(6) * (-t132 * t18 + (t127 * t25 + t129 * t24) * t131); t7 * t159 + t8 * t156 + m(6) * (t18 ^ 2 + t24 ^ 2 + t25 ^ 2) - t132 * (t171 * t83 + (t33 * t129 + t32 * t127 - (-t117 * t84 + t118 * t85) * t132) * t131) + m(5) * (t36 ^ 2 + t41 ^ 2 + t42 ^ 2) + t143; m(6) * t35; m(6) * (t16 * t35 + t37 * t40 + t38 * t39) + t142; m(6) * (-t132 * t35 + (t127 * t40 + t129 * t39) * t131); m(6) * (t18 * t35 + t24 * t39 + t25 * t40) + t143; m(6) * (t35 ^ 2 + t39 ^ 2 + t40 ^ 2) + t143;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
