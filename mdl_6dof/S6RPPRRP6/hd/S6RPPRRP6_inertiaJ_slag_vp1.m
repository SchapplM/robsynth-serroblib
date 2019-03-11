% Calculate joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:10:22
% EndTime: 2019-03-09 02:10:25
% DurationCPUTime: 1.75s
% Computational Cost: add. (2046->301), mult. (4896->443), div. (0->0), fcn. (5204->6), ass. (0->137)
t181 = rSges(7,1) + pkin(5);
t180 = rSges(7,3) + qJ(6);
t129 = cos(qJ(4));
t183 = Icges(5,5) * t129;
t182 = t183 / 0.2e1;
t126 = sin(qJ(4));
t128 = cos(qJ(5));
t125 = sin(qJ(5));
t156 = t125 * t129;
t71 = Icges(7,6) * t126 + (Icges(7,5) * t128 + Icges(7,3) * t125) * t129;
t72 = Icges(6,3) * t126 + (Icges(6,5) * t128 - Icges(6,6) * t125) * t129;
t75 = Icges(7,2) * t126 + (Icges(7,4) * t128 + Icges(7,6) * t125) * t129;
t79 = Icges(7,4) * t126 + (Icges(7,1) * t128 + Icges(7,5) * t125) * t129;
t80 = Icges(6,5) * t126 + (Icges(6,1) * t128 - Icges(6,4) * t125) * t129;
t179 = t71 * t156 + (t79 + t80) * t128 * t129 + (t72 + t75) * t126;
t130 = cos(qJ(1));
t149 = t130 * t128;
t127 = sin(qJ(1));
t154 = t127 * t125;
t92 = t126 * t154 - t149;
t153 = t127 * t128;
t93 = t125 * t130 + t126 * t153;
t178 = -t180 * t92 - t181 * t93;
t152 = t127 * t129;
t42 = Icges(7,5) * t93 - Icges(7,6) * t152 + Icges(7,3) * t92;
t46 = Icges(7,4) * t93 - Icges(7,2) * t152 + Icges(7,6) * t92;
t50 = Icges(7,1) * t93 - Icges(7,4) * t152 + Icges(7,5) * t92;
t11 = -t46 * t152 + t42 * t92 + t50 * t93;
t150 = t129 * t130;
t155 = t126 * t130;
t94 = t125 * t155 + t153;
t95 = t126 * t149 - t154;
t43 = Icges(7,5) * t95 - Icges(7,6) * t150 + Icges(7,3) * t94;
t47 = Icges(7,4) * t95 - Icges(7,2) * t150 + Icges(7,6) * t94;
t51 = Icges(7,1) * t95 - Icges(7,4) * t150 + Icges(7,5) * t94;
t12 = -t47 * t152 + t43 * t92 + t51 * t93;
t44 = Icges(6,5) * t93 - Icges(6,6) * t92 - Icges(6,3) * t152;
t48 = Icges(6,4) * t93 - Icges(6,2) * t92 - Icges(6,6) * t152;
t52 = Icges(6,1) * t93 - Icges(6,4) * t92 - Icges(6,5) * t152;
t13 = -t44 * t152 - t48 * t92 + t52 * t93;
t45 = Icges(6,5) * t95 - Icges(6,6) * t94 - Icges(6,3) * t150;
t49 = Icges(6,4) * t95 - Icges(6,2) * t94 - Icges(6,6) * t150;
t53 = Icges(6,1) * t95 - Icges(6,4) * t94 - Icges(6,5) * t150;
t14 = -t45 * t152 - t49 * t92 + t53 * t93;
t28 = -t75 * t152 + t71 * t92 + t79 * t93;
t76 = Icges(6,6) * t126 + (Icges(6,4) * t128 - Icges(6,2) * t125) * t129;
t29 = -t72 * t152 - t76 * t92 + t80 * t93;
t177 = ((-t12 - t14) * t130 + (-t11 - t13) * t127) * t129 + (t28 + t29) * t126;
t15 = -t46 * t150 + t94 * t42 + t95 * t50;
t16 = -t47 * t150 + t94 * t43 + t95 * t51;
t17 = -t44 * t150 - t94 * t48 + t95 * t52;
t18 = -t45 * t150 - t94 * t49 + t95 * t53;
t30 = -t75 * t150 + t94 * t71 + t95 * t79;
t31 = -t72 * t150 - t94 * t76 + t95 * t80;
t176 = ((-t16 - t18) * t130 + (-t15 - t17) * t127) * t129 + (t30 + t31) * t126;
t19 = t126 * t46 + (t125 * t42 + t128 * t50) * t129;
t21 = t126 * t44 + (-t125 * t48 + t128 * t52) * t129;
t175 = t19 + t21;
t20 = t126 * t47 + (t125 * t43 + t128 * t51) * t129;
t22 = t126 * t45 + (-t125 * t49 + t128 * t53) * t129;
t174 = -t20 - t22;
t173 = t180 * t94 + t181 * t95;
t123 = t127 ^ 2;
t124 = t130 ^ 2;
t172 = t126 / 0.2e1;
t170 = t130 / 0.2e1;
t169 = -rSges(5,3) - pkin(7);
t110 = rSges(5,1) * t129 - rSges(5,2) * t126;
t168 = m(5) * t110;
t113 = t123 + t124;
t103 = m(5) * t113;
t167 = pkin(4) * t126;
t166 = -pkin(1) - qJ(3);
t165 = (-t76 * t156 + t179) * t126;
t164 = rSges(7,2) * t152 + t178;
t163 = -rSges(7,2) * t150 + t173;
t161 = t126 * rSges(7,2) + (t180 * t125 + t128 * t181) * t129;
t159 = t95 * rSges(6,1) - t94 * rSges(6,2);
t158 = Icges(5,4) * t126;
t148 = rSges(5,1) * t155 + rSges(5,2) * t150;
t147 = t130 * pkin(1) + t127 * qJ(2);
t145 = t130 * qJ(3) + t147;
t144 = t161 * t130;
t143 = t103 + (m(4) + m(6) + m(7)) * t113;
t142 = t166 - t167;
t116 = pkin(8) * t152;
t121 = t130 * qJ(2);
t141 = -pkin(7) * t130 + t116 + t121;
t140 = -t93 * rSges(6,1) + t92 * rSges(6,2);
t139 = -rSges(5,1) * t126 - rSges(5,2) * t129;
t135 = Icges(5,2) * t129 + t158;
t134 = Icges(5,5) * t126 + Icges(5,6) * t129;
t133 = -t22 / 0.2e1 - t20 / 0.2e1 - t31 / 0.2e1 - t30 / 0.2e1;
t132 = -t29 / 0.2e1 - t28 / 0.2e1 - t21 / 0.2e1 - t19 / 0.2e1;
t117 = pkin(4) * t155;
t131 = -t127 * pkin(7) + t117 + t145;
t112 = pkin(4) * t129 + pkin(8) * t126;
t111 = rSges(2,1) * t130 - t127 * rSges(2,2);
t109 = -t127 * rSges(2,1) - rSges(2,2) * t130;
t106 = -Icges(5,6) * t126 + t183;
t100 = t130 * t112;
t99 = t127 * t112;
t98 = -pkin(8) * t150 + t117;
t97 = t127 * t167 - t116;
t86 = -rSges(3,2) * t130 + t127 * rSges(3,3) + t147;
t85 = rSges(3,3) * t130 + t121 + (rSges(3,2) - pkin(1)) * t127;
t84 = t126 * rSges(6,3) + (rSges(6,1) * t128 - rSges(6,2) * t125) * t129;
t74 = -Icges(5,3) * t127 + t134 * t130;
t73 = Icges(5,3) * t130 + t134 * t127;
t70 = t127 * rSges(4,2) + rSges(4,3) * t130 + t145;
t69 = rSges(4,2) * t130 + t121 + (-rSges(4,3) + t166) * t127;
t61 = t130 * t84 + t100;
t60 = t127 * t84 + t99;
t59 = t169 * t127 + t145 + t148;
t58 = t121 + t169 * t130 + (t139 + t166) * t127;
t57 = -rSges(6,3) * t150 + t159;
t55 = -rSges(6,3) * t152 - t140;
t41 = t139 * t123 - t130 * t148;
t40 = t100 + t144;
t39 = t161 * t127 + t99;
t38 = t126 * t57 + t84 * t150;
t37 = -t126 * t55 - t84 * t152;
t36 = (-rSges(6,3) - pkin(8)) * t150 + t131 + t159;
t35 = (rSges(6,3) * t129 + t142) * t127 + t140 + t141;
t32 = (t127 * t57 - t130 * t55) * t129;
t27 = (-rSges(7,2) - pkin(8)) * t150 + t131 + t173;
t26 = (rSges(7,2) * t129 + t142) * t127 + t141 + t178;
t25 = (-t57 - t98) * t130 + (-t55 - t97) * t127;
t24 = t163 * t126 + t129 * t144;
t23 = t164 * t126 - t161 * t152;
t10 = (t163 * t127 + t164 * t130) * t129;
t9 = (-t98 - t163) * t130 + (-t97 + t164) * t127;
t8 = -t18 * t127 + t130 * t17;
t7 = -t16 * t127 + t130 * t15;
t6 = -t14 * t127 + t13 * t130;
t5 = t11 * t130 - t12 * t127;
t1 = [-t126 * (Icges(5,4) * t129 - Icges(5,2) * t126) + Icges(3,1) + Icges(4,1) + Icges(2,3) + (Icges(5,1) * t129 - t125 * t76 - t158) * t129 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(3) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t69 ^ 2 + t70 ^ 2) + m(2) * (t109 ^ 2 + t111 ^ 2) + t179; m(7) * (t127 * t26 - t130 * t27) + m(6) * (t127 * t35 - t130 * t36) + m(5) * (t127 * t58 - t130 * t59) + m(3) * (t127 * t85 - t130 * t86) + m(4) * (t127 * t69 - t130 * t70); m(3) * t113 + t143; m(7) * (t127 * t27 + t130 * t26) + m(6) * (t127 * t36 + t130 * t35) + m(5) * (t127 * t59 + t130 * t58) + m(4) * (t127 * t70 + t130 * t69); 0; t143; m(7) * (t26 * t40 + t27 * t39) + m(6) * (t35 * t61 + t36 * t60) + (-t126 * (Icges(5,6) * t130 + t135 * t127) / 0.2e1 + t130 * t182 + t58 * t168 + t106 * t170 - t132) * t130 + (t135 * t130 * t172 + t59 * t168 + t133 + (-Icges(5,6) * t172 + t182 + t106 / 0.2e1) * t127) * t127; m(6) * (t61 * t127 - t130 * t60) + m(7) * (t40 * t127 - t130 * t39); m(6) * (t60 * t127 + t130 * t61) + m(7) * (t39 * t127 + t130 * t40) + t110 * t103; m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t113 * t110 ^ 2 + t41 ^ 2) + (-t123 * t74 - t7 - t8) * t127 + (t124 * t73 + t5 + t6 + (t127 * t73 - t130 * t74) * t127) * t130; m(7) * (t23 * t26 + t24 * t27) + m(6) * (t35 * t37 + t36 * t38) + (t132 * t127 + t133 * t130) * t129 + t165; m(6) * (t37 * t127 - t130 * t38) + m(7) * (t23 * t127 - t130 * t24); m(6) * (t38 * t127 + t130 * t37) + m(7) * (t24 * t127 + t130 * t23); m(7) * (t10 * t9 + t23 * t40 + t24 * t39) + m(6) * (t25 * t32 + t37 * t61 + t38 * t60) + ((-t8 / 0.2e1 - t7 / 0.2e1) * t130 + (-t6 / 0.2e1 - t5 / 0.2e1) * t127) * t129 + (t174 * t127 + t175 * t130) * t172 - t176 * t127 / 0.2e1 + t177 * t170; t165 * t126 + m(7) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t32 ^ 2 + t37 ^ 2 + t38 ^ 2) + (-t176 * t130 - t177 * t127 + (-t175 * t127 + t174 * t130) * t126) * t129; m(7) * (t26 * t94 + t27 * t92); m(7) * (t94 * t127 - t130 * t92); m(7) * (t92 * t127 + t130 * t94); m(7) * (t9 * t156 + t39 * t92 + t40 * t94); m(7) * (t10 * t156 + t23 * t94 + t24 * t92); m(7) * (t125 ^ 2 * t129 ^ 2 + t92 ^ 2 + t94 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
