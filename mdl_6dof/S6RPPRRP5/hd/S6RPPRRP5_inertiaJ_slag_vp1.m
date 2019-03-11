% Calculate joint inertia matrix for
% S6RPPRRP5
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:42
% EndTime: 2019-03-09 02:07:46
% DurationCPUTime: 1.80s
% Computational Cost: add. (2189->310), mult. (4898->441), div. (0->0), fcn. (5121->6), ass. (0->142)
t126 = -qJ(6) - pkin(8);
t187 = rSges(7,3) - t126;
t131 = cos(qJ(4));
t189 = Icges(5,5) * t131;
t188 = t189 / 0.2e1;
t128 = sin(qJ(4));
t130 = cos(qJ(5));
t127 = sin(qJ(5));
t71 = Icges(7,3) * t128 + (Icges(7,5) * t130 - Icges(7,6) * t127) * t131;
t72 = Icges(6,3) * t128 + (Icges(6,5) * t130 - Icges(6,6) * t127) * t131;
t79 = Icges(7,5) * t128 + (Icges(7,1) * t130 - Icges(7,4) * t127) * t131;
t80 = Icges(6,5) * t128 + (Icges(6,1) * t130 - Icges(6,4) * t127) * t131;
t186 = (t79 + t80) * t130 * t131 + (t71 + t72) * t128;
t75 = Icges(7,6) * t128 + (Icges(7,4) * t130 - Icges(7,2) * t127) * t131;
t76 = Icges(6,6) * t128 + (Icges(6,4) * t130 - Icges(6,2) * t127) * t131;
t185 = (-t75 - t76) * t127;
t129 = sin(qJ(1));
t156 = t129 * t131;
t132 = cos(qJ(1));
t154 = t130 * t132;
t158 = t129 * t127;
t91 = -t128 * t158 + t154;
t157 = t129 * t130;
t160 = t127 * t132;
t92 = t128 * t157 + t160;
t42 = Icges(7,5) * t92 + Icges(7,6) * t91 - Icges(7,3) * t156;
t46 = Icges(7,4) * t92 + Icges(7,2) * t91 - Icges(7,6) * t156;
t50 = Icges(7,1) * t92 + Icges(7,4) * t91 - Icges(7,5) * t156;
t11 = -t42 * t156 + t46 * t91 + t50 * t92;
t153 = t131 * t132;
t159 = t128 * t132;
t93 = -t127 * t159 - t157;
t94 = t128 * t154 - t158;
t43 = Icges(7,5) * t94 + Icges(7,6) * t93 - Icges(7,3) * t153;
t47 = Icges(7,4) * t94 + Icges(7,2) * t93 - Icges(7,6) * t153;
t51 = Icges(7,1) * t94 + Icges(7,4) * t93 - Icges(7,5) * t153;
t12 = -t43 * t156 + t47 * t91 + t51 * t92;
t44 = Icges(6,5) * t92 + Icges(6,6) * t91 - Icges(6,3) * t156;
t48 = Icges(6,4) * t92 + Icges(6,2) * t91 - Icges(6,6) * t156;
t52 = Icges(6,1) * t92 + Icges(6,4) * t91 - Icges(6,5) * t156;
t13 = -t44 * t156 + t48 * t91 + t52 * t92;
t45 = Icges(6,5) * t94 + Icges(6,6) * t93 - Icges(6,3) * t153;
t49 = Icges(6,4) * t94 + Icges(6,2) * t93 - Icges(6,6) * t153;
t53 = Icges(6,1) * t94 + Icges(6,4) * t93 - Icges(6,5) * t153;
t14 = -t45 * t156 + t49 * t91 + t53 * t92;
t26 = -t71 * t156 + t75 * t91 + t79 * t92;
t27 = -t72 * t156 + t76 * t91 + t80 * t92;
t184 = ((-t12 - t14) * t132 + (-t11 - t13) * t129) * t131 + (t26 + t27) * t128;
t15 = -t42 * t153 + t93 * t46 + t94 * t50;
t16 = -t43 * t153 + t93 * t47 + t94 * t51;
t17 = -t44 * t153 + t93 * t48 + t94 * t52;
t18 = -t45 * t153 + t93 * t49 + t94 * t53;
t28 = -t71 * t153 + t93 * t75 + t94 * t79;
t29 = -t72 * t153 + t93 * t76 + t94 * t80;
t183 = ((-t16 - t18) * t132 + (-t15 - t17) * t129) * t131 + (t28 + t29) * t128;
t21 = t128 * t42 + (-t127 * t46 + t130 * t50) * t131;
t23 = t128 * t44 + (-t127 * t48 + t130 * t52) * t131;
t182 = t21 + t23;
t22 = t128 * t43 + (-t127 * t47 + t130 * t51) * t131;
t24 = t128 * t45 + (-t127 * t49 + t130 * t53) * t131;
t181 = -t22 - t24;
t118 = pkin(5) * t130 + pkin(4);
t180 = t94 * rSges(7,1) + t93 * rSges(7,2) + t118 * t159 - t153 * t187;
t123 = t129 ^ 2;
t125 = t132 ^ 2;
t179 = t128 / 0.2e1;
t177 = t132 / 0.2e1;
t176 = -rSges(5,3) - pkin(7);
t109 = rSges(5,1) * t131 - rSges(5,2) * t128;
t175 = m(5) * t109;
t113 = t123 + t125;
t101 = m(5) * t113;
t174 = m(7) * t131;
t173 = pkin(4) * t128;
t172 = -pkin(1) - qJ(3);
t171 = -pkin(4) + t118;
t170 = (t131 * t185 + t186) * t128;
t116 = pkin(8) * t156;
t144 = -t92 * rSges(7,1) - t91 * rSges(7,2);
t169 = -pkin(5) * t160 - t116 - (t126 * t131 + t171 * t128) * t129 + rSges(7,3) * t156 + t144;
t117 = pkin(4) * t159;
t96 = -pkin(8) * t153 + t117;
t168 = -pkin(5) * t158 + t180 - t96;
t165 = (rSges(7,1) * t130 - rSges(7,2) * t127 + t171) * t131 + (-pkin(8) + t187) * t128;
t164 = t94 * rSges(6,1) + t93 * rSges(6,2);
t163 = Icges(5,4) * t128;
t151 = rSges(5,1) * t159 + rSges(5,2) * t153;
t150 = t132 * pkin(1) + t129 * qJ(2);
t149 = t132 * qJ(3) + t150;
t148 = -pkin(5) * t127 - pkin(7);
t147 = t165 * t132;
t146 = t101 + (m(4) + m(6) + m(7)) * t113;
t145 = -t92 * rSges(6,1) - t91 * rSges(6,2);
t143 = -rSges(5,1) * t128 - rSges(5,2) * t131;
t19 = t169 * t128 - t165 * t156;
t20 = t168 * t128 + t131 * t147;
t140 = t20 * t129 + t132 * t19;
t121 = t132 * qJ(2);
t31 = t121 + t148 * t132 + (-t118 * t128 + t187 * t131 + t172) * t129 + t144;
t32 = t148 * t129 + t149 + t180;
t139 = t129 * t32 + t132 * t31;
t111 = pkin(4) * t131 + pkin(8) * t128;
t97 = t129 * t111;
t39 = t165 * t129 + t97;
t98 = t132 * t111;
t40 = t98 + t147;
t138 = t39 * t129 + t132 * t40;
t136 = Icges(5,2) * t131 + t163;
t135 = Icges(5,5) * t128 + Icges(5,6) * t131;
t134 = -t27 / 0.2e1 - t26 / 0.2e1 - t23 / 0.2e1 - t21 / 0.2e1;
t133 = -t29 / 0.2e1 - t28 / 0.2e1 - t24 / 0.2e1 - t22 / 0.2e1;
t110 = rSges(2,1) * t132 - t129 * rSges(2,2);
t108 = -t129 * rSges(2,1) - rSges(2,2) * t132;
t105 = -Icges(5,6) * t128 + t189;
t95 = t129 * t173 - t116;
t86 = -rSges(3,2) * t132 + t129 * rSges(3,3) + t150;
t85 = rSges(3,3) * t132 + t121 + (rSges(3,2) - pkin(1)) * t129;
t84 = t128 * rSges(6,3) + (rSges(6,1) * t130 - rSges(6,2) * t127) * t131;
t74 = -Icges(5,3) * t129 + t135 * t132;
t73 = Icges(5,3) * t132 + t135 * t129;
t70 = t129 * rSges(4,2) + rSges(4,3) * t132 + t149;
t69 = rSges(4,2) * t132 + t121 + (-rSges(4,3) + t172) * t129;
t63 = t132 * t84 + t98;
t62 = t129 * t84 + t97;
t61 = t176 * t129 + t149 + t151;
t60 = t121 + t176 * t132 + (t143 + t172) * t129;
t59 = -rSges(6,3) * t153 + t164;
t57 = -rSges(6,3) * t156 - t145;
t41 = t143 * t123 - t132 * t151;
t38 = t128 * t59 + t84 * t153;
t37 = -t128 * t57 - t84 * t156;
t36 = -t129 * pkin(7) + t117 + (-rSges(6,3) - pkin(8)) * t153 + t149 + t164;
t35 = -pkin(7) * t132 + t116 + t121 + (rSges(6,3) * t131 + t172 - t173) * t129 + t145;
t30 = (t129 * t59 - t132 * t57) * t131;
t25 = (-t59 - t96) * t132 + (-t57 - t95) * t129;
t10 = (t168 * t129 + t169 * t132) * t131;
t9 = (-t96 - t168) * t132 + (-t95 + t169) * t129;
t8 = -t18 * t129 + t132 * t17;
t7 = -t16 * t129 + t132 * t15;
t6 = -t14 * t129 + t13 * t132;
t5 = t11 * t132 - t12 * t129;
t1 = [-t128 * (Icges(5,4) * t131 - Icges(5,2) * t128) + Icges(3,1) + Icges(4,1) + Icges(2,3) + (Icges(5,1) * t131 - t163 + t185) * t131 + m(6) * (t35 ^ 2 + t36 ^ 2) + m(7) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t69 ^ 2 + t70 ^ 2) + m(2) * (t108 ^ 2 + t110 ^ 2) + t186; m(6) * (t129 * t35 - t132 * t36) + m(7) * (t129 * t31 - t132 * t32) + m(5) * (t129 * t60 - t132 * t61) + m(3) * (t129 * t85 - t132 * t86) + m(4) * (t129 * t69 - t132 * t70); m(3) * t113 + t146; m(6) * (t129 * t36 + t132 * t35) + m(7) * t139 + m(5) * (t129 * t61 + t132 * t60) + m(4) * (t129 * t70 + t132 * t69); 0; t146; m(6) * (t35 * t63 + t36 * t62) + m(7) * (t31 * t40 + t32 * t39) + (-t128 * (Icges(5,6) * t132 + t136 * t129) / 0.2e1 + t132 * t188 + t60 * t175 + t105 * t177 - t134) * t132 + (t136 * t132 * t179 + t61 * t175 + t133 + (-Icges(5,6) * t179 + t188 + t105 / 0.2e1) * t129) * t129; m(6) * (t63 * t129 - t132 * t62) + m(7) * (t40 * t129 - t132 * t39); m(6) * (t62 * t129 + t132 * t63) + m(7) * t138 + t109 * t101; m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t113 * t109 ^ 2 + t41 ^ 2) + (-t123 * t74 - t7 - t8) * t129 + (t125 * t73 + t5 + t6 + (t129 * t73 - t132 * t74) * t129) * t132; m(6) * (t35 * t37 + t36 * t38) + m(7) * (t19 * t31 + t20 * t32) + (t134 * t129 + t133 * t132) * t131 + t170; m(6) * (t37 * t129 - t132 * t38) + m(7) * (t19 * t129 - t132 * t20); m(6) * (t38 * t129 + t132 * t37) + m(7) * t140; m(7) * (t10 * t9 + t19 * t40 + t20 * t39) + m(6) * (t25 * t30 + t37 * t63 + t38 * t62) + ((-t8 / 0.2e1 - t7 / 0.2e1) * t132 + (-t6 / 0.2e1 - t5 / 0.2e1) * t129) * t131 + (t181 * t129 + t182 * t132) * t179 - t183 * t129 / 0.2e1 + t184 * t177; t170 * t128 + m(7) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(6) * (t30 ^ 2 + t37 ^ 2 + t38 ^ 2) + (-t183 * t132 - t184 * t129 + (-t182 * t129 + t181 * t132) * t128) * t131; -t139 * t174; 0; -t113 * t174; m(7) * (t128 * t9 - t138 * t131); m(7) * (t128 * t10 - t140 * t131); m(7) * (t113 * t131 ^ 2 + t128 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
