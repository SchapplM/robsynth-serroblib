% Calculate joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:17
% EndTime: 2019-03-09 02:15:21
% DurationCPUTime: 1.86s
% Computational Cost: add. (3726->300), mult. (4950->447), div. (0->0), fcn. (5258->8), ass. (0->142)
t128 = pkin(9) + qJ(4);
t123 = cos(t128);
t197 = Icges(5,5) * t123;
t122 = sin(t128);
t196 = Icges(5,6) * t122;
t193 = rSges(7,1) + pkin(5);
t192 = rSges(7,3) + qJ(6);
t135 = sin(qJ(1));
t136 = cos(qJ(5));
t160 = t135 * t136;
t134 = sin(qJ(5));
t137 = cos(qJ(1));
t162 = t134 * t137;
t97 = t122 * t162 + t160;
t159 = t137 * t136;
t161 = t135 * t134;
t99 = -t122 * t159 + t161;
t195 = -t192 * t97 + t193 * t99;
t194 = t197 / 0.2e1 - t196 / 0.2e1;
t166 = t123 * t134;
t71 = Icges(7,6) * t122 + (Icges(7,5) * t136 + Icges(7,3) * t134) * t123;
t72 = Icges(6,3) * t122 + (Icges(6,5) * t136 - Icges(6,6) * t134) * t123;
t73 = Icges(7,2) * t122 + (Icges(7,4) * t136 + Icges(7,6) * t134) * t123;
t75 = Icges(7,4) * t122 + (Icges(7,1) * t136 + Icges(7,5) * t134) * t123;
t76 = Icges(6,5) * t122 + (Icges(6,1) * t136 - Icges(6,4) * t134) * t123;
t191 = t71 * t166 + (t75 + t76) * t123 * t136 + (t72 + t73) * t122;
t165 = t123 * t135;
t95 = t122 * t161 - t159;
t96 = t122 * t160 + t162;
t42 = Icges(7,5) * t96 - Icges(7,6) * t165 + Icges(7,3) * t95;
t46 = Icges(7,4) * t96 - Icges(7,2) * t165 + Icges(7,6) * t95;
t50 = Icges(7,1) * t96 - Icges(7,4) * t165 + Icges(7,5) * t95;
t11 = -t46 * t165 + t42 * t95 + t50 * t96;
t163 = t123 * t137;
t43 = Icges(7,5) * t99 + Icges(7,6) * t163 - Icges(7,3) * t97;
t47 = Icges(7,4) * t99 + Icges(7,2) * t163 - Icges(7,6) * t97;
t51 = Icges(7,1) * t99 + Icges(7,4) * t163 - Icges(7,5) * t97;
t12 = -t47 * t165 + t43 * t95 + t51 * t96;
t44 = Icges(6,5) * t96 - Icges(6,6) * t95 - Icges(6,3) * t165;
t48 = Icges(6,4) * t96 - Icges(6,2) * t95 - Icges(6,6) * t165;
t52 = Icges(6,1) * t96 - Icges(6,4) * t95 - Icges(6,5) * t165;
t13 = -t44 * t165 - t48 * t95 + t52 * t96;
t45 = Icges(6,5) * t99 + Icges(6,6) * t97 + Icges(6,3) * t163;
t49 = Icges(6,4) * t99 + Icges(6,2) * t97 + Icges(6,6) * t163;
t53 = Icges(6,1) * t99 + Icges(6,4) * t97 + Icges(6,5) * t163;
t14 = -t45 * t165 - t49 * t95 + t53 * t96;
t28 = -t73 * t165 + t71 * t95 + t75 * t96;
t74 = Icges(6,6) * t122 + (Icges(6,4) * t136 - Icges(6,2) * t134) * t123;
t29 = -t72 * t165 - t74 * t95 + t76 * t96;
t190 = ((t12 + t14) * t137 + (-t11 - t13) * t135) * t123 + (t28 + t29) * t122;
t15 = t46 * t163 - t97 * t42 + t99 * t50;
t16 = t47 * t163 - t97 * t43 + t99 * t51;
t17 = t44 * t163 + t97 * t48 + t99 * t52;
t18 = t45 * t163 + t97 * t49 + t99 * t53;
t30 = t73 * t163 - t97 * t71 + t99 * t75;
t31 = t72 * t163 + t97 * t74 + t99 * t76;
t189 = ((t16 + t18) * t137 + (-t15 - t17) * t135) * t123 + (t30 + t31) * t122;
t19 = t122 * t46 + (t134 * t42 + t136 * t50) * t123;
t21 = t122 * t44 + (-t134 * t48 + t136 * t52) * t123;
t188 = t19 + t21;
t20 = t122 * t47 + (t134 * t43 + t136 * t51) * t123;
t22 = t122 * t45 + (-t134 * t49 + t136 * t53) * t123;
t187 = t20 + t22;
t186 = t192 * t95 + t193 * t96;
t185 = (rSges(5,1) * t122 + rSges(5,2) * t123) * t137;
t129 = t135 ^ 2;
t130 = t137 ^ 2;
t181 = t135 / 0.2e1;
t180 = t137 / 0.2e1;
t105 = rSges(5,1) * t123 - rSges(5,2) * t122;
t179 = m(5) * t105;
t117 = t129 + t130;
t109 = m(5) * t117;
t131 = sin(pkin(9));
t178 = pkin(3) * t131;
t177 = (-t74 * t166 + t191) * t122;
t176 = -rSges(7,2) * t165 + t186;
t175 = rSges(7,2) * t163 + t195;
t173 = rSges(7,2) * t122 + (t192 * t134 + t193 * t136) * t123;
t171 = t96 * rSges(6,1) - t95 * rSges(6,2);
t170 = rSges(4,3) + qJ(3);
t167 = t122 * t135;
t158 = t137 * pkin(1) + t135 * qJ(2);
t156 = rSges(5,1) * t167 + rSges(5,2) * t165 + t137 * rSges(5,3);
t125 = t137 * qJ(2);
t133 = -pkin(7) - qJ(3);
t155 = t135 * t133 + t137 * t178 + t125;
t154 = t123 * (-rSges(7,2) - pkin(8));
t153 = t123 * (-rSges(6,3) - pkin(8));
t152 = t173 * t135;
t151 = t109 + (m(4) + m(6) + m(7)) * t117;
t150 = -t99 * rSges(6,1) - t97 * rSges(6,2);
t132 = cos(pkin(9));
t149 = rSges(4,1) * t131 + rSges(4,2) * t132;
t143 = Icges(5,5) * t122 + Icges(5,6) * t123;
t142 = -t133 * t137 + t135 * t178 + t158;
t141 = -t29 / 0.2e1 - t28 / 0.2e1 - t21 / 0.2e1 - t19 / 0.2e1;
t140 = t31 / 0.2e1 + t30 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t116 = t137 * t122 * pkin(4);
t139 = -t135 * pkin(1) + t116 + t155;
t115 = pkin(4) * t167;
t138 = t115 + t142;
t112 = rSges(2,1) * t137 - t135 * rSges(2,2);
t111 = -t135 * rSges(2,1) - rSges(2,2) * t137;
t106 = pkin(4) * t123 + pkin(8) * t122;
t102 = -t196 + t197;
t100 = t135 * t106;
t94 = -rSges(3,2) * t137 + t135 * rSges(3,3) + t158;
t93 = rSges(3,3) * t137 + t125 + (rSges(3,2) - pkin(1)) * t135;
t92 = -pkin(8) * t165 + t115;
t85 = t137 * (pkin(8) * t163 - t116);
t80 = Icges(5,3) * t135 - t137 * t143;
t79 = Icges(5,3) * t137 + t135 * t143;
t78 = rSges(6,3) * t122 + (rSges(6,1) * t136 - rSges(6,2) * t134) * t123;
t68 = t149 * t135 + t170 * t137 + t158;
t67 = t125 + t149 * t137 + (-pkin(1) - t170) * t135;
t61 = t142 + t156;
t60 = t185 + (-rSges(5,3) - pkin(1)) * t135 + t155;
t59 = (-t106 - t78) * t137;
t58 = t135 * t78 + t100;
t57 = rSges(6,3) * t163 - t150;
t55 = -rSges(6,3) * t165 + t171;
t41 = -t135 * t156 + (t135 * rSges(5,3) - t185) * t137;
t40 = (-t106 - t173) * t137;
t39 = t100 + t152;
t38 = t135 * t153 + t138 + t171;
t37 = t137 * t153 + t139 + t150;
t36 = -t122 * t57 + t78 * t163;
t35 = t122 * t55 + t78 * t165;
t32 = (-t135 * t57 - t137 * t55) * t123;
t27 = t135 * t154 + t138 + t186;
t26 = t137 * t154 + t139 - t195;
t25 = t137 * t57 + t85 + (-t55 - t92) * t135;
t24 = -t175 * t122 + t173 * t163;
t23 = t176 * t122 + t123 * t152;
t10 = (-t175 * t135 - t176 * t137) * t123;
t9 = t85 + t175 * t137 + (-t92 - t176) * t135;
t8 = t18 * t135 + t137 * t17;
t7 = t16 * t135 + t137 * t15;
t6 = t13 * t137 + t14 * t135;
t5 = t11 * t137 + t12 * t135;
t1 = [Icges(4,1) * t132 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t132 + Icges(4,2) * t131) * t131 + (Icges(5,1) * t123 - t134 * t74) * t123 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t67 ^ 2 + t68 ^ 2) + m(3) * (t93 ^ 2 + t94 ^ 2) + m(2) * (t111 ^ 2 + t112 ^ 2) + t191 + (-0.2e1 * Icges(5,4) * t123 + Icges(5,2) * t122) * t122; m(7) * (t135 * t26 - t137 * t27) + m(6) * (t135 * t37 - t137 * t38) + m(5) * (t135 * t60 - t137 * t61) + m(4) * (t135 * t67 - t137 * t68) + m(3) * (t135 * t93 - t137 * t94); m(3) * t117 + t151; m(7) * (t135 * t27 + t137 * t26) + m(6) * (t135 * t38 + t137 * t37) + m(5) * (t135 * t61 + t137 * t60) + m(4) * (t135 * t68 + t137 * t67); 0; t151; m(7) * (t26 * t39 + t27 * t40) + m(6) * (t37 * t58 + t38 * t59) + (t102 * t180 + t137 * t194 - t61 * t179 - t141) * t137 + (t102 * t181 + t135 * t194 + t60 * t179 + t140) * t135; m(6) * (t58 * t135 - t137 * t59) + m(7) * (t39 * t135 - t137 * t40) + t105 * t109; m(6) * (t59 * t135 + t137 * t58) + m(7) * (t40 * t135 + t137 * t39); m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t105 ^ 2 * t117 + t41 ^ 2) + (t129 * t80 + t7 + t8) * t135 + (t130 * t79 + t5 + t6 + (t135 * t79 + t137 * t80) * t135) * t137; m(7) * (t23 * t27 + t24 * t26) + m(6) * (t35 * t38 + t36 * t37) + (t135 * t141 + t137 * t140) * t123 + t177; m(6) * (t36 * t135 - t137 * t35) + m(7) * (t24 * t135 - t137 * t23); m(6) * (t35 * t135 + t137 * t36) + m(7) * (t23 * t135 + t137 * t24); m(7) * (t10 * t9 + t23 * t40 + t24 * t39) + m(6) * (t25 * t32 + t35 * t59 + t36 * t58) + ((t8 / 0.2e1 + t7 / 0.2e1) * t137 + (-t6 / 0.2e1 - t5 / 0.2e1) * t135) * t123 + (t187 * t135 + t188 * t137) * t122 / 0.2e1 + t189 * t181 + t190 * t180; t177 * t122 + m(7) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + (t189 * t137 - t190 * t135 + (-t188 * t135 + t187 * t137) * t122) * t123; m(7) * (t26 * t95 - t27 * t97); m(7) * (t95 * t135 + t137 * t97); m(7) * (-t97 * t135 + t137 * t95); m(7) * (t9 * t166 + t39 * t95 - t40 * t97); m(7) * (t10 * t166 - t23 * t97 + t24 * t95); m(7) * (t123 ^ 2 * t134 ^ 2 + t95 ^ 2 + t97 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
