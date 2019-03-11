% Calculate joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:34
% EndTime: 2019-03-09 01:41:38
% DurationCPUTime: 1.47s
% Computational Cost: add. (3500->248), mult. (2921->371), div. (0->0), fcn. (2959->10), ass. (0->124)
t172 = -Icges(6,4) + Icges(5,5);
t171 = Icges(6,5) - Icges(5,6);
t170 = Icges(6,1) + Icges(5,3);
t100 = pkin(10) + qJ(4);
t95 = sin(t100);
t97 = cos(t100);
t169 = t171 * t95 + t172 * t97;
t168 = (Icges(5,6) / 0.2e1 - Icges(6,5) / 0.2e1) * t97 + (Icges(5,5) / 0.2e1 - Icges(6,4) / 0.2e1) * t95;
t101 = qJ(1) + pkin(9);
t96 = sin(t101);
t167 = -t96 / 0.2e1;
t155 = t96 / 0.2e1;
t98 = cos(t101);
t166 = -t98 / 0.2e1;
t165 = t98 / 0.2e1;
t160 = -t169 * t96 + t170 * t98;
t159 = t169 * t98 + t170 * t96;
t150 = t97 * t98;
t107 = cos(qJ(6));
t133 = t107 * t98;
t105 = sin(qJ(6));
t136 = t105 * t96;
t62 = t133 * t95 - t136;
t134 = t107 * t96;
t135 = t105 * t98;
t63 = t135 * t95 + t134;
t32 = t63 * rSges(7,1) + t62 * rSges(7,2) + rSges(7,3) * t150;
t158 = t96 * pkin(5) + pkin(8) * t150 + t32;
t93 = t96 ^ 2;
t94 = t98 ^ 2;
t157 = m(6) / 0.2e1;
t156 = m(7) / 0.2e1;
t154 = -m(6) - m(7);
t106 = sin(qJ(1));
t153 = pkin(1) * t106;
t152 = t95 * t98;
t151 = t96 * t97;
t137 = qJ(5) * t95;
t147 = pkin(4) * t150 + t98 * t137;
t149 = t93 * (pkin(4) * t97 + t137) + t98 * t147;
t75 = pkin(4) * t95 - qJ(5) * t97;
t148 = rSges(6,2) * t95 + rSges(6,3) * t97 - t75;
t145 = t93 + t94;
t144 = Icges(5,4) * t95;
t143 = Icges(5,4) * t97;
t142 = Icges(7,5) * t97;
t141 = Icges(6,6) * t95;
t140 = Icges(6,6) * t97;
t139 = Icges(7,6) * t97;
t138 = Icges(7,3) * t97;
t132 = rSges(4,3) + qJ(3);
t131 = t157 + t156;
t64 = t134 * t95 + t135;
t65 = t136 * t95 - t133;
t27 = Icges(7,5) * t65 + Icges(7,6) * t64 + t138 * t96;
t29 = Icges(7,4) * t65 + Icges(7,2) * t64 + t139 * t96;
t31 = Icges(7,1) * t65 + Icges(7,4) * t64 + t142 * t96;
t11 = t27 * t95 + (-t105 * t31 - t107 * t29) * t97;
t56 = Icges(7,3) * t95 + (-Icges(7,5) * t105 - Icges(7,6) * t107) * t97;
t57 = Icges(7,6) * t95 + (-Icges(7,4) * t105 - Icges(7,2) * t107) * t97;
t58 = Icges(7,5) * t95 + (-Icges(7,1) * t105 - Icges(7,4) * t107) * t97;
t14 = t151 * t56 + t57 * t64 + t58 * t65;
t130 = t11 / 0.2e1 + t14 / 0.2e1;
t26 = Icges(7,5) * t63 + Icges(7,6) * t62 + t138 * t98;
t28 = Icges(7,4) * t63 + Icges(7,2) * t62 + t139 * t98;
t30 = Icges(7,1) * t63 + Icges(7,4) * t62 + t142 * t98;
t10 = t26 * t95 + (-t105 * t30 - t107 * t28) * t97;
t13 = t150 * t56 + t57 * t62 + t58 * t63;
t129 = t13 / 0.2e1 + t10 / 0.2e1;
t59 = t95 * rSges(7,3) + (-rSges(7,1) * t105 - rSges(7,2) * t107) * t97;
t128 = -pkin(8) * t95 - t59 - t75;
t104 = -pkin(7) - qJ(3);
t103 = cos(pkin(10));
t91 = pkin(3) * t103 + pkin(2);
t108 = cos(qJ(1));
t99 = t108 * pkin(1);
t126 = -t104 * t96 + t98 * t91 + t99;
t125 = rSges(5,1) * t97 - rSges(5,2) * t95;
t124 = -rSges(7,1) * t65 - rSges(7,2) * t64;
t119 = Icges(5,1) * t97 - t144;
t118 = -Icges(5,2) * t95 + t143;
t115 = -Icges(6,2) * t97 + t141;
t114 = Icges(6,3) * t95 - t140;
t113 = -t105 * t58 - t107 * t57;
t112 = rSges(5,1) * t150 - rSges(5,2) * t152 + t96 * rSges(5,3);
t111 = t96 * rSges(6,1) - rSges(6,2) * t150 + rSges(6,3) * t152;
t102 = sin(pkin(10));
t110 = rSges(4,1) * t103 - rSges(4,2) * t102 + pkin(2);
t109 = t126 + t147;
t87 = rSges(2,1) * t108 - t106 * rSges(2,2);
t86 = -t106 * rSges(2,1) - rSges(2,2) * t108;
t77 = rSges(5,1) * t95 + rSges(5,2) * t97;
t67 = rSges(3,1) * t98 - rSges(3,2) * t96 + t99;
t66 = -rSges(3,1) * t96 - rSges(3,2) * t98 - t153;
t40 = t95 * t56;
t39 = t148 * t98;
t38 = t148 * t96;
t37 = t110 * t98 + t132 * t96 + t99;
t36 = -t110 * t96 + t132 * t98 - t153;
t35 = t112 + t126;
t34 = -t153 + (rSges(5,3) - t104) * t98 + (-t125 - t91) * t96;
t33 = rSges(7,3) * t151 - t124;
t25 = t128 * t98;
t24 = t128 * t96;
t23 = t98 * t112 + (-t98 * rSges(5,3) + t125 * t96) * t96;
t22 = t109 + t111;
t21 = -t153 + (rSges(6,1) - t104) * t98 + (-t91 + (rSges(6,2) - pkin(4)) * t97 + (-rSges(6,3) - qJ(5)) * t95) * t96;
t20 = -t150 * t59 + t32 * t95;
t19 = t151 * t59 - t33 * t95;
t18 = (t113 * t97 + t40) * t95;
t17 = t109 + t158;
t16 = -t153 + (pkin(5) - t104) * t98 + (-t137 - t91 + (-rSges(7,3) - pkin(4) - pkin(8)) * t97) * t96 + t124;
t15 = t98 * t111 + (-t98 * rSges(6,1) + (-rSges(6,2) * t97 + rSges(6,3) * t95) * t96) * t96 + t149;
t12 = (-t32 * t96 + t33 * t98) * t97;
t9 = t151 * t27 + t29 * t64 + t31 * t65;
t8 = t151 * t26 + t28 * t64 + t30 * t65;
t7 = t150 * t27 + t29 * t62 + t31 * t63;
t6 = t150 * t26 + t28 * t62 + t30 * t63;
t5 = t158 * t98 + (-pkin(5) * t98 + pkin(8) * t151 + t33) * t96 + t149;
t4 = t8 * t96 - t9 * t98;
t3 = t6 * t96 - t7 * t98;
t2 = t14 * t95 + (t8 * t98 + t9 * t96) * t97;
t1 = t13 * t95 + (t6 * t98 + t7 * t96) * t97;
t41 = [Icges(4,2) * t103 ^ 2 + Icges(2,3) + Icges(3,3) + t40 + (Icges(4,1) * t102 + 0.2e1 * Icges(4,4) * t103) * t102 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t21 ^ 2 + t22 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t66 ^ 2 + t67 ^ 2) + m(2) * (t86 ^ 2 + t87 ^ 2) + (t113 + t141 + t144 + (Icges(5,2) + Icges(6,3)) * t97) * t97 + (t140 + t143 + (Icges(5,1) + Icges(6,2)) * t95) * t95; 0; m(3) + m(4) + m(5) - t154; m(7) * (t16 * t96 - t17 * t98) + m(6) * (t21 * t96 - t22 * t98) + m(5) * (t34 * t96 - t35 * t98) + m(4) * (t36 * t96 - t37 * t98); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t131) * t145; ((t114 * t155 + t118 * t167) * t97 + (t115 * t155 + t119 * t167) * t95 - t130 + t168 * t98) * t98 + ((t114 * t166 + t118 * t165) * t97 + (t115 * t166 + t119 * t165) * t95 + t129 + t168 * t96) * t96 + m(7) * (t16 * t25 + t17 * t24) + m(6) * (t21 * t39 + t22 * t38) + m(5) * (-t34 * t98 - t35 * t96) * t77 + (-t171 * t97 + t172 * t95) * (t94 / 0.2e1 + t93 / 0.2e1); m(5) * t23 + m(6) * t15 + m(7) * t5; m(6) * (-t38 * t98 + t39 * t96) + m(7) * (-t24 * t98 + t25 * t96); m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t15 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t145 * t77 ^ 2 + t23 ^ 2) + (t159 * t93 + t3) * t96 + (-t4 + t160 * t94 + (t159 * t98 + t160 * t96) * t96) * t98; 0.2e1 * ((t16 * t98 + t17 * t96) * t156 + (t21 * t98 + t22 * t96) * t157) * t95; t154 * t97; 0; m(7) * (-t5 * t97 + (t24 * t96 + t25 * t98) * t95) + m(6) * (-t15 * t97 + (t38 * t96 + t39 * t98) * t95); 0.2e1 * t131 * (t145 * t95 ^ 2 + t97 ^ 2); m(7) * (t16 * t19 + t17 * t20) + t18 + (t129 * t98 + t130 * t96) * t97; m(7) * t12; m(7) * (t19 * t96 - t20 * t98); t2 * t166 + t95 * (t10 * t96 - t11 * t98) / 0.2e1 + m(7) * (t12 * t5 + t19 * t25 + t20 * t24) + t1 * t155 + (t4 * t155 + t3 * t165) * t97; m(7) * (-t12 * t97 + (t19 * t98 + t20 * t96) * t95); t95 * t18 + m(7) * (t12 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t98 * t1 + t96 * t2 + t95 * (t10 * t98 + t11 * t96)) * t97;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t41(1) t41(2) t41(4) t41(7) t41(11) t41(16); t41(2) t41(3) t41(5) t41(8) t41(12) t41(17); t41(4) t41(5) t41(6) t41(9) t41(13) t41(18); t41(7) t41(8) t41(9) t41(10) t41(14) t41(19); t41(11) t41(12) t41(13) t41(14) t41(15) t41(20); t41(16) t41(17) t41(18) t41(19) t41(20) t41(21);];
Mq  = res;
