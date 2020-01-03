% Calculate joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:53:32
% DurationCPUTime: 1.54s
% Computational Cost: add. (3414->271), mult. (4492->414), div. (0->0), fcn. (4844->8), ass. (0->133)
t115 = pkin(8) + qJ(3);
t112 = sin(t115);
t177 = Icges(4,5) * t112;
t176 = t177 / 0.2e1;
t172 = rSges(6,3) + qJ(5);
t174 = rSges(6,1) + pkin(4);
t113 = cos(t115);
t123 = cos(qJ(4));
t124 = cos(qJ(1));
t142 = t123 * t124;
t121 = sin(qJ(4));
t122 = sin(qJ(1));
t144 = t122 * t121;
t91 = t113 * t144 + t142;
t143 = t122 * t123;
t145 = t121 * t124;
t92 = t113 * t143 - t145;
t175 = -t172 * t91 - t174 * t92;
t70 = -Icges(5,3) * t113 + (Icges(5,5) * t123 - Icges(5,6) * t121) * t112;
t71 = -Icges(6,2) * t113 + (Icges(6,4) * t123 + Icges(6,6) * t121) * t112;
t173 = -t70 - t71;
t150 = t112 * t122;
t42 = Icges(6,5) * t92 + Icges(6,6) * t150 + Icges(6,3) * t91;
t46 = Icges(6,4) * t92 + Icges(6,2) * t150 + Icges(6,6) * t91;
t50 = Icges(6,1) * t92 + Icges(6,4) * t150 + Icges(6,5) * t91;
t11 = t150 * t46 + t42 * t91 + t50 * t92;
t148 = t112 * t124;
t93 = t113 * t145 - t143;
t94 = t113 * t142 + t144;
t43 = Icges(6,5) * t94 + Icges(6,6) * t148 + Icges(6,3) * t93;
t47 = Icges(6,4) * t94 + Icges(6,2) * t148 + Icges(6,6) * t93;
t51 = Icges(6,1) * t94 + Icges(6,4) * t148 + Icges(6,5) * t93;
t12 = t150 * t47 + t43 * t91 + t51 * t92;
t44 = Icges(5,5) * t92 - Icges(5,6) * t91 + Icges(5,3) * t150;
t48 = Icges(5,4) * t92 - Icges(5,2) * t91 + Icges(5,6) * t150;
t52 = Icges(5,1) * t92 - Icges(5,4) * t91 + Icges(5,5) * t150;
t13 = t150 * t44 - t48 * t91 + t52 * t92;
t45 = Icges(5,5) * t94 - Icges(5,6) * t93 + Icges(5,3) * t148;
t49 = Icges(5,4) * t94 - Icges(5,2) * t93 + Icges(5,6) * t148;
t53 = Icges(5,1) * t94 - Icges(5,4) * t93 + Icges(5,5) * t148;
t14 = t150 * t45 - t49 * t91 + t53 * t92;
t69 = -Icges(6,6) * t113 + (Icges(6,5) * t123 + Icges(6,3) * t121) * t112;
t73 = -Icges(6,4) * t113 + (Icges(6,1) * t123 + Icges(6,5) * t121) * t112;
t26 = t150 * t71 + t69 * t91 + t73 * t92;
t72 = -Icges(5,6) * t113 + (Icges(5,4) * t123 - Icges(5,2) * t121) * t112;
t74 = -Icges(5,5) * t113 + (Icges(5,1) * t123 - Icges(5,4) * t121) * t112;
t27 = t150 * t70 - t72 * t91 + t74 * t92;
t171 = (-t26 - t27) * t113 + ((t12 + t14) * t124 + (t11 + t13) * t122) * t112;
t15 = t148 * t46 + t93 * t42 + t94 * t50;
t16 = t148 * t47 + t93 * t43 + t94 * t51;
t17 = t148 * t44 - t93 * t48 + t94 * t52;
t18 = t148 * t45 - t93 * t49 + t94 * t53;
t28 = t148 * t71 + t93 * t69 + t94 * t73;
t29 = t148 * t70 - t93 * t72 + t94 * t74;
t170 = (-t28 - t29) * t113 + ((t16 + t18) * t124 + (t15 + t17) * t122) * t112;
t19 = -t113 * t46 + (t121 * t42 + t123 * t50) * t112;
t21 = -t113 * t44 + (-t121 * t48 + t123 * t52) * t112;
t169 = -t19 - t21;
t20 = -t113 * t47 + (t121 * t43 + t123 * t51) * t112;
t22 = -t113 * t45 + (-t121 * t49 + t123 * t53) * t112;
t168 = t20 + t22;
t151 = t112 * t121;
t167 = t69 * t151 + (t73 + t74) * t112 * t123;
t116 = t122 ^ 2;
t117 = t124 ^ 2;
t99 = rSges(4,1) * t112 + rSges(4,2) * t113;
t166 = m(4) * t99;
t165 = -t113 / 0.2e1;
t164 = t122 / 0.2e1;
t162 = pkin(3) * t113;
t161 = t173 * t113 - t151 * t72 + t167;
t160 = rSges(6,2) * t150 - t175;
t159 = rSges(6,2) * t148 + t172 * t93 + t174 * t94;
t157 = -rSges(6,2) * t113 + (t172 * t121 + t174 * t123) * t112;
t147 = t113 * t124;
t141 = pkin(3) * t147 + pkin(7) * t148;
t156 = t116 * (pkin(7) * t112 + t162) + t124 * t141;
t155 = rSges(3,3) + qJ(2);
t100 = pkin(3) * t112 - pkin(7) * t113;
t76 = -rSges(5,3) * t113 + (rSges(5,1) * t123 - rSges(5,2) * t121) * t112;
t154 = -t100 - t76;
t152 = Icges(4,4) * t113;
t120 = -pkin(6) - qJ(2);
t146 = t120 * t124;
t140 = t116 + t117;
t139 = -t100 - t157;
t57 = t94 * rSges(5,1) - t93 * rSges(5,2) + rSges(5,3) * t148;
t119 = cos(pkin(8));
t110 = pkin(2) * t119 + pkin(1);
t138 = -t110 - t162;
t137 = t124 * t110 - t122 * t120;
t136 = -t92 * rSges(5,1) + t91 * rSges(5,2);
t135 = rSges(4,1) * t113 - rSges(4,2) * t112;
t131 = -Icges(4,2) * t112 + t152;
t130 = Icges(4,5) * t113 - Icges(4,6) * t112;
t129 = rSges(4,1) * t147 - rSges(4,2) * t148 + t122 * rSges(4,3);
t118 = sin(pkin(8));
t128 = rSges(3,1) * t119 - rSges(3,2) * t118 + pkin(1);
t127 = t137 + t141;
t126 = t20 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1 + t22 / 0.2e1;
t125 = t27 / 0.2e1 + t19 / 0.2e1 + t26 / 0.2e1 + t21 / 0.2e1;
t102 = rSges(2,1) * t124 - t122 * rSges(2,2);
t101 = -t122 * rSges(2,1) - rSges(2,2) * t124;
t96 = Icges(4,6) * t113 + t177;
t78 = Icges(4,3) * t122 + t124 * t130;
t77 = -Icges(4,3) * t124 + t122 * t130;
t68 = t122 * t155 + t124 * t128;
t67 = -t122 * t128 + t124 * t155;
t61 = t129 + t137;
t60 = (rSges(4,3) - t120) * t124 + (-t110 - t135) * t122;
t59 = t154 * t124;
t58 = t154 * t122;
t55 = rSges(5,3) * t150 - t136;
t41 = t124 * t129 + (-t124 * rSges(4,3) + t122 * t135) * t122;
t40 = t139 * t124;
t39 = t139 * t122;
t38 = t127 + t57;
t37 = -t146 + ((-rSges(5,3) - pkin(7)) * t112 + t138) * t122 + t136;
t36 = -t113 * t57 - t148 * t76;
t35 = t113 * t55 + t150 * t76;
t32 = (-t122 * t57 + t124 * t55) * t112;
t31 = t127 + t159;
t30 = -t146 + ((-rSges(6,2) - pkin(7)) * t112 + t138) * t122 + t175;
t25 = t122 * t55 + t124 * t57 + t156;
t24 = -t113 * t159 - t148 * t157;
t23 = t113 * t160 + t150 * t157;
t10 = (-t122 * t159 + t124 * t160) * t112;
t9 = t122 * t160 + t124 * t159 + t156;
t8 = t18 * t122 - t124 * t17;
t7 = t16 * t122 - t124 * t15;
t6 = t14 * t122 - t124 * t13;
t5 = -t11 * t124 + t12 * t122;
t1 = [Icges(3,2) * t119 ^ 2 + Icges(2,3) + (Icges(3,1) * t118 + 0.2e1 * Icges(3,4) * t119) * t118 + (Icges(4,1) * t112 - t121 * t72 + t152) * t112 + (Icges(4,4) * t112 + Icges(4,2) * t113 + t173) * t113 + m(6) * (t30 ^ 2 + t31 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t60 ^ 2 + t61 ^ 2) + m(3) * (t67 ^ 2 + t68 ^ 2) + m(2) * (t101 ^ 2 + t102 ^ 2) + t167; m(6) * (t122 * t30 - t124 * t31) + m(5) * (t122 * t37 - t124 * t38) + m(4) * (t122 * t60 - t124 * t61) + m(3) * (t122 * t67 - t124 * t68); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t140; m(6) * (t30 * t40 + t31 * t39) + m(5) * (t37 * t59 + t38 * t58) + (t122 * t131 * t165 - t60 * t166 - t125 + (t176 - Icges(4,6) * t165 + t96 / 0.2e1) * t124) * t124 + (-t61 * t166 + t122 * t176 + t113 * (Icges(4,6) * t122 + t124 * t131) / 0.2e1 + t96 * t164 + t126) * t122; m(5) * (t59 * t122 - t124 * t58) + m(6) * (t40 * t122 - t124 * t39); m(5) * (t25 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(6) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(4) * (t140 * t99 ^ 2 + t41 ^ 2) + (t116 * t78 + t7 + t8) * t122 + (-t117 * t77 - t5 - t6 + (-t122 * t77 + t124 * t78) * t122) * t124; -t161 * t113 + m(6) * (t23 * t30 + t24 * t31) + m(5) * (t35 * t37 + t36 * t38) + (t122 * t125 + t124 * t126) * t112; m(5) * (t35 * t122 - t124 * t36) + m(6) * (t23 * t122 - t124 * t24); m(5) * (t25 * t32 + t35 * t59 + t36 * t58) + m(6) * (t10 * t9 + t23 * t40 + t24 * t39) + ((t7 / 0.2e1 + t8 / 0.2e1) * t124 + (t5 / 0.2e1 + t6 / 0.2e1) * t122) * t112 + (t168 * t122 + t169 * t124) * t165 + t170 * t164 - t171 * t124 / 0.2e1; m(6) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + t161 * t113 ^ 2 + (t170 * t124 + t171 * t122 + (t169 * t122 - t168 * t124) * t113) * t112; m(6) * (t30 * t93 + t31 * t91); m(6) * (t93 * t122 - t124 * t91); m(6) * (t151 * t9 + t39 * t91 + t40 * t93); m(6) * (t10 * t151 + t23 * t93 + t24 * t91); m(6) * (t112 ^ 2 * t121 ^ 2 + t91 ^ 2 + t93 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
