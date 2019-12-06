% Calculate joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:44
% EndTime: 2019-12-05 16:08:52
% DurationCPUTime: 2.45s
% Computational Cost: add. (3536->275), mult. (5545->401), div. (0->0), fcn. (6026->8), ass. (0->133)
t185 = -Icges(6,4) - Icges(5,5);
t184 = -Icges(6,6) + Icges(5,6);
t183 = -Icges(5,4) + Icges(6,5);
t123 = qJ(3) + pkin(8);
t119 = sin(t123);
t120 = cos(t123);
t128 = sin(qJ(3));
t129 = sin(qJ(2));
t130 = cos(qJ(3));
t131 = cos(qJ(2));
t182 = (Icges(5,3) + Icges(6,2) + Icges(4,3)) * t131 + (-Icges(4,5) * t130 + Icges(4,6) * t128 + t184 * t119 + t185 * t120) * t129;
t181 = t184 * t131 + (t183 * t120 + (Icges(5,2) + Icges(6,3)) * t119) * t129;
t180 = t185 * t131 + ((Icges(5,1) + Icges(6,1)) * t120 + t183 * t119) * t129;
t179 = t182 * t131;
t178 = rSges(6,1) + pkin(4);
t177 = rSges(6,3) + qJ(5);
t125 = sin(pkin(7));
t121 = t125 ^ 2;
t126 = cos(pkin(7));
t122 = t126 ^ 2;
t145 = t121 + t122;
t103 = -Icges(4,6) * t131 + (Icges(4,4) * t130 - Icges(4,2) * t128) * t129;
t104 = -Icges(4,5) * t131 + (Icges(4,1) * t130 - Icges(4,4) * t128) * t129;
t149 = t128 * t131;
t110 = -t125 * t149 - t126 * t130;
t148 = t130 * t131;
t152 = t126 * t128;
t111 = t125 * t148 - t152;
t154 = t125 * t129;
t153 = t125 * t131;
t98 = t119 * t153 + t126 * t120;
t99 = -t126 * t119 + t120 * t153;
t47 = Icges(6,5) * t99 + Icges(6,6) * t154 + Icges(6,3) * t98;
t51 = Icges(6,4) * t99 + Icges(6,2) * t154 + Icges(6,6) * t98;
t55 = Icges(6,1) * t99 + Icges(6,4) * t154 + Icges(6,5) * t98;
t19 = t51 * t154 + t47 * t98 + t55 * t99;
t150 = t126 * t131;
t100 = t119 * t150 - t125 * t120;
t101 = t125 * t119 + t120 * t150;
t151 = t126 * t129;
t48 = Icges(6,5) * t101 + Icges(6,6) * t151 + Icges(6,3) * t100;
t52 = Icges(6,4) * t101 + Icges(6,2) * t151 + Icges(6,6) * t100;
t56 = Icges(6,1) * t101 + Icges(6,4) * t151 + Icges(6,5) * t100;
t20 = t52 * t154 + t48 * t98 + t56 * t99;
t49 = Icges(5,5) * t99 - Icges(5,6) * t98 + Icges(5,3) * t154;
t53 = Icges(5,4) * t99 - Icges(5,2) * t98 + Icges(5,6) * t154;
t57 = Icges(5,1) * t99 - Icges(5,4) * t98 + Icges(5,5) * t154;
t21 = t49 * t154 - t53 * t98 + t57 * t99;
t50 = Icges(5,5) * t101 - Icges(5,6) * t100 + Icges(5,3) * t151;
t54 = Icges(5,4) * t101 - Icges(5,2) * t100 + Icges(5,6) * t151;
t58 = Icges(5,1) * t101 - Icges(5,4) * t100 + Icges(5,5) * t151;
t22 = t50 * t154 - t54 * t98 + t58 * t99;
t70 = Icges(4,5) * t111 + Icges(4,6) * t110 + Icges(4,3) * t154;
t72 = Icges(4,4) * t111 + Icges(4,2) * t110 + Icges(4,6) * t154;
t74 = Icges(4,1) * t111 + Icges(4,4) * t110 + Icges(4,5) * t154;
t31 = t110 * t72 + t111 * t74 + t70 * t154;
t112 = t125 * t130 - t126 * t149;
t155 = t125 * t128;
t113 = t126 * t148 + t155;
t71 = Icges(4,5) * t113 + Icges(4,6) * t112 + Icges(4,3) * t151;
t73 = Icges(4,4) * t113 + Icges(4,2) * t112 + Icges(4,6) * t151;
t75 = Icges(4,1) * t113 + Icges(4,4) * t112 + Icges(4,5) * t151;
t32 = t110 * t73 + t111 * t75 + t71 * t154;
t176 = (-t103 * t110 - t104 * t111 - t180 * t99 - t181 * t98) * t131 + ((t20 + t22 + t32) * t126 + (t19 + t21 + t31 + t179) * t125) * t129;
t23 = t100 * t47 + t101 * t55 + t51 * t151;
t24 = t100 * t48 + t101 * t56 + t52 * t151;
t25 = -t100 * t53 + t101 * t57 + t49 * t151;
t26 = -t100 * t54 + t101 * t58 + t50 * t151;
t33 = t112 * t72 + t113 * t74 + t70 * t151;
t34 = t112 * t73 + t113 * t75 + t71 * t151;
t175 = (-t181 * t100 - t180 * t101 - t103 * t112 - t104 * t113) * t131 + ((t24 + t26 + t34 + t179) * t126 + (t23 + t25 + t33) * t125) * t129;
t174 = (t51 + t49 + t70) * t131 + (t128 * t72 - t130 * t74 + (-t55 - t57) * t120 + (-t47 + t53) * t119) * t129;
t173 = (-t52 - t50 - t71) * t131 + (-t128 * t73 + t130 * t75 + (t56 + t58) * t120 + (t48 - t54) * t119) * t129;
t172 = t131 ^ 2;
t171 = -m(5) - m(6);
t167 = pkin(3) * t130;
t165 = rSges(6,2) * t154 + t177 * t98 + t178 * t99;
t164 = rSges(6,2) * t151 + t177 * t100 + t178 * t101;
t62 = rSges(5,1) * t101 - rSges(5,2) * t100 + rSges(5,3) * t151;
t132 = qJ(4) * t129 + t167 * t131;
t77 = pkin(3) * t155 + t132 * t126;
t163 = -t62 - t77;
t76 = -pkin(3) * t152 + t132 * t125;
t83 = -qJ(4) * t131 + t167 * t129;
t162 = t131 * t76 + t83 * t154;
t91 = -rSges(5,3) * t131 + (rSges(5,1) * t120 - rSges(5,2) * t119) * t129;
t161 = -t83 - t91;
t158 = -rSges(6,2) * t131 + (t177 * t119 + t178 * t120) * t129;
t156 = t119 * t129;
t106 = -rSges(4,3) * t131 + (rSges(4,1) * t130 - rSges(4,2) * t128) * t129;
t116 = t129 * pkin(2) - pkin(6) * t131;
t147 = -t106 - t116;
t146 = t145 * (pkin(2) * t131 + pkin(6) * t129);
t144 = -t77 - t164;
t143 = -t83 - t158;
t142 = -t116 + t161;
t141 = -t116 + t143;
t140 = t125 * t76 + t126 * t77 + t146;
t133 = Icges(3,5) * t131 - Icges(3,6) * t129;
t124 = t129 ^ 2;
t115 = t129 * rSges(3,1) + rSges(3,2) * t131;
t93 = Icges(3,3) * t125 + t133 * t126;
t92 = -Icges(3,3) * t126 + t133 * t125;
t81 = t147 * t126;
t80 = t147 * t125;
t79 = rSges(4,1) * t113 + rSges(4,2) * t112 + rSges(4,3) * t151;
t78 = rSges(4,1) * t111 + rSges(4,2) * t110 + rSges(4,3) * t154;
t64 = t145 * (rSges(3,1) * t131 - rSges(3,2) * t129);
t63 = t76 * t151;
t60 = rSges(5,1) * t99 - rSges(5,2) * t98 + rSges(5,3) * t154;
t46 = t142 * t126;
t45 = t142 * t125;
t44 = -t106 * t151 - t131 * t79;
t43 = t106 * t154 + t131 * t78;
t42 = t141 * t126;
t41 = t141 * t125;
t40 = (-t125 * t79 + t126 * t78) * t129;
t39 = t125 * t78 + t126 * t79 + t146;
t36 = t163 * t131 + t161 * t151;
t35 = t131 * t60 + t91 * t154 + t162;
t18 = t63 + (t163 * t125 + t126 * t60) * t129;
t17 = t144 * t131 + t143 * t151;
t16 = t165 * t131 + t158 * t154 + t162;
t15 = t125 * t60 + t126 * t62 + t140;
t14 = t63 + (t144 * t125 + t165 * t126) * t129;
t13 = t165 * t125 + t164 * t126 + t140;
t12 = t125 * t34 - t126 * t33;
t11 = t125 * t32 - t126 * t31;
t10 = t125 * t26 - t126 * t25;
t9 = t125 * t24 - t126 * t23;
t8 = t125 * t22 - t126 * t21;
t7 = t125 * t20 - t126 * t19;
t1 = [m(2) + m(3) + m(4) - t171; m(3) * t64 + m(4) * t39 + m(5) * t15 + m(6) * t13; m(6) * (t13 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t15 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(4) * (t39 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(3) * (t145 * t115 ^ 2 + t64 ^ 2) + (-t122 * t92 - t11 - t7 - t8) * t126 + (t121 * t93 + t10 + t12 + t9 + (-t125 * t92 + t126 * t93) * t126) * t125; m(4) * t40 + m(5) * t18 + m(6) * t14; m(6) * (t13 * t14 + t16 * t42 + t17 * t41) + m(5) * (t15 * t18 + t35 * t46 + t36 * t45) + m(4) * (t39 * t40 + t43 * t81 + t44 * t80) + ((t9 / 0.2e1 + t10 / 0.2e1 + t12 / 0.2e1) * t126 + (t8 / 0.2e1 + t11 / 0.2e1 + t7 / 0.2e1) * t125) * t129 + t175 * t125 / 0.2e1 - t176 * t126 / 0.2e1 - (t173 * t125 + t174 * t126) * t131 / 0.2e1; m(5) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(4) * (t40 ^ 2 + t43 ^ 2 + t44 ^ 2) + (t182 * t172 + ((-t103 * t128 + t104 * t130 + t181 * t119 + t180 * t120) * t131 - t173 * t126 + t174 * t125) * t129) * t131 + t176 * t154 + t175 * t151; t171 * t131; m(6) * (-t13 * t131 + (t125 * t41 + t126 * t42) * t129) + m(5) * (-t131 * t15 + (t125 * t45 + t126 * t46) * t129); m(5) * (-t131 * t18 + (t125 * t36 + t126 * t35) * t129) + m(6) * (-t131 * t14 + (t125 * t17 + t126 * t16) * t129); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t145 * t124 + t172); m(6) * t156; m(6) * (t100 * t42 + t13 * t156 + t41 * t98); m(6) * (t100 * t16 + t14 * t156 + t17 * t98); m(6) * (t100 * t126 - t119 * t131 + t125 * t98) * t129; m(6) * (t119 ^ 2 * t124 + t100 ^ 2 + t98 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
