% Calculate joint inertia matrix for
% S5PRRPP3
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:01
% EndTime: 2019-12-05 16:12:07
% DurationCPUTime: 1.94s
% Computational Cost: add. (3014->277), mult. (7690->399), div. (0->0), fcn. (9152->8), ass. (0->128)
t176 = -Icges(5,4) + Icges(6,5);
t175 = Icges(6,4) + Icges(5,5);
t174 = -Icges(5,6) + Icges(6,6);
t122 = sin(pkin(8));
t126 = sin(qJ(2));
t127 = cos(qJ(3));
t128 = cos(qJ(2));
t150 = cos(pkin(8));
t111 = t126 * t127 * t122 + t128 * t150;
t138 = t126 * t150;
t112 = -t128 * t122 + t127 * t138;
t125 = sin(qJ(3));
t147 = t125 * t126;
t173 = t174 * t147 + t176 * t112 + (Icges(5,2) + Icges(6,3)) * t111;
t172 = (Icges(6,2) + Icges(5,3)) * t147 + t175 * t112 + t174 * t111;
t171 = t175 * t147 + (Icges(5,1) + Icges(6,1)) * t112 + t176 * t111;
t104 = -Icges(4,6) * t128 + (Icges(4,4) * t127 - Icges(4,2) * t125) * t126;
t170 = t104 - t172;
t169 = rSges(6,1) + pkin(4);
t168 = rSges(6,3) + qJ(5);
t124 = cos(pkin(7));
t161 = t124 ^ 2;
t123 = sin(pkin(7));
t162 = t123 ^ 2;
t167 = t161 + t162;
t144 = t128 * (-Icges(4,3) * t128 + (Icges(4,5) * t127 - Icges(4,6) * t125) * t126);
t105 = -Icges(4,5) * t128 + (Icges(4,1) * t127 - Icges(4,4) * t125) * t126;
t146 = t125 * t128;
t109 = t123 * t146 + t124 * t127;
t145 = t127 * t128;
t110 = t123 * t145 - t124 * t125;
t92 = t110 * t122 - t123 * t138;
t149 = t123 * t126;
t93 = t110 * t150 + t122 * t149;
t43 = Icges(6,5) * t93 + Icges(6,6) * t109 + Icges(6,3) * t92;
t47 = Icges(6,4) * t93 + Icges(6,2) * t109 + Icges(6,6) * t92;
t51 = Icges(6,1) * t93 + Icges(6,4) * t109 + Icges(6,5) * t92;
t15 = t109 * t47 + t92 * t43 + t93 * t51;
t113 = -t123 * t127 + t124 * t146;
t114 = t123 * t125 + t124 * t145;
t94 = t114 * t122 - t124 * t138;
t148 = t124 * t126;
t95 = t114 * t150 + t122 * t148;
t44 = Icges(6,5) * t95 + Icges(6,6) * t113 + Icges(6,3) * t94;
t48 = Icges(6,4) * t95 + Icges(6,2) * t113 + Icges(6,6) * t94;
t52 = Icges(6,1) * t95 + Icges(6,4) * t113 + Icges(6,5) * t94;
t16 = t109 * t48 + t92 * t44 + t93 * t52;
t45 = Icges(5,5) * t93 - Icges(5,6) * t92 + Icges(5,3) * t109;
t49 = Icges(5,4) * t93 - Icges(5,2) * t92 + Icges(5,6) * t109;
t53 = Icges(5,1) * t93 - Icges(5,4) * t92 + Icges(5,5) * t109;
t17 = t109 * t45 - t92 * t49 + t93 * t53;
t46 = Icges(5,5) * t95 - Icges(5,6) * t94 + Icges(5,3) * t113;
t50 = Icges(5,4) * t95 - Icges(5,2) * t94 + Icges(5,6) * t113;
t54 = Icges(5,1) * t95 - Icges(5,4) * t94 + Icges(5,5) * t113;
t18 = t109 * t46 - t92 * t50 + t93 * t54;
t66 = Icges(4,5) * t110 - Icges(4,6) * t109 + Icges(4,3) * t149;
t68 = Icges(4,4) * t110 - Icges(4,2) * t109 + Icges(4,6) * t149;
t70 = Icges(4,1) * t110 - Icges(4,4) * t109 + Icges(4,5) * t149;
t33 = -t109 * t68 + t110 * t70 + t66 * t149;
t67 = Icges(4,5) * t114 - Icges(4,6) * t113 + Icges(4,3) * t148;
t69 = Icges(4,4) * t114 - Icges(4,2) * t113 + Icges(4,6) * t148;
t71 = Icges(4,1) * t114 - Icges(4,4) * t113 + Icges(4,5) * t148;
t34 = -t109 * t69 + t110 * t71 + t67 * t149;
t166 = (-t110 * t105 + t170 * t109 - t171 * t93 - t173 * t92) * t128 + ((t16 + t18 + t34) * t124 + (t15 + t17 + t33 - t144) * t123) * t126;
t19 = t113 * t47 + t94 * t43 + t95 * t51;
t20 = t113 * t48 + t94 * t44 + t95 * t52;
t21 = t113 * t45 - t94 * t49 + t95 * t53;
t22 = t113 * t46 - t94 * t50 + t95 * t54;
t35 = -t113 * t68 + t114 * t70 + t66 * t148;
t36 = -t113 * t69 + t114 * t71 + t67 * t148;
t165 = (-t114 * t105 + t170 * t113 - t171 * t95 - t173 * t94) * t128 + ((t36 - t144 + t20 + t22) * t124 + (t35 + t19 + t21) * t123) * t126;
t164 = t128 * t66 - (-t125 * t68 + t127 * t70) * t126 + (-t45 - t47) * t147 + (-t51 - t53) * t112 + (-t43 + t49) * t111;
t163 = -t128 * t67 + (-t125 * t69 + t127 * t71) * t126 + (t46 + t48) * t147 + (t52 + t54) * t112 + (t44 - t50) * t111;
t160 = m(5) + m(6);
t156 = t109 * rSges(6,2) + t168 * t92 + t169 * t93;
t155 = t113 * rSges(6,2) + t168 * t94 + t169 * t95;
t58 = t95 * rSges(5,1) - t94 * rSges(5,2) + t113 * rSges(5,3);
t91 = t114 * pkin(3) + t113 * qJ(4);
t154 = -t58 - t91;
t153 = rSges(6,2) * t147 + t168 * t111 + t169 * t112;
t115 = (pkin(3) * t127 + qJ(4) * t125) * t126;
t89 = t110 * pkin(3) + t109 * qJ(4);
t152 = t115 * t149 + t128 * t89;
t81 = t112 * rSges(5,1) - t111 * rSges(5,2) + rSges(5,3) * t147;
t151 = -t115 - t81;
t106 = -t128 * rSges(4,3) + (rSges(4,1) * t127 - rSges(4,2) * t125) * t126;
t120 = t126 * pkin(2) - t128 * pkin(6);
t143 = -t106 - t120;
t142 = t167 * (pkin(2) * t128 + pkin(6) * t126);
t141 = -t91 - t155;
t140 = -t115 - t153;
t139 = -t120 + t151;
t137 = t123 * t89 + t124 * t91 + t142;
t136 = -t120 + t140;
t130 = Icges(3,5) * t128 - Icges(3,6) * t126;
t117 = t126 * rSges(3,1) + t128 * rSges(3,2);
t98 = Icges(3,3) * t123 + t130 * t124;
t97 = -Icges(3,3) * t124 + t130 * t123;
t88 = t143 * t124;
t87 = t143 * t123;
t82 = t89 * t148;
t79 = t114 * rSges(4,1) - t113 * rSges(4,2) + rSges(4,3) * t148;
t78 = t110 * rSges(4,1) - t109 * rSges(4,2) + rSges(4,3) * t149;
t65 = t167 * (rSges(3,1) * t128 - rSges(3,2) * t126);
t62 = -t106 * t148 - t128 * t79;
t61 = t106 * t149 + t128 * t78;
t60 = t139 * t124;
t59 = t139 * t123;
t56 = t93 * rSges(5,1) - t92 * rSges(5,2) + t109 * rSges(5,3);
t42 = (-t123 * t79 + t124 * t78) * t126;
t41 = t136 * t124;
t40 = t136 * t123;
t39 = t123 * t78 + t124 * t79 + t142;
t32 = t154 * t128 + t151 * t148;
t31 = t128 * t56 + t81 * t149 + t152;
t30 = t82 + (t154 * t123 + t124 * t56) * t126;
t25 = t123 * t56 + t124 * t58 + t137;
t24 = t141 * t128 + t140 * t148;
t23 = t156 * t128 + t153 * t149 + t152;
t14 = t36 * t123 - t35 * t124;
t13 = t34 * t123 - t33 * t124;
t12 = t82 + (t141 * t123 + t156 * t124) * t126;
t11 = t156 * t123 + t155 * t124 + t137;
t8 = t22 * t123 - t21 * t124;
t7 = t20 * t123 - t19 * t124;
t6 = t18 * t123 - t17 * t124;
t5 = t16 * t123 - t15 * t124;
t1 = [m(2) + m(3) + m(4) + t160; m(3) * t65 + m(4) * t39 + m(5) * t25 + m(6) * t11; m(6) * (t11 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(5) * (t25 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t39 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(3) * (t167 * t117 ^ 2 + t65 ^ 2) + (-t161 * t97 - t13 - t5 - t6) * t124 + (t162 * t98 + t14 + t7 + t8 + (-t123 * t97 + t124 * t98) * t124) * t123; m(4) * t42 + m(5) * t30 + m(6) * t12; m(6) * (t11 * t12 + t23 * t41 + t24 * t40) + m(5) * (t30 * t25 + t31 * t60 + t32 * t59) + m(4) * (t42 * t39 + t61 * t88 + t62 * t87) + ((t8 / 0.2e1 + t7 / 0.2e1 + t14 / 0.2e1) * t124 + (t5 / 0.2e1 + t6 / 0.2e1 + t13 / 0.2e1) * t123) * t126 + t165 * t123 / 0.2e1 - t166 * t124 / 0.2e1 - (t163 * t123 + t164 * t124) * t128 / 0.2e1; m(6) * (t12 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(4) * (t42 ^ 2 + t61 ^ 2 + t62 ^ 2) + t166 * t149 + t165 * t148 + ((t173 * t111 + t171 * t112 - t144) * t128 + ((-t104 * t125 + t105 * t127) * t128 + t172 * t146 - t163 * t124 + t164 * t123) * t126) * t128; t160 * t147; m(6) * (t109 * t40 + t11 * t147 + t113 * t41) + m(5) * (t109 * t59 + t113 * t60 + t25 * t147); m(6) * (t109 * t24 + t113 * t23 + t12 * t147) + m(5) * (t109 * t32 + t113 * t31 + t30 * t147); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t126 ^ 2 * t125 ^ 2 + t109 ^ 2 + t113 ^ 2); m(6) * t111; m(6) * (t111 * t11 + t92 * t40 + t94 * t41); m(6) * (t111 * t12 + t94 * t23 + t92 * t24); m(6) * (t92 * t109 + t111 * t147 + t94 * t113); m(6) * (t111 ^ 2 + t92 ^ 2 + t94 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
