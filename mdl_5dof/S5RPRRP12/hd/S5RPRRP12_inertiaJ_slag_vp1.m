% Calculate joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:17
% EndTime: 2019-12-31 18:56:21
% DurationCPUTime: 1.67s
% Computational Cost: add. (1974->277), mult. (4500->401), div. (0->0), fcn. (4769->6), ass. (0->135)
t126 = cos(qJ(3));
t187 = Icges(4,5) * t126;
t123 = sin(qJ(3));
t186 = Icges(4,6) * t123;
t121 = -qJ(5) - pkin(7);
t185 = -rSges(6,3) + t121;
t184 = t187 / 0.2e1 - t186 / 0.2e1;
t125 = cos(qJ(4));
t122 = sin(qJ(4));
t69 = Icges(6,3) * t123 + (Icges(6,5) * t125 - Icges(6,6) * t122) * t126;
t70 = Icges(5,3) * t123 + (Icges(5,5) * t125 - Icges(5,6) * t122) * t126;
t77 = Icges(6,5) * t123 + (Icges(6,1) * t125 - Icges(6,4) * t122) * t126;
t78 = Icges(5,5) * t123 + (Icges(5,1) * t125 - Icges(5,4) * t122) * t126;
t183 = (t77 + t78) * t125 * t126 + (t69 + t70) * t123;
t73 = Icges(6,6) * t123 + (Icges(6,4) * t125 - Icges(6,2) * t122) * t126;
t74 = Icges(5,6) * t123 + (Icges(5,4) * t125 - Icges(5,2) * t122) * t126;
t182 = (-t73 - t74) * t122;
t124 = sin(qJ(1));
t152 = t124 * t126;
t127 = cos(qJ(1));
t148 = t127 * t125;
t154 = t124 * t122;
t90 = -t123 * t154 + t148;
t149 = t127 * t122;
t153 = t124 * t125;
t91 = t123 * t153 + t149;
t42 = Icges(6,5) * t91 + Icges(6,6) * t90 - Icges(6,3) * t152;
t46 = Icges(6,4) * t91 + Icges(6,2) * t90 - Icges(6,6) * t152;
t50 = Icges(6,1) * t91 + Icges(6,4) * t90 - Icges(6,5) * t152;
t11 = -t42 * t152 + t90 * t46 + t91 * t50;
t150 = t126 * t127;
t92 = t123 * t149 + t153;
t93 = -t123 * t148 + t154;
t43 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t150;
t47 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t150;
t51 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t150;
t12 = -t43 * t152 + t90 * t47 + t91 * t51;
t44 = Icges(5,5) * t91 + Icges(5,6) * t90 - Icges(5,3) * t152;
t48 = Icges(5,4) * t91 + Icges(5,2) * t90 - Icges(5,6) * t152;
t52 = Icges(5,1) * t91 + Icges(5,4) * t90 - Icges(5,5) * t152;
t13 = -t44 * t152 + t90 * t48 + t91 * t52;
t45 = Icges(5,5) * t93 + Icges(5,6) * t92 + Icges(5,3) * t150;
t49 = Icges(5,4) * t93 + Icges(5,2) * t92 + Icges(5,6) * t150;
t53 = Icges(5,1) * t93 + Icges(5,4) * t92 + Icges(5,5) * t150;
t14 = -t45 * t152 + t90 * t49 + t91 * t53;
t26 = -t69 * t152 + t90 * t73 + t91 * t77;
t27 = -t70 * t152 + t90 * t74 + t91 * t78;
t181 = ((t12 + t14) * t127 + (-t11 - t13) * t124) * t126 + (t26 + t27) * t123;
t15 = t42 * t150 + t92 * t46 + t93 * t50;
t16 = t43 * t150 + t92 * t47 + t93 * t51;
t17 = t44 * t150 + t92 * t48 + t93 * t52;
t18 = t45 * t150 + t92 * t49 + t93 * t53;
t28 = t69 * t150 + t92 * t73 + t93 * t77;
t29 = t70 * t150 + t92 * t74 + t93 * t78;
t180 = ((t16 + t18) * t127 + (-t15 - t17) * t124) * t126 + (t28 + t29) * t123;
t21 = t123 * t42 + (-t122 * t46 + t125 * t50) * t126;
t23 = t123 * t44 + (-t122 * t48 + t125 * t52) * t126;
t179 = t21 + t23;
t22 = t123 * t43 + (-t122 * t47 + t125 * t51) * t126;
t24 = t123 * t45 + (-t122 * t49 + t125 * t53) * t126;
t178 = t22 + t24;
t112 = t125 * pkin(4) + pkin(3);
t155 = t123 * t124;
t177 = t91 * rSges(6,1) + t90 * rSges(6,2) + pkin(4) * t149 + t112 * t155 + t185 * t152;
t176 = (rSges(4,1) * t123 + rSges(4,2) * t126) * t127;
t118 = t124 ^ 2;
t120 = t127 ^ 2;
t175 = -pkin(1) - pkin(6);
t172 = t124 / 0.2e1;
t170 = t127 / 0.2e1;
t102 = t126 * rSges(4,1) - t123 * rSges(4,2);
t169 = m(4) * t102;
t168 = m(6) * t126;
t167 = -pkin(7) - t121;
t166 = (t126 * t182 + t183) * t123;
t109 = pkin(3) * t155;
t94 = -pkin(7) * t152 + t109;
t165 = -t94 + t177;
t111 = t127 * t123 * pkin(3);
t139 = -t93 * rSges(6,1) - t92 * rSges(6,2);
t157 = t112 * t123;
t164 = pkin(4) * t154 + t111 + (t167 * t126 - t157) * t127 + rSges(6,3) * t150 - t139;
t161 = (rSges(6,1) * t125 - rSges(6,2) * t122 - pkin(3) + t112) * t126 + (t167 + rSges(6,3)) * t123;
t160 = t91 * rSges(5,1) + t90 * rSges(5,2);
t147 = t127 * pkin(1) + t124 * qJ(2);
t146 = t118 + t120;
t144 = rSges(4,1) * t155 + rSges(4,2) * t152 + t127 * rSges(4,3);
t143 = t127 * pkin(6) + t147;
t142 = (-rSges(5,3) - pkin(7)) * t126;
t141 = t161 * t126;
t140 = -t93 * rSges(5,1) - t92 * rSges(5,2);
t19 = t165 * t123 + t124 * t141;
t20 = -t164 * t123 + t127 * t141;
t135 = t20 * t124 - t19 * t127;
t114 = t127 * qJ(2);
t31 = t114 + (t185 * t126 + t157) * t127 + (-pkin(4) * t122 + t175) * t124 + t139;
t32 = t143 + t177;
t134 = t124 * t31 - t127 * t32;
t104 = t126 * pkin(3) + t123 * pkin(7);
t95 = t124 * t104;
t39 = t161 * t124 + t95;
t40 = (-t104 - t161) * t127;
t133 = t39 * t124 - t40 * t127;
t130 = Icges(4,5) * t123 + Icges(4,6) * t126;
t129 = -t27 / 0.2e1 - t26 / 0.2e1 - t23 / 0.2e1 - t21 / 0.2e1;
t128 = t29 / 0.2e1 + t28 / 0.2e1 + t24 / 0.2e1 + t22 / 0.2e1;
t103 = t127 * rSges(2,1) - t124 * rSges(2,2);
t101 = -t124 * rSges(2,1) - t127 * rSges(2,2);
t98 = -t186 + t187;
t85 = t127 * (pkin(7) * t150 - t111);
t84 = -t127 * rSges(3,2) + t124 * rSges(3,3) + t147;
t83 = t127 * rSges(3,3) + t114 + (rSges(3,2) - pkin(1)) * t124;
t82 = t123 * rSges(5,3) + (rSges(5,1) * t125 - rSges(5,2) * t122) * t126;
t72 = Icges(4,3) * t124 - t130 * t127;
t71 = Icges(4,3) * t127 + t130 * t124;
t63 = t143 + t144;
t62 = t114 + t176 + (-rSges(4,3) + t175) * t124;
t61 = (-t104 - t82) * t127;
t60 = t124 * t82 + t95;
t59 = rSges(5,3) * t150 - t140;
t57 = -rSges(5,3) * t152 + t160;
t41 = -t124 * t144 + (t124 * rSges(4,3) - t176) * t127;
t38 = t124 * t142 + t109 + t143 + t160;
t37 = t175 * t124 + t127 * t142 + t111 + t114 + t140;
t36 = -t123 * t59 + t82 * t150;
t35 = t123 * t57 + t82 * t152;
t30 = (-t124 * t59 - t127 * t57) * t126;
t25 = t127 * t59 + t85 + (-t57 - t94) * t124;
t10 = (-t164 * t124 - t165 * t127) * t126;
t9 = t85 + t164 * t127 + (-t94 - t165) * t124;
t8 = t18 * t124 + t17 * t127;
t7 = t16 * t124 + t15 * t127;
t6 = t14 * t124 + t13 * t127;
t5 = t11 * t127 + t12 * t124;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t126 + t182) * t126 + m(6) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t83 ^ 2 + t84 ^ 2) + m(2) * (t101 ^ 2 + t103 ^ 2) + t183 + (-0.2e1 * Icges(4,4) * t126 + Icges(4,2) * t123) * t123; m(6) * t134 + m(5) * (t124 * t37 - t127 * t38) + m(4) * (t124 * t62 - t127 * t63) + m(3) * (t124 * t83 - t127 * t84); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t146; m(6) * (t39 * t31 + t40 * t32) + m(5) * (t60 * t37 + t61 * t38) + (t184 * t127 - t63 * t169 + t98 * t170 - t129) * t127 + (t184 * t124 + t62 * t169 + t98 * t172 + t128) * t124; m(5) * (t60 * t124 - t61 * t127) + m(6) * t133 + t146 * t169; m(6) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(5) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t146 * t102 ^ 2 + t41 ^ 2) + (t118 * t72 + t7 + t8) * t124 + (t120 * t71 + t5 + t6 + (t124 * t71 + t127 * t72) * t124) * t127; m(6) * (t19 * t32 + t20 * t31) + m(5) * (t35 * t38 + t36 * t37) + (t129 * t124 + t128 * t127) * t126 + t166; m(5) * (t36 * t124 - t35 * t127) + m(6) * t135; m(6) * (t10 * t9 + t19 * t40 + t20 * t39) + m(5) * (t30 * t25 + t35 * t61 + t36 * t60) + ((t8 / 0.2e1 + t7 / 0.2e1) * t127 + (-t5 / 0.2e1 - t6 / 0.2e1) * t124) * t126 + (t178 * t124 + t179 * t127) * t123 / 0.2e1 + t180 * t172 + t181 * t170; t166 * t123 + m(6) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(5) * (t30 ^ 2 + t35 ^ 2 + t36 ^ 2) + (t180 * t127 - t181 * t124 + (-t179 * t124 + t178 * t127) * t123) * t126; -t134 * t168; -t146 * t168; m(6) * (t123 * t9 - t133 * t126); m(6) * (t123 * t10 - t135 * t126); m(6) * (t146 * t126 ^ 2 + t123 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
