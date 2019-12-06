% Calculate joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:30
% EndTime: 2019-12-05 16:21:38
% DurationCPUTime: 2.50s
% Computational Cost: add. (5874->308), mult. (7162->485), div. (0->0), fcn. (7771->10), ass. (0->159)
t147 = qJ(3) + pkin(9);
t140 = sin(t147);
t141 = cos(t147);
t153 = sin(qJ(2));
t155 = cos(qJ(2));
t103 = -Icges(5,3) * t155 + (Icges(5,5) * t141 - Icges(5,6) * t140) * t153;
t152 = sin(qJ(3));
t154 = cos(qJ(3));
t121 = -Icges(4,3) * t155 + (Icges(4,5) * t154 - Icges(4,6) * t152) * t153;
t203 = (-t121 - t103) * t155;
t149 = sin(pkin(8));
t145 = t149 ^ 2;
t150 = cos(pkin(8));
t146 = t150 ^ 2;
t173 = t145 + t146;
t104 = -Icges(5,6) * t155 + (Icges(5,4) * t141 - Icges(5,2) * t140) * t153;
t105 = -Icges(5,5) * t155 + (Icges(5,1) * t141 - Icges(5,4) * t140) * t153;
t186 = t149 * t155;
t117 = -t140 * t186 - t150 * t141;
t118 = -t150 * t140 + t141 * t186;
t122 = -Icges(4,6) * t155 + (Icges(4,4) * t154 - Icges(4,2) * t152) * t153;
t123 = -Icges(4,5) * t155 + (Icges(4,1) * t154 - Icges(4,4) * t152) * t153;
t182 = t152 * t155;
t128 = -t149 * t182 - t150 * t154;
t181 = t154 * t155;
t185 = t150 * t152;
t129 = t149 * t181 - t185;
t187 = t149 * t153;
t70 = Icges(5,5) * t118 + Icges(5,6) * t117 + Icges(5,3) * t187;
t72 = Icges(5,4) * t118 + Icges(5,2) * t117 + Icges(5,6) * t187;
t74 = Icges(5,1) * t118 + Icges(5,4) * t117 + Icges(5,5) * t187;
t31 = t117 * t72 + t118 * t74 + t70 * t187;
t183 = t150 * t155;
t119 = -t140 * t183 + t149 * t141;
t120 = t149 * t140 + t141 * t183;
t184 = t150 * t153;
t71 = Icges(5,5) * t120 + Icges(5,6) * t119 + Icges(5,3) * t184;
t73 = Icges(5,4) * t120 + Icges(5,2) * t119 + Icges(5,6) * t184;
t75 = Icges(5,1) * t120 + Icges(5,4) * t119 + Icges(5,5) * t184;
t32 = t117 * t73 + t118 * t75 + t71 * t187;
t83 = Icges(4,5) * t129 + Icges(4,6) * t128 + Icges(4,3) * t187;
t85 = Icges(4,4) * t129 + Icges(4,2) * t128 + Icges(4,6) * t187;
t87 = Icges(4,1) * t129 + Icges(4,4) * t128 + Icges(4,5) * t187;
t39 = t128 * t85 + t129 * t87 + t83 * t187;
t130 = t149 * t154 - t150 * t182;
t188 = t149 * t152;
t131 = t150 * t181 + t188;
t84 = Icges(4,5) * t131 + Icges(4,6) * t130 + Icges(4,3) * t184;
t86 = Icges(4,4) * t131 + Icges(4,2) * t130 + Icges(4,6) * t184;
t88 = Icges(4,1) * t131 + Icges(4,4) * t130 + Icges(4,5) * t184;
t40 = t128 * t86 + t129 * t88 + t84 * t187;
t202 = (-t117 * t104 - t118 * t105 - t128 * t122 - t129 * t123) * t155 + ((t40 + t32) * t150 + (t39 + t31 + t203) * t149) * t153;
t33 = t119 * t72 + t120 * t74 + t70 * t184;
t34 = t119 * t73 + t120 * t75 + t71 * t184;
t41 = t130 * t85 + t131 * t87 + t83 * t184;
t42 = t130 * t86 + t131 * t88 + t84 * t184;
t201 = (-t119 * t104 - t120 * t105 - t130 * t122 - t131 * t123) * t155 + ((t42 + t34 + t203) * t150 + (t41 + t33) * t149) * t153;
t200 = t155 ^ 2;
t199 = -m(5) - m(6);
t198 = t149 / 0.2e1;
t197 = -t150 / 0.2e1;
t196 = -t155 / 0.2e1;
t194 = t154 * pkin(3);
t175 = pkin(4) * t141;
t156 = pkin(7) * t153 + t175 * t155;
t169 = pkin(4) * t140;
t142 = qJ(5) + t147;
t137 = sin(t142);
t138 = cos(t142);
t109 = -t137 * t183 + t149 * t138;
t110 = t149 * t137 + t138 * t183;
t69 = t110 * rSges(6,1) + t109 * rSges(6,2) + rSges(6,3) * t184;
t192 = t169 * t149 + t156 * t150 + t69;
t77 = t120 * rSges(5,1) + t119 * rSges(5,2) + rSges(5,3) * t184;
t157 = qJ(4) * t153 + t194 * t155;
t90 = pkin(3) * t188 + t157 * t150;
t191 = -t77 - t90;
t102 = -qJ(4) * t155 + t194 * t153;
t89 = -pkin(3) * t185 + t157 * t149;
t190 = t102 * t187 + t155 * t89;
t101 = -t155 * rSges(6,3) + (rSges(6,1) * t138 - rSges(6,2) * t137) * t153;
t107 = -t137 * t186 - t150 * t138;
t108 = -t150 * t137 + t138 * t186;
t68 = t108 * rSges(6,1) + t107 * rSges(6,2) + rSges(6,3) * t187;
t52 = t101 * t187 + t155 * t68;
t98 = -Icges(6,3) * t155 + (Icges(6,5) * t138 - Icges(6,6) * t137) * t153;
t189 = t155 * t98;
t106 = -t155 * rSges(5,3) + (rSges(5,1) * t141 - rSges(5,2) * t140) * t153;
t178 = -t102 - t106;
t124 = -t155 * rSges(4,3) + (rSges(4,1) * t154 - rSges(4,2) * t152) * t153;
t136 = t153 * pkin(2) - t155 * pkin(6);
t177 = -t124 - t136;
t176 = t173 * (pkin(2) * t155 + pkin(6) * t153);
t172 = -t90 - t192;
t93 = -pkin(7) * t155 + t175 * t153;
t171 = -t101 - t102 - t93;
t170 = -t136 + t178;
t168 = t149 * t89 + t150 * t90 + t176;
t100 = -Icges(6,5) * t155 + (Icges(6,1) * t138 - Icges(6,4) * t137) * t153;
t62 = Icges(6,5) * t108 + Icges(6,6) * t107 + Icges(6,3) * t187;
t64 = Icges(6,4) * t108 + Icges(6,2) * t107 + Icges(6,6) * t187;
t66 = Icges(6,1) * t108 + Icges(6,4) * t107 + Icges(6,5) * t187;
t35 = -t155 * t62 + (-t137 * t64 + t138 * t66) * t153;
t63 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t184;
t65 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t184;
t67 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t184;
t36 = -t155 * t63 + (-t137 * t65 + t138 * t67) * t153;
t26 = t107 * t64 + t108 * t66 + t62 * t187;
t27 = t107 * t65 + t108 * t67 + t63 * t187;
t99 = -Icges(6,6) * t155 + (Icges(6,4) * t138 - Icges(6,2) * t137) * t153;
t5 = -(t108 * t100 + t107 * t99) * t155 + (t27 * t150 + (t26 - t189) * t149) * t153;
t28 = t109 * t64 + t110 * t66 + t62 * t184;
t29 = t109 * t65 + t110 * t67 + t63 * t184;
t6 = -(t110 * t100 + t109 * t99) * t155 + (t28 * t149 + (t29 - t189) * t150) * t153;
t167 = -t155 * (t200 * t98 + (t36 * t150 + t35 * t149 - (t100 * t138 - t137 * t99) * t155) * t153) + t6 * t184 + t5 * t187;
t166 = -t136 + t171;
t12 = t27 * t149 - t26 * t150;
t13 = t29 * t149 - t28 * t150;
t165 = t12 * t187 / 0.2e1 + t5 * t197 + t6 * t198 + t13 * t184 / 0.2e1 + (t36 * t149 - t35 * t150) * t196;
t160 = Icges(3,5) * t155 - Icges(3,6) * t153;
t135 = t153 * rSges(3,1) + t155 * rSges(3,2);
t112 = Icges(3,3) * t149 + t160 * t150;
t111 = -Icges(3,3) * t150 + t160 * t149;
t95 = t177 * t150;
t94 = t177 * t149;
t92 = t131 * rSges(4,1) + t130 * rSges(4,2) + rSges(4,3) * t184;
t91 = t129 * rSges(4,1) + t128 * rSges(4,2) + rSges(4,3) * t187;
t79 = t173 * (rSges(3,1) * t155 - rSges(3,2) * t153);
t78 = t89 * t184;
t76 = t118 * rSges(5,1) + t117 * rSges(5,2) + rSges(5,3) * t187;
t60 = t68 * t184;
t59 = t170 * t150;
t58 = t170 * t149;
t56 = t156 * t149 - t169 * t150;
t55 = -t124 * t184 - t155 * t92;
t54 = t124 * t187 + t155 * t91;
t53 = -t101 * t184 - t155 * t69;
t51 = (-t149 * t92 + t150 * t91) * t153;
t50 = t166 * t150;
t49 = t166 * t149;
t48 = -t69 * t187 + t60;
t47 = t149 * t91 + t150 * t92 + t176;
t46 = -t155 * t84 + (-t152 * t86 + t154 * t88) * t153;
t45 = -t155 * t83 + (-t152 * t85 + t154 * t87) * t153;
t44 = t191 * t155 + t178 * t184;
t43 = t106 * t187 + t155 * t76 + t190;
t38 = -t155 * t71 + (-t140 * t73 + t141 * t75) * t153;
t37 = -t155 * t70 + (-t140 * t72 + t141 * t74) * t153;
t30 = t78 + (t191 * t149 + t150 * t76) * t153;
t25 = t149 * t76 + t150 * t77 + t168;
t24 = t172 * t155 + t171 * t184;
t23 = t155 * t56 + t93 * t187 + t190 + t52;
t22 = t42 * t149 - t41 * t150;
t21 = t40 * t149 - t39 * t150;
t20 = t60 + t78 + (t172 * t149 + t150 * t56) * t153;
t18 = t192 * t150 + (t56 + t68) * t149 + t168;
t17 = t34 * t149 - t33 * t150;
t16 = t32 * t149 - t31 * t150;
t1 = [m(2) + m(3) + m(4) - t199; m(3) * t79 + m(4) * t47 + m(5) * t25 + m(6) * t18; m(6) * (t18 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t25 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t47 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(3) * (t173 * t135 ^ 2 + t79 ^ 2) + (-t146 * t111 - t12 - t16 - t21) * t150 + (t145 * t112 + t13 + t17 + t22 + (-t149 * t111 + t150 * t112) * t150) * t149; m(4) * t51 + m(5) * t30 + m(6) * t20; m(6) * (t20 * t18 + t23 * t50 + t24 * t49) + m(5) * (t30 * t25 + t43 * t59 + t44 * t58) + m(4) * (t51 * t47 + t54 * t95 + t55 * t94) + ((t17 / 0.2e1 + t22 / 0.2e1) * t150 + (t21 / 0.2e1 + t16 / 0.2e1) * t149) * t153 + t165 + t201 * t198 + t202 * t197 + ((-t37 - t45) * t150 + (t38 + t46) * t149) * t196; m(6) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(5) * (t30 ^ 2 + t43 ^ 2 + t44 ^ 2) - t155 * (t200 * t121 + (t46 * t150 + t45 * t149 - (-t122 * t152 + t123 * t154) * t155) * t153) - t155 * (t200 * t103 + (t38 * t150 + t37 * t149 - (-t104 * t140 + t105 * t141) * t155) * t153) + m(4) * (t51 ^ 2 + t54 ^ 2 + t55 ^ 2) + t167 + t202 * t187 + t201 * t184; t199 * t155; m(6) * (-t155 * t18 + (t149 * t49 + t150 * t50) * t153) + m(5) * (-t155 * t25 + (t149 * t58 + t150 * t59) * t153); m(6) * (-t155 * t20 + (t149 * t24 + t150 * t23) * t153) + m(5) * (-t155 * t30 + (t149 * t44 + t150 * t43) * t153); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t173 * t153 ^ 2 + t200); m(6) * t48; m(6) * (t48 * t18 + t53 * t49 + t52 * t50) + t165; m(6) * (t48 * t20 + t52 * t23 + t53 * t24) + t167; m(6) * (-t48 * t155 + (t149 * t53 + t150 * t52) * t153); m(6) * (t48 ^ 2 + t52 ^ 2 + t53 ^ 2) + t167;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
