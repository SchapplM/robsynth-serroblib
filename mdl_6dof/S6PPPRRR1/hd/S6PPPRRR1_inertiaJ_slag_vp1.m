% Calculate joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPPRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:15
% EndTime: 2019-03-08 18:39:26
% DurationCPUTime: 5.71s
% Computational Cost: add. (51403->358), mult. (146347->557), div. (0->0), fcn. (195576->18), ass. (0->172)
t159 = m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t187 = 0.2e1 * t159;
t137 = sin(pkin(12));
t139 = cos(pkin(12));
t140 = cos(pkin(6));
t173 = cos(pkin(13));
t166 = t140 * t173;
t169 = sin(pkin(13));
t153 = -t137 * t166 - t139 * t169;
t175 = cos(pkin(7));
t151 = t153 * t175;
t165 = t140 * t169;
t154 = -t137 * t165 + t139 * t173;
t138 = sin(pkin(6));
t168 = sin(pkin(14));
t171 = sin(pkin(7));
t160 = t171 * t168;
t157 = t138 * t160;
t172 = cos(pkin(14));
t125 = t137 * t157 + t151 * t168 + t154 * t172;
t167 = t138 * t175;
t132 = t137 * t167 - t153 * t171;
t143 = sin(qJ(4));
t161 = t172 * t171;
t158 = t138 * t161;
t147 = -t137 * t158 - t151 * t172 + t154 * t168;
t174 = cos(pkin(8));
t145 = t147 * t174;
t170 = sin(pkin(8));
t180 = cos(qJ(4));
t112 = t125 * t180 + (t132 * t170 - t145) * t143;
t119 = t132 * t174 + t147 * t170;
t142 = sin(qJ(5));
t179 = cos(qJ(5));
t100 = t112 * t142 - t119 * t179;
t162 = t175 * t173;
t130 = t140 * t160 + (t162 * t168 + t169 * t172) * t138;
t135 = -t138 * t171 * t173 + t140 * t175;
t150 = t140 * t161 + (t162 * t172 - t168 * t169) * t138;
t149 = t150 * t174;
t116 = t130 * t180 + (t135 * t170 + t149) * t143;
t126 = t135 * t174 - t150 * t170;
t107 = t116 * t142 - t126 * t179;
t155 = -t137 * t169 + t139 * t166;
t152 = t155 * t175;
t156 = t137 * t173 + t139 * t165;
t124 = -t139 * t157 + t152 * t168 + t156 * t172;
t131 = -t139 * t167 - t155 * t171;
t148 = t139 * t158 - t152 * t172 + t156 * t168;
t146 = t148 * t174;
t163 = t180 * t170;
t109 = t124 * t143 - t131 * t163 + t146 * t180;
t141 = sin(qJ(6));
t144 = cos(qJ(6));
t110 = t124 * t180 + (t131 * t170 - t146) * t143;
t118 = t131 * t174 + t148 * t170;
t99 = t110 * t179 + t118 * t142;
t74 = t109 * t144 - t99 * t141;
t75 = t109 * t141 + t99 * t144;
t98 = t110 * t142 - t118 * t179;
t49 = Icges(7,5) * t75 + Icges(7,6) * t74 + Icges(7,3) * t98;
t51 = Icges(7,4) * t75 + Icges(7,2) * t74 + Icges(7,6) * t98;
t53 = Icges(7,1) * t75 + Icges(7,4) * t74 + Icges(7,5) * t98;
t17 = t98 * t49 + t74 * t51 + t75 * t53;
t101 = t112 * t179 + t119 * t142;
t111 = t125 * t143 - t132 * t163 + t145 * t180;
t76 = -t101 * t141 + t111 * t144;
t77 = t101 * t144 + t111 * t141;
t50 = Icges(7,5) * t77 + Icges(7,6) * t76 + Icges(7,3) * t100;
t52 = Icges(7,4) * t77 + Icges(7,2) * t76 + Icges(7,6) * t100;
t54 = Icges(7,1) * t77 + Icges(7,4) * t76 + Icges(7,5) * t100;
t18 = t98 * t50 + t74 * t52 + t75 * t54;
t108 = t116 * t179 + t126 * t142;
t115 = t130 * t143 - t135 * t163 - t149 * t180;
t96 = -t108 * t141 + t115 * t144;
t97 = t108 * t144 + t115 * t141;
t60 = Icges(7,5) * t97 + Icges(7,6) * t96 + Icges(7,3) * t107;
t61 = Icges(7,4) * t97 + Icges(7,2) * t96 + Icges(7,6) * t107;
t62 = Icges(7,1) * t97 + Icges(7,4) * t96 + Icges(7,5) * t107;
t26 = t98 * t60 + t74 * t61 + t75 * t62;
t1 = t18 * t100 + t26 * t107 + t17 * t98;
t186 = t1 / 0.2e1;
t19 = t100 * t49 + t76 * t51 + t77 * t53;
t20 = t100 * t50 + t76 * t52 + t77 * t54;
t27 = t100 * t60 + t76 * t61 + t77 * t62;
t2 = t20 * t100 + t27 * t107 + t19 * t98;
t185 = t2 / 0.2e1;
t21 = t107 * t49 + t96 * t51 + t97 * t53;
t22 = t107 * t50 + t96 * t52 + t97 * t54;
t30 = t107 * t60 + t96 * t61 + t97 * t62;
t7 = t22 * t100 + t30 * t107 + t21 * t98;
t184 = t7 / 0.2e1;
t183 = t98 / 0.2e1;
t182 = t100 / 0.2e1;
t181 = t107 / 0.2e1;
t55 = t75 * rSges(7,1) + t74 * rSges(7,2) + t98 * rSges(7,3);
t178 = t99 * pkin(5) + t98 * pkin(11) + t55;
t56 = t77 * rSges(7,1) + t76 * rSges(7,2) + t100 * rSges(7,3);
t177 = t101 * pkin(5) + t100 * pkin(11) + t56;
t63 = t97 * rSges(7,1) + t96 * rSges(7,2) + t107 * rSges(7,3);
t176 = t108 * pkin(5) + t107 * pkin(11) + t63;
t164 = m(3) + m(4) + m(5) + m(6) + m(7);
t106 = t116 * pkin(4) + t115 * pkin(10);
t105 = t116 * rSges(5,1) - t115 * rSges(5,2) + t126 * rSges(5,3);
t104 = Icges(5,1) * t116 - Icges(5,4) * t115 + Icges(5,5) * t126;
t103 = Icges(5,4) * t116 - Icges(5,2) * t115 + Icges(5,6) * t126;
t102 = Icges(5,5) * t116 - Icges(5,6) * t115 + Icges(5,3) * t126;
t95 = t118 * t106;
t94 = t112 * pkin(4) + t111 * pkin(10);
t93 = t110 * pkin(4) + t109 * pkin(10);
t91 = t126 * t94;
t90 = t119 * t93;
t89 = t112 * rSges(5,1) - t111 * rSges(5,2) + t119 * rSges(5,3);
t88 = t110 * rSges(5,1) - t109 * rSges(5,2) + t118 * rSges(5,3);
t87 = Icges(5,1) * t112 - Icges(5,4) * t111 + Icges(5,5) * t119;
t86 = Icges(5,1) * t110 - Icges(5,4) * t109 + Icges(5,5) * t118;
t85 = Icges(5,4) * t112 - Icges(5,2) * t111 + Icges(5,6) * t119;
t84 = Icges(5,4) * t110 - Icges(5,2) * t109 + Icges(5,6) * t118;
t83 = Icges(5,5) * t112 - Icges(5,6) * t111 + Icges(5,3) * t119;
t82 = Icges(5,5) * t110 - Icges(5,6) * t109 + Icges(5,3) * t118;
t81 = t108 * rSges(6,1) - t107 * rSges(6,2) + t115 * rSges(6,3);
t80 = Icges(6,1) * t108 - Icges(6,4) * t107 + Icges(6,5) * t115;
t79 = Icges(6,4) * t108 - Icges(6,2) * t107 + Icges(6,6) * t115;
t78 = Icges(6,5) * t108 - Icges(6,6) * t107 + Icges(6,3) * t115;
t71 = t101 * rSges(6,1) - t100 * rSges(6,2) + t111 * rSges(6,3);
t70 = t99 * rSges(6,1) - t98 * rSges(6,2) + t109 * rSges(6,3);
t69 = Icges(6,1) * t101 - Icges(6,4) * t100 + Icges(6,5) * t111;
t68 = Icges(6,1) * t99 - Icges(6,4) * t98 + Icges(6,5) * t109;
t67 = Icges(6,4) * t101 - Icges(6,2) * t100 + Icges(6,6) * t111;
t66 = Icges(6,4) * t99 - Icges(6,2) * t98 + Icges(6,6) * t109;
t65 = Icges(6,5) * t101 - Icges(6,6) * t100 + Icges(6,3) * t111;
t64 = Icges(6,5) * t99 - Icges(6,6) * t98 + Icges(6,3) * t109;
t59 = -t119 * t105 + t126 * t89;
t58 = t118 * t105 - t126 * t88;
t57 = -t118 * t89 + t119 * t88;
t48 = -t111 * t81 + t115 * t71;
t47 = t109 * t81 - t115 * t70;
t46 = -t107 * t79 + t108 * t80 + t115 * t78;
t45 = -t109 * t71 + t111 * t70;
t44 = t126 * t71 + t91 + (-t106 - t81) * t119;
t43 = t118 * t81 + t95 + (-t70 - t93) * t126;
t42 = -t100 * t79 + t101 * t80 + t111 * t78;
t41 = t109 * t78 - t98 * t79 + t99 * t80;
t40 = -t100 * t63 + t107 * t56;
t39 = -t107 * t55 + t98 * t63;
t38 = t119 * t70 + t90 + (-t71 - t94) * t118;
t37 = -t107 * t67 + t108 * t69 + t115 * t65;
t36 = -t107 * t66 + t108 * t68 + t115 * t64;
t35 = -t100 * t67 + t101 * t69 + t111 * t65;
t34 = -t100 * t66 + t101 * t68 + t111 * t64;
t33 = t109 * t65 - t98 * t67 + t99 * t69;
t32 = t109 * t64 - t98 * t66 + t99 * t68;
t31 = t100 * t55 - t98 * t56;
t29 = -t111 * t176 + t115 * t177;
t28 = t109 * t176 - t115 * t178;
t25 = t91 + t177 * t126 + (-t106 - t176) * t119;
t24 = t95 + t176 * t118 + (-t93 - t178) * t126;
t23 = -t109 * t177 + t111 * t178;
t16 = t90 + t178 * t119 + (-t94 - t177) * t118;
t15 = t36 * t118 + t37 * t119 + t46 * t126;
t14 = t36 * t109 + t37 * t111 + t46 * t115;
t13 = t34 * t118 + t35 * t119 + t42 * t126;
t12 = t32 * t118 + t33 * t119 + t41 * t126;
t11 = t34 * t109 + t35 * t111 + t42 * t115;
t10 = t32 * t109 + t33 * t111 + t41 * t115;
t9 = t21 * t118 + t22 * t119 + t30 * t126;
t8 = t21 * t109 + t22 * t111 + t30 * t115;
t6 = t19 * t118 + t20 * t119 + t27 * t126;
t5 = t17 * t118 + t18 * t119 + t26 * t126;
t4 = t19 * t109 + t20 * t111 + t27 * t115;
t3 = t17 * t109 + t18 * t111 + t26 * t115;
t72 = [m(2) + t164; t164 * t140; 0.2e1 * (m(3) / 0.2e1 + t159) * (t140 ^ 2 + (t137 ^ 2 + t139 ^ 2) * t138 ^ 2); t135 * t187; (t135 * t140 + (-t131 * t139 + t132 * t137) * t138) * t187; (t131 ^ 2 + t132 ^ 2 + t135 ^ 2) * t187; m(5) * t57 + m(6) * t38 + m(7) * t16; m(5) * (t57 * t140 + (t137 * t58 - t139 * t59) * t138) + m(6) * (t38 * t140 + (t137 * t43 - t139 * t44) * t138) + m(7) * (t16 * t140 + (t137 * t24 - t139 * t25) * t138); m(5) * (t59 * t131 + t58 * t132 + t57 * t135) + m(6) * (t44 * t131 + t43 * t132 + t38 * t135) + m(7) * (t25 * t131 + t24 * t132 + t16 * t135); (t9 + t15 + (t126 * t102 - t115 * t103 + t116 * t104) * t126) * t126 + (t6 + t13 + (-t111 * t85 + t112 * t87 + t119 * t83) * t119 + (t119 * t102 - t111 * t103 + t112 * t104 - t115 * t85 + t116 * t87 + t126 * t83) * t126) * t119 + (t5 + t12 + (-t109 * t84 + t110 * t86 + t118 * t82) * t118 + (t118 * t102 - t109 * t103 + t110 * t104 - t115 * t84 + t116 * t86 + t126 * t82) * t126 + (-t109 * t85 + t110 * t87 - t111 * t84 + t112 * t86 + t118 * t83 + t119 * t82) * t119) * t118 + m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2); m(6) * t45 + m(7) * t23; m(6) * (t45 * t140 + (t137 * t47 - t139 * t48) * t138) + m(7) * (t23 * t140 + (t137 * t28 - t139 * t29) * t138); m(6) * (t48 * t131 + t47 * t132 + t45 * t135) + m(7) * (t29 * t131 + t28 * t132 + t23 * t135); (t8 / 0.2e1 + t14 / 0.2e1) * t126 + (t4 / 0.2e1 + t11 / 0.2e1) * t119 + (t3 / 0.2e1 + t10 / 0.2e1) * t118 + (t9 / 0.2e1 + t15 / 0.2e1) * t115 + (t6 / 0.2e1 + t13 / 0.2e1) * t111 + (t5 / 0.2e1 + t12 / 0.2e1) * t109 + m(7) * (t16 * t23 + t24 * t28 + t25 * t29) + m(6) * (t45 * t38 + t47 * t43 + t48 * t44); (t8 + t14) * t115 + (t4 + t11) * t111 + (t3 + t10) * t109 + m(7) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t45 ^ 2 + t47 ^ 2 + t48 ^ 2); m(7) * t31; m(7) * (t31 * t140 + (t137 * t39 - t139 * t40) * t138); m(7) * (t40 * t131 + t39 * t132 + t31 * t135); m(7) * (t31 * t16 + t39 * t24 + t40 * t25) + t118 * t186 + t126 * t184 + t119 * t185 + t9 * t181 + t5 * t183 + t6 * t182; t4 * t182 + t111 * t185 + t3 * t183 + m(7) * (t31 * t23 + t39 * t28 + t40 * t29) + t8 * t181 + t109 * t186 + t115 * t184; t98 * t1 + t107 * t7 + t100 * t2 + m(7) * (t31 ^ 2 + t39 ^ 2 + t40 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t72(1) t72(2) t72(4) t72(7) t72(11) t72(16); t72(2) t72(3) t72(5) t72(8) t72(12) t72(17); t72(4) t72(5) t72(6) t72(9) t72(13) t72(18); t72(7) t72(8) t72(9) t72(10) t72(14) t72(19); t72(11) t72(12) t72(13) t72(14) t72(15) t72(20); t72(16) t72(17) t72(18) t72(19) t72(20) t72(21);];
Mq  = res;
