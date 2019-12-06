% Calculate joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:37
% DurationCPUTime: 3.02s
% Computational Cost: add. (6673->336), mult. (17238->514), div. (0->0), fcn. (21768->10), ass. (0->156)
t188 = Icges(3,1) + Icges(4,2);
t186 = Icges(3,4) + Icges(4,6);
t185 = Icges(4,4) - Icges(3,5);
t187 = Icges(3,2) + Icges(4,3);
t184 = Icges(3,6) - Icges(4,5);
t183 = Icges(3,3) + Icges(4,1);
t138 = sin(pkin(9));
t140 = cos(pkin(9));
t144 = sin(qJ(2));
t141 = cos(pkin(5));
t146 = cos(qJ(2));
t157 = t141 * t146;
t128 = t138 * t144 - t140 * t157;
t158 = t141 * t144;
t129 = t138 * t146 + t140 * t158;
t139 = sin(pkin(5));
t162 = t139 * t140;
t182 = t187 * t128 - t186 * t129 + t184 * t162;
t181 = t186 * t128 - t188 * t129 - t185 * t162;
t130 = t138 * t157 + t140 * t144;
t131 = -t138 * t158 + t140 * t146;
t163 = t138 * t139;
t180 = t187 * t130 - t186 * t131 - t184 * t163;
t179 = -t186 * t130 + t188 * t131 - t185 * t163;
t178 = t184 * t128 + t185 * t129 + t183 * t162;
t177 = -t184 * t130 - t185 * t131 + t183 * t163;
t176 = t183 * t141 + (-t185 * t144 + t184 * t146) * t139;
t175 = t184 * t141 + (t186 * t144 + t187 * t146) * t139;
t174 = t185 * t141 + (-t188 * t144 - t186 * t146) * t139;
t143 = sin(qJ(4));
t161 = t139 * t143;
t167 = cos(qJ(4));
t111 = -t130 * t167 + t138 * t161;
t113 = t128 * t167 + t140 * t161;
t151 = t139 * t167;
t132 = t141 * t143 + t146 * t151;
t112 = t130 * t143 + t138 * t151;
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t77 = -t112 * t142 + t131 * t145;
t78 = t112 * t145 + t131 * t142;
t50 = Icges(6,5) * t78 + Icges(6,6) * t77 + Icges(6,3) * t111;
t52 = Icges(6,4) * t78 + Icges(6,2) * t77 + Icges(6,6) * t111;
t54 = Icges(6,1) * t78 + Icges(6,4) * t77 + Icges(6,5) * t111;
t17 = t111 * t50 + t52 * t77 + t54 * t78;
t114 = t128 * t143 - t140 * t151;
t79 = -t114 * t142 + t129 * t145;
t80 = t114 * t145 + t129 * t142;
t51 = Icges(6,5) * t80 + Icges(6,6) * t79 - Icges(6,3) * t113;
t53 = Icges(6,4) * t80 + Icges(6,2) * t79 - Icges(6,6) * t113;
t55 = Icges(6,1) * t80 + Icges(6,4) * t79 - Icges(6,5) * t113;
t18 = t111 * t51 + t53 * t77 + t55 * t78;
t159 = t139 * t146;
t133 = t141 * t167 - t143 * t159;
t160 = t139 * t144;
t115 = -t133 * t142 + t145 * t160;
t116 = t133 * t145 + t142 * t160;
t69 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t132;
t70 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t132;
t71 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t132;
t26 = t111 * t69 + t70 * t77 + t71 * t78;
t1 = t111 * t17 - t113 * t18 + t132 * t26;
t172 = t1 / 0.2e1;
t21 = t115 * t52 + t116 * t54 + t132 * t50;
t22 = t115 * t53 + t116 * t55 + t132 * t51;
t36 = t115 * t70 + t116 * t71 + t132 * t69;
t7 = t111 * t21 - t113 * t22 + t132 * t36;
t171 = t7 / 0.2e1;
t170 = t111 / 0.2e1;
t169 = -t113 / 0.2e1;
t168 = t132 / 0.2e1;
t56 = rSges(6,1) * t78 + rSges(6,2) * t77 + rSges(6,3) * t111;
t166 = pkin(4) * t112 + pkin(8) * t111 + t56;
t57 = rSges(6,1) * t80 + rSges(6,2) * t79 - rSges(6,3) * t113;
t165 = pkin(4) * t114 - pkin(8) * t113 + t57;
t72 = rSges(6,1) * t116 + rSges(6,2) * t115 + rSges(6,3) * t132;
t164 = pkin(4) * t133 + pkin(8) * t132 + t72;
t105 = pkin(2) * t129 + qJ(3) * t128;
t106 = pkin(2) * t131 + qJ(3) * t130;
t156 = t105 * t163 + t106 * t162;
t104 = t141 * t106;
t117 = pkin(3) * t163 + pkin(7) * t131;
t155 = t141 * t117 + t104;
t118 = -pkin(3) * t162 + pkin(7) * t129;
t154 = -t105 - t118;
t134 = (pkin(2) * t144 - qJ(3) * t146) * t139;
t153 = -pkin(3) * t141 - pkin(7) * t160 - t134;
t152 = -m(4) - m(5) - m(6);
t150 = (-t141 * rSges(4,1) - (-rSges(4,2) * t144 - rSges(4,3) * t146) * t139 - t134) * t139;
t149 = t117 * t162 + t118 * t163 + t156;
t100 = rSges(5,1) * t133 - rSges(5,2) * t132 + rSges(5,3) * t160;
t148 = (-t100 + t153) * t139;
t147 = (t153 - t164) * t139;
t125 = t141 * rSges(3,3) + (rSges(3,1) * t144 + rSges(3,2) * t146) * t139;
t99 = Icges(5,1) * t133 - Icges(5,4) * t132 + Icges(5,5) * t160;
t98 = Icges(5,4) * t133 - Icges(5,2) * t132 + Icges(5,6) * t160;
t97 = Icges(5,5) * t133 - Icges(5,6) * t132 + Icges(5,3) * t160;
t96 = rSges(3,1) * t131 - rSges(3,2) * t130 + rSges(3,3) * t163;
t95 = rSges(3,1) * t129 - rSges(3,2) * t128 - rSges(3,3) * t162;
t94 = -rSges(4,1) * t162 - rSges(4,2) * t129 + rSges(4,3) * t128;
t93 = rSges(4,1) * t163 - rSges(4,2) * t131 + rSges(4,3) * t130;
t74 = -t125 * t162 - t141 * t95;
t73 = -t125 * t163 + t141 * t96;
t68 = rSges(5,1) * t114 + rSges(5,2) * t113 + rSges(5,3) * t129;
t67 = rSges(5,1) * t112 - rSges(5,2) * t111 + rSges(5,3) * t131;
t66 = Icges(5,1) * t114 + Icges(5,4) * t113 + Icges(5,5) * t129;
t65 = Icges(5,1) * t112 - Icges(5,4) * t111 + Icges(5,5) * t131;
t64 = Icges(5,4) * t114 + Icges(5,2) * t113 + Icges(5,6) * t129;
t63 = Icges(5,4) * t112 - Icges(5,2) * t111 + Icges(5,6) * t131;
t62 = Icges(5,5) * t114 + Icges(5,6) * t113 + Icges(5,3) * t129;
t61 = Icges(5,5) * t112 - Icges(5,6) * t111 + Icges(5,3) * t131;
t60 = (t138 * t95 + t140 * t96) * t139;
t59 = (-t105 - t94) * t141 + t140 * t150;
t58 = t138 * t150 + t141 * t93 + t104;
t49 = -t100 * t131 + t67 * t160;
t48 = t100 * t129 - t68 * t160;
t47 = -t132 * t98 + t133 * t99 + t97 * t160;
t46 = (t138 * t94 + t140 * t93) * t139 + t156;
t45 = -t129 * t67 + t131 * t68;
t44 = t113 * t98 + t114 * t99 + t129 * t97;
t43 = -t111 * t98 + t112 * t99 + t131 * t97;
t42 = (-t68 + t154) * t141 + t140 * t148;
t41 = t138 * t148 + t141 * t67 + t155;
t40 = -t113 * t72 - t132 * t57;
t39 = -t111 * t72 + t132 * t56;
t38 = -t132 * t64 + t133 * t66 + t62 * t160;
t37 = -t132 * t63 + t133 * t65 + t61 * t160;
t35 = t113 * t64 + t114 * t66 + t129 * t62;
t34 = t113 * t63 + t114 * t65 + t129 * t61;
t33 = -t111 * t64 + t112 * t66 + t131 * t62;
t32 = -t111 * t63 + t112 * t65 + t131 * t61;
t31 = (t138 * t68 + t140 * t67) * t139 + t149;
t30 = t111 * t57 + t113 * t56;
t29 = -t164 * t131 + t166 * t160;
t28 = t164 * t129 - t165 * t160;
t27 = -t113 * t69 + t70 * t79 + t71 * t80;
t25 = (t154 - t165) * t141 + t140 * t147;
t24 = t138 * t147 + t166 * t141 + t155;
t23 = -t166 * t129 + t165 * t131;
t20 = -t113 * t51 + t53 * t79 + t55 * t80;
t19 = -t113 * t50 + t52 * t79 + t54 * t80;
t16 = (t165 * t138 + t166 * t140) * t139 + t149;
t15 = t141 * t47 + (t138 * t37 - t140 * t38) * t139;
t14 = t129 * t38 + t131 * t37 + t47 * t160;
t13 = t141 * t44 + (t138 * t34 - t140 * t35) * t139;
t12 = t141 * t43 + (t138 * t32 - t140 * t33) * t139;
t11 = t129 * t35 + t131 * t34 + t44 * t160;
t10 = t129 * t33 + t131 * t32 + t43 * t160;
t9 = t141 * t36 + (t138 * t21 - t140 * t22) * t139;
t8 = t129 * t22 + t131 * t21 + t36 * t160;
t6 = t141 * t27 + (t138 * t19 - t140 * t20) * t139;
t5 = t141 * t26 + (t138 * t17 - t140 * t18) * t139;
t4 = t129 * t20 + t131 * t19 + t27 * t160;
t3 = t129 * t18 + t131 * t17 + t26 * t160;
t2 = t111 * t19 - t113 * t20 + t132 * t27;
t75 = [m(2) + m(3) - t152; m(3) * t60 + m(4) * t46 + m(5) * t31 + m(6) * t16; m(6) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t31 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(4) * (t46 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(3) * (t60 ^ 2 + t73 ^ 2 + t74 ^ 2) + (t15 + t9 + t176 * t141 ^ 2 + ((-t174 * t144 + t175 * t146) * t141 + (t178 * t141 + (t181 * t144 + t182 * t146) * t139) * t140 + (t177 * t141 + (t179 * t144 - t180 * t146) * t139) * t138) * t139) * t141 + (t5 + t12 + (t180 * t130 + t179 * t131 + t177 * t163) * t163 + (-t175 * t130 - t174 * t131 + t176 * t163) * t141) * t163 + (-t6 - t13 + (t182 * t128 - t181 * t129 + t178 * t162) * t162 + (t175 * t128 + t174 * t129 + t176 * t162) * t141 + (-t180 * t128 - t179 * t129 - t182 * t130 + t181 * t131 + t177 * t162 + t178 * t163) * t163) * t162; t152 * t159; m(6) * (t128 * t24 + t130 * t25 - t16 * t159) + m(5) * (t128 * t41 + t130 * t42 - t31 * t159) + m(4) * (t128 * t58 + t130 * t59 - t46 * t159); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t139 ^ 2 * t146 ^ 2 + t128 ^ 2 + t130 ^ 2); m(5) * t45 + m(6) * t23; (t8 / 0.2e1 + t14 / 0.2e1) * t141 + (t5 / 0.2e1 + t12 / 0.2e1) * t131 + (t6 / 0.2e1 + t13 / 0.2e1) * t129 + m(6) * (t16 * t23 + t24 * t29 + t25 * t28) + m(5) * (t31 * t45 + t41 * t49 + t42 * t48) + ((t9 / 0.2e1 + t15 / 0.2e1) * t144 + (-t4 / 0.2e1 - t11 / 0.2e1) * t140 + (t3 / 0.2e1 + t10 / 0.2e1) * t138) * t139; m(5) * (t49 * t128 + t48 * t130 - t45 * t159) + m(6) * (t29 * t128 + t28 * t130 - t23 * t159); (t14 + t8) * t160 + (t3 + t10) * t131 + (t4 + t11) * t129 + m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t45 ^ 2 + t48 ^ 2 + t49 ^ 2); m(6) * t30; m(6) * (t16 * t30 + t24 * t39 + t25 * t40) + t9 * t168 + t141 * t171 + t6 * t169 + t5 * t170 + (t138 * t172 - t140 * t2 / 0.2e1) * t139; m(6) * (t39 * t128 + t40 * t130 - t30 * t159); m(6) * (t23 * t30 + t28 * t40 + t29 * t39) + t4 * t169 + t3 * t170 + t131 * t172 + t129 * t2 / 0.2e1 + t8 * t168 + t160 * t171; m(6) * (t30 ^ 2 + t39 ^ 2 + t40 ^ 2) + t111 * t1 - t113 * t2 + t132 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t75(1), t75(2), t75(4), t75(7), t75(11); t75(2), t75(3), t75(5), t75(8), t75(12); t75(4), t75(5), t75(6), t75(9), t75(13); t75(7), t75(8), t75(9), t75(10), t75(14); t75(11), t75(12), t75(13), t75(14), t75(15);];
Mq = res;
