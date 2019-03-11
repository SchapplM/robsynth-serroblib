% Calculate joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:57
% EndTime: 2019-03-09 02:47:01
% DurationCPUTime: 2.48s
% Computational Cost: add. (4552->343), mult. (6334->512), div. (0->0), fcn. (7259->10), ass. (0->159)
t208 = Icges(5,1) + Icges(6,1);
t207 = -Icges(5,4) + Icges(6,5);
t206 = Icges(6,4) + Icges(5,5);
t205 = Icges(5,2) + Icges(6,3);
t204 = Icges(5,6) - Icges(6,6);
t137 = sin(qJ(1));
t203 = -t137 / 0.2e1;
t190 = t137 / 0.2e1;
t139 = cos(qJ(1));
t202 = t139 / 0.2e1;
t128 = pkin(9) + qJ(3);
t126 = cos(t128);
t133 = cos(pkin(10));
t171 = t139 * t133;
t131 = sin(pkin(10));
t174 = t137 * t131;
t102 = t126 * t174 + t171;
t172 = t139 * t131;
t173 = t137 * t133;
t103 = t126 * t173 - t172;
t125 = sin(t128);
t177 = t125 * t137;
t198 = t205 * t102 + t207 * t103 - t204 * t177;
t104 = t126 * t172 - t173;
t105 = t126 * t171 + t174;
t176 = t125 * t139;
t197 = t205 * t104 + t207 * t105 - t204 * t176;
t49 = Icges(5,5) * t103 - Icges(5,6) * t102 + Icges(5,3) * t177;
t51 = Icges(6,4) * t103 + Icges(6,2) * t177 + Icges(6,6) * t102;
t201 = t49 + t51;
t50 = Icges(5,5) * t105 - Icges(5,6) * t104 + Icges(5,3) * t176;
t52 = Icges(6,4) * t105 + Icges(6,2) * t176 + Icges(6,6) * t104;
t200 = t52 + t50;
t196 = t207 * t102 + t208 * t103 + t206 * t177;
t195 = t207 * t104 + t208 * t105 + t206 * t176;
t191 = m(7) / 0.2e1;
t192 = m(6) / 0.2e1;
t167 = t192 + t191;
t199 = 0.2e1 * t167;
t129 = t137 ^ 2;
t130 = t139 ^ 2;
t194 = 0.2e1 * t125;
t193 = m(5) / 0.2e1;
t189 = -t139 / 0.2e1;
t188 = -rSges(7,3) - pkin(8);
t187 = pkin(3) * t126;
t136 = sin(qJ(6));
t138 = cos(qJ(6));
t69 = t104 * t138 - t105 * t136;
t70 = t104 * t136 + t105 * t138;
t186 = t70 * rSges(7,1) + t69 * rSges(7,2);
t175 = t126 * t139;
t169 = pkin(3) * t175 + qJ(4) * t176;
t185 = t129 * (qJ(4) * t125 + t187) + t139 * t169;
t184 = t102 * rSges(6,3);
t183 = rSges(3,3) + qJ(2);
t111 = t125 * pkin(3) - t126 * qJ(4);
t182 = -t111 + t126 * rSges(5,3) - (rSges(5,1) * t133 - rSges(5,2) * t131) * t125;
t181 = -t111 - (pkin(4) * t133 + qJ(5) * t131) * t125;
t180 = Icges(4,4) * t125;
t179 = Icges(4,4) * t126;
t178 = t125 * t131;
t135 = -pkin(7) - qJ(2);
t170 = t139 * t135;
t168 = t129 + t130;
t89 = (t131 * t138 - t133 * t136) * t125;
t90 = (t131 * t136 + t133 * t138) * t125;
t42 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t126;
t43 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t126;
t44 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t126;
t166 = t126 * t42 + t89 * t43 + t90 * t44;
t165 = t126 * rSges(6,2) - (rSges(6,1) * t133 + rSges(6,3) * t131) * t125 + t181;
t164 = t105 * rSges(6,1) + rSges(6,2) * t176 + t104 * rSges(6,3);
t163 = t105 * rSges(5,1) - t104 * rSges(5,2) + rSges(5,3) * t176;
t66 = t102 * t138 - t103 * t136;
t67 = t102 * t136 + t103 * t138;
t26 = Icges(7,5) * t67 + Icges(7,6) * t66 - Icges(7,3) * t177;
t28 = Icges(7,4) * t67 + Icges(7,2) * t66 - Icges(7,6) * t177;
t30 = Icges(7,1) * t67 + Icges(7,4) * t66 - Icges(7,5) * t177;
t10 = t126 * t26 + t89 * t28 + t90 * t30;
t12 = -t177 * t42 + t66 * t43 + t67 * t44;
t162 = -t10 / 0.2e1 - t12 / 0.2e1;
t27 = Icges(7,5) * t70 + Icges(7,6) * t69 - Icges(7,3) * t176;
t29 = Icges(7,4) * t70 + Icges(7,2) * t69 - Icges(7,6) * t176;
t31 = Icges(7,1) * t70 + Icges(7,4) * t69 - Icges(7,5) * t176;
t11 = t126 * t27 + t89 * t29 + t90 * t31;
t13 = -t176 * t42 + t69 * t43 + t70 * t44;
t161 = -t13 / 0.2e1 - t11 / 0.2e1;
t75 = -Icges(6,6) * t126 + (Icges(6,5) * t133 + Icges(6,3) * t131) * t125;
t78 = -Icges(5,6) * t126 + (Icges(5,4) * t133 - Icges(5,2) * t131) * t125;
t160 = t78 / 0.2e1 - t75 / 0.2e1;
t79 = -Icges(6,4) * t126 + (Icges(6,1) * t133 + Icges(6,5) * t131) * t125;
t80 = -Icges(5,5) * t126 + (Icges(5,1) * t133 - Icges(5,4) * t131) * t125;
t159 = t80 / 0.2e1 + t79 / 0.2e1;
t134 = cos(pkin(9));
t123 = t134 * pkin(2) + pkin(1);
t158 = -t123 - t187;
t155 = t105 * pkin(4) + t104 * qJ(5);
t94 = t102 * qJ(5);
t157 = t137 * (t103 * pkin(4) + t94) + t139 * t155 + t185;
t156 = -t94 - t170;
t154 = t139 * t123 - t137 * t135;
t45 = t90 * rSges(7,1) + t89 * rSges(7,2) + t126 * rSges(7,3);
t153 = -t125 * t133 * pkin(5) - t126 * pkin(8) + t181 - t45;
t152 = t193 + t167;
t151 = -t67 * rSges(7,1) - t66 * rSges(7,2);
t150 = rSges(4,1) * t126 - rSges(4,2) * t125;
t149 = -t103 * rSges(5,1) + t102 * rSges(5,2);
t146 = Icges(4,1) * t126 - t180;
t145 = -Icges(4,2) * t125 + t179;
t144 = Icges(4,5) * t126 - Icges(4,6) * t125;
t143 = rSges(4,1) * t175 - rSges(4,2) * t176 + t137 * rSges(4,3);
t132 = sin(pkin(9));
t142 = rSges(3,1) * t134 - rSges(3,2) * t132 + pkin(1);
t141 = t154 + t169;
t140 = t141 + t155;
t124 = t125 ^ 2;
t114 = t139 * rSges(2,1) - t137 * rSges(2,2);
t113 = -t137 * rSges(2,1) - t139 * rSges(2,2);
t112 = t125 * rSges(4,1) + t126 * rSges(4,2);
t108 = Icges(4,5) * t125 + Icges(4,6) * t126;
t100 = t105 * pkin(5);
t84 = Icges(4,3) * t137 + t139 * t144;
t83 = -Icges(4,3) * t139 + t137 * t144;
t74 = t137 * t183 + t139 * t142;
t73 = -t137 * t142 + t139 * t183;
t72 = t143 + t154;
t71 = (rSges(4,3) - t135) * t139 + (-t123 - t150) * t137;
t60 = t182 * t139;
t59 = t182 * t137;
t46 = t139 * t143 + (-t139 * rSges(4,3) + t137 * t150) * t137;
t39 = t165 * t139;
t38 = t165 * t137;
t35 = t141 + t163;
t34 = -t170 + ((-rSges(5,3) - qJ(4)) * t125 + t158) * t137 + t149;
t33 = -rSges(7,3) * t176 + t186;
t32 = -rSges(7,3) * t177 - t151;
t25 = t153 * t139;
t24 = t153 * t137;
t23 = t140 + t164;
t22 = -t184 + (-rSges(6,1) - pkin(4)) * t103 + ((-rSges(6,2) - qJ(4)) * t125 + t158) * t137 + t156;
t21 = t137 * (rSges(5,3) * t177 - t149) + t139 * t163 + t185;
t20 = t126 * t33 + t176 * t45;
t19 = -t126 * t32 - t177 * t45;
t18 = t176 * t188 + t100 + t140 + t186;
t17 = (-pkin(4) - pkin(5)) * t103 + ((-qJ(4) - t188) * t125 + t158) * t137 + t151 + t156;
t16 = (t137 * t33 - t139 * t32) * t125;
t15 = t166 * t126;
t14 = t137 * (t103 * rSges(6,1) + rSges(6,2) * t177 + t184) + t139 * t164 + t157;
t9 = -t176 * t27 + t69 * t29 + t70 * t31;
t8 = -t176 * t26 + t69 * t28 + t70 * t30;
t7 = -t177 * t27 + t66 * t29 + t67 * t31;
t6 = -t177 * t26 + t66 * t28 + t67 * t30;
t5 = (-pkin(8) * t176 + t100 + t33) * t139 + (t103 * pkin(5) - pkin(8) * t177 + t32) * t137 + t157;
t4 = t9 * t137 - t8 * t139;
t3 = t7 * t137 - t6 * t139;
t2 = t13 * t126 + (-t137 * t8 - t139 * t9) * t125;
t1 = t12 * t126 + (-t137 * t6 - t139 * t7) * t125;
t36 = [Icges(3,2) * t134 ^ 2 + Icges(2,3) + (Icges(3,1) * t132 + 0.2e1 * Icges(3,4) * t134) * t132 + (t180 + (t204 * t131 - t206 * t133) * t125 + (Icges(4,2) + Icges(5,3) + Icges(6,2)) * t126) * t126 + (Icges(4,1) * t125 + t179 + (t79 + t80) * t133 + (t75 - t78) * t131) * t125 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2) + m(4) * (t71 ^ 2 + t72 ^ 2) + m(3) * (t73 ^ 2 + t74 ^ 2) + m(2) * (t113 ^ 2 + t114 ^ 2) + t166; m(7) * (t137 * t17 - t139 * t18) + m(6) * (t137 * t22 - t139 * t23) + m(5) * (t137 * t34 - t139 * t35) + m(4) * (t137 * t71 - t139 * t72) + m(3) * (t137 * t73 - t139 * t74); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t152) * t168; (t160 * t102 - t159 * t103 + t108 * t202 + t162) * t139 + (-t160 * t104 + t159 * t105 + t108 * t190 - t161) * t137 + m(7) * (t17 * t25 + t18 * t24) + m(6) * (t22 * t39 + t23 * t38) + m(5) * (t34 * t60 + t35 * t59) + m(4) * (-t137 * t72 - t139 * t71) * t112 + ((Icges(4,6) * t202 + t145 * t203 + t49 / 0.2e1 + t51 / 0.2e1) * t139 + (Icges(4,6) * t190 + t145 * t202 - t50 / 0.2e1 - t52 / 0.2e1) * t137) * t126 + ((Icges(4,5) * t137 + t197 * t131 + t195 * t133 + t139 * t146) * t190 + (-Icges(4,5) * t139 + t198 * t131 + t196 * t133 + t137 * t146) * t189) * t125; m(5) * (t60 * t137 - t59 * t139) + m(6) * (t39 * t137 - t38 * t139) + m(7) * (t25 * t137 - t24 * t139); m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t14 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t21 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(4) * (t112 ^ 2 * t168 + t46 ^ 2) + (-t130 * t83 - t3 + (t198 * t102 + t196 * t103 + t201 * t177) * t139) * t139 + (t4 + t129 * t84 + (t197 * t104 + t195 * t105 + t200 * t176) * t137 + (-t197 * t102 - t195 * t103 - t198 * t104 - t196 * t105 - t137 * t83 + t139 * t84 - t201 * t176 - t200 * t177) * t139) * t137; ((t137 * t18 + t139 * t17) * t191 + (t137 * t23 + t139 * t22) * t192 + (t137 * t35 + t139 * t34) * t193) * t194; 0; m(7) * (-t126 * t5 + (t137 * t24 + t139 * t25) * t125) + m(6) * (-t126 * t14 + (t137 * t38 + t139 * t39) * t125) + m(5) * (-t126 * t21 + (t137 * t59 + t139 * t60) * t125); 0.2e1 * t152 * (t124 * t168 + t126 ^ 2); m(7) * (t102 * t18 + t104 * t17) + m(6) * (t102 * t23 + t104 * t22); (-t102 * t139 + t104 * t137) * t199; m(7) * (t102 * t24 + t104 * t25 + t178 * t5) + m(6) * (t102 * t38 + t104 * t39 + t14 * t178); t167 * (t102 * t137 + t104 * t139 - t126 * t131) * t194; (t124 * t131 ^ 2 + t102 ^ 2 + t104 ^ 2) * t199; m(7) * (t17 * t19 + t18 * t20) + t15 + (t137 * t162 + t139 * t161) * t125; m(7) * (t19 * t137 - t20 * t139); t126 * (-t10 * t139 + t11 * t137) / 0.2e1 + m(7) * (t16 * t5 + t19 * t25 + t20 * t24) + t1 * t189 + t2 * t190 + (t4 * t189 + t3 * t203) * t125; m(7) * (-t16 * t126 + (t137 * t20 + t139 * t19) * t125); m(7) * (t20 * t102 + t19 * t104 + t16 * t178); t126 * t15 + m(7) * (t16 ^ 2 + t19 ^ 2 + t20 ^ 2) + (-t139 * t2 - t137 * t1 + t126 * (-t10 * t137 - t11 * t139)) * t125;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t36(1) t36(2) t36(4) t36(7) t36(11) t36(16); t36(2) t36(3) t36(5) t36(8) t36(12) t36(17); t36(4) t36(5) t36(6) t36(9) t36(13) t36(18); t36(7) t36(8) t36(9) t36(10) t36(14) t36(19); t36(11) t36(12) t36(13) t36(14) t36(15) t36(20); t36(16) t36(17) t36(18) t36(19) t36(20) t36(21);];
Mq  = res;
