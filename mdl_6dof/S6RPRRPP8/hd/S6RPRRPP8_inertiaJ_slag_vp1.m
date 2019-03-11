% Calculate joint inertia matrix for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:54:00
% EndTime: 2019-03-09 04:54:07
% DurationCPUTime: 3.15s
% Computational Cost: add. (3129->402), mult. (7593->576), div. (0->0), fcn. (8237->6), ass. (0->183)
t175 = cos(qJ(3));
t251 = Icges(4,5) * t175;
t172 = sin(qJ(3));
t250 = Icges(4,6) * t172;
t249 = rSges(7,3) + qJ(6);
t248 = t251 / 0.2e1 - t250 / 0.2e1;
t247 = rSges(7,1) + pkin(5);
t171 = sin(qJ(4));
t174 = cos(qJ(4));
t108 = Icges(5,3) * t172 + (Icges(5,5) * t174 - Icges(5,6) * t171) * t175;
t114 = Icges(5,5) * t172 + (Icges(5,1) * t174 - Icges(5,4) * t171) * t175;
t117 = Icges(7,5) * t172 + (Icges(7,6) * t171 + Icges(7,3) * t174) * t175;
t118 = Icges(6,5) * t172 + (-Icges(6,6) * t174 + Icges(6,3) * t171) * t175;
t119 = Icges(7,4) * t172 + (Icges(7,2) * t171 + Icges(7,6) * t174) * t175;
t121 = Icges(7,1) * t172 + (Icges(7,4) * t171 + Icges(7,5) * t174) * t175;
t122 = Icges(6,1) * t172 + (-Icges(6,4) * t174 + Icges(6,5) * t171) * t175;
t212 = t174 * t175;
t217 = t171 * t175;
t245 = (t118 + t119) * t217 + (t114 + t117) * t212 + (t108 + t121 + t122) * t172;
t200 = m(6) / 0.2e1 + m(7) / 0.2e1;
t244 = 0.2e1 * t200;
t176 = cos(qJ(1));
t243 = (rSges(4,1) * t172 + rSges(4,2) * t175) * t176;
t209 = t176 * t174;
t173 = sin(qJ(1));
t215 = t173 * t171;
t138 = t172 * t215 - t209;
t210 = t176 * t171;
t214 = t173 * t174;
t139 = t172 * t214 + t210;
t242 = t138 * rSges(7,2) + t139 * t249;
t213 = t173 * t175;
t63 = -Icges(7,5) * t213 + Icges(7,6) * t138 + Icges(7,3) * t139;
t67 = -Icges(7,4) * t213 + Icges(7,2) * t138 + Icges(7,6) * t139;
t71 = -Icges(7,1) * t213 + Icges(7,4) * t138 + Icges(7,5) * t139;
t19 = t138 * t67 + t139 * t63 - t213 * t71;
t140 = t172 * t210 + t214;
t142 = t172 * t209 - t215;
t211 = t175 * t176;
t64 = Icges(7,5) * t211 - Icges(7,6) * t140 - Icges(7,3) * t142;
t68 = Icges(7,4) * t211 - Icges(7,2) * t140 - Icges(7,6) * t142;
t72 = Icges(7,1) * t211 - Icges(7,4) * t140 - Icges(7,5) * t142;
t20 = t138 * t68 + t139 * t64 - t213 * t72;
t65 = -Icges(6,5) * t213 - Icges(6,6) * t139 + Icges(6,3) * t138;
t69 = -Icges(6,4) * t213 - Icges(6,2) * t139 + Icges(6,6) * t138;
t73 = -Icges(6,1) * t213 - Icges(6,4) * t139 + Icges(6,5) * t138;
t21 = t138 * t65 - t139 * t69 - t213 * t73;
t66 = Icges(6,5) * t211 + Icges(6,6) * t142 - Icges(6,3) * t140;
t70 = Icges(6,4) * t211 + Icges(6,2) * t142 - Icges(6,6) * t140;
t74 = Icges(6,1) * t211 + Icges(6,4) * t142 - Icges(6,5) * t140;
t22 = t138 * t66 - t139 * t70 - t213 * t74;
t75 = Icges(5,5) * t139 - Icges(5,6) * t138 - Icges(5,3) * t213;
t77 = Icges(5,4) * t139 - Icges(5,2) * t138 - Icges(5,6) * t213;
t79 = Icges(5,1) * t139 - Icges(5,4) * t138 - Icges(5,5) * t213;
t27 = -t138 * t77 + t139 * t79 - t213 * t75;
t76 = -Icges(5,5) * t142 + Icges(5,6) * t140 + Icges(5,3) * t211;
t78 = -Icges(5,4) * t142 + Icges(5,2) * t140 + Icges(5,6) * t211;
t80 = -Icges(5,1) * t142 + Icges(5,4) * t140 + Icges(5,5) * t211;
t28 = -t138 * t78 + t139 * t80 - t213 * t76;
t44 = t139 * t117 + t138 * t119 - t121 * t213;
t120 = Icges(6,4) * t172 + (-Icges(6,2) * t174 + Icges(6,6) * t171) * t175;
t45 = t138 * t118 - t139 * t120 - t122 * t213;
t111 = Icges(5,6) * t172 + (Icges(5,4) * t174 - Icges(5,2) * t171) * t175;
t48 = -t108 * t213 - t138 * t111 + t139 * t114;
t241 = ((t20 + t22 + t28) * t176 + (-t19 - t21 - t27) * t173) * t175 + (t44 + t45 + t48) * t172;
t23 = -t140 * t67 - t142 * t63 + t211 * t71;
t24 = -t140 * t68 - t142 * t64 + t211 * t72;
t25 = -t140 * t65 + t142 * t69 + t211 * t73;
t26 = -t140 * t66 + t142 * t70 + t211 * t74;
t29 = t140 * t77 - t142 * t79 + t211 * t75;
t30 = t140 * t78 - t142 * t80 + t211 * t76;
t46 = -t142 * t117 - t140 * t119 + t121 * t211;
t47 = -t140 * t118 + t142 * t120 + t122 * t211;
t49 = t108 * t211 + t140 * t111 - t142 * t114;
t240 = ((t24 + t26 + t30) * t176 + (-t23 - t25 - t29) * t173) * t175 + (t46 + t47 + t49) * t172;
t31 = t172 * t75 + (-t171 * t77 + t174 * t79) * t175;
t33 = t172 * t71 + (t171 * t67 + t174 * t63) * t175;
t35 = t172 * t73 + (t171 * t65 - t174 * t69) * t175;
t239 = t31 + t33 + t35;
t32 = t172 * t76 + (-t171 * t78 + t174 * t80) * t175;
t34 = t172 * t72 + (t171 * t68 + t174 * t64) * t175;
t36 = t172 * t74 + (t171 * t66 - t174 * t70) * t175;
t238 = t32 + t34 + t36;
t168 = t173 ^ 2;
t170 = t176 ^ 2;
t237 = -pkin(1) - pkin(7);
t234 = t173 / 0.2e1;
t232 = t176 / 0.2e1;
t153 = t175 * rSges(4,1) - t172 * rSges(4,2);
t231 = m(4) * t153;
t206 = -t139 * rSges(6,2) + t138 * rSges(6,3);
t82 = -rSges(6,1) * t213 + t206;
t93 = t139 * pkin(4) + t138 * qJ(5);
t230 = -t82 - t93;
t227 = t140 * rSges(6,3);
t84 = rSges(6,1) * t211 + t142 * rSges(6,2) - t227;
t130 = t140 * qJ(5);
t94 = -t142 * pkin(4) - t130;
t229 = -t84 - t94;
t228 = t140 * rSges(7,2);
t226 = -t247 * t213 + t242;
t225 = -t142 * t249 + t247 * t211 - t228;
t144 = (pkin(4) * t174 + qJ(5) * t171) * t175;
t222 = t144 * t213 + t172 * t93;
t160 = t176 * t172 * pkin(3);
t129 = t176 * (pkin(8) * t211 - t160);
t221 = t176 * t94 + t129;
t216 = t172 * t173;
t208 = (rSges(7,2) * t171 + rSges(7,3) * t174) * t175 + qJ(6) * t212 + t247 * t172;
t155 = t175 * pkin(3) + t172 * pkin(8);
t146 = t173 * t155;
t207 = t173 * t144 + t146;
t204 = t139 * rSges(5,1) - t138 * rSges(5,2);
t203 = -t144 - t155;
t202 = t176 * pkin(1) + t173 * qJ(2);
t201 = t168 + t170;
t199 = (-t111 * t217 - t120 * t212 + t245) * t172;
t198 = -t93 - t226;
t197 = -t94 - t225;
t195 = rSges(4,1) * t216 + rSges(4,2) * t213 + t176 * rSges(4,3);
t194 = t176 * pkin(7) + t202;
t193 = (-rSges(6,1) - pkin(8)) * t175;
t192 = (-rSges(5,3) - pkin(8)) * t175;
t191 = t208 * t175;
t190 = (-pkin(8) - t247) * t175;
t159 = pkin(3) * t216;
t189 = t159 + t194;
t187 = t142 * rSges(5,1) - t140 * rSges(5,2);
t184 = Icges(4,5) * t172 + Icges(4,6) * t175;
t164 = t176 * qJ(2);
t181 = t173 * t237 + t160 + t164;
t180 = t130 + t181;
t179 = t189 + t93;
t178 = t32 / 0.2e1 + t49 / 0.2e1 + t36 / 0.2e1 + t47 / 0.2e1 + t46 / 0.2e1 + t34 / 0.2e1;
t177 = -t44 / 0.2e1 - t35 / 0.2e1 - t33 / 0.2e1 - t31 / 0.2e1 - t48 / 0.2e1 - t45 / 0.2e1;
t169 = t175 ^ 2;
t154 = t176 * rSges(2,1) - t173 * rSges(2,2);
t152 = -t173 * rSges(2,1) - t176 * rSges(2,2);
t149 = -t250 + t251;
t145 = -pkin(8) * t213 + t159;
t127 = -t176 * rSges(3,2) + t173 * rSges(3,3) + t202;
t126 = t176 * rSges(3,3) + t164 + (rSges(3,2) - pkin(1)) * t173;
t125 = t172 * rSges(6,1) + (-rSges(6,2) * t174 + rSges(6,3) * t171) * t175;
t123 = t172 * rSges(5,3) + (rSges(5,1) * t174 - rSges(5,2) * t171) * t175;
t110 = Icges(4,3) * t173 - t176 * t184;
t109 = Icges(4,3) * t176 + t173 * t184;
t107 = t144 * t211;
t96 = t194 + t195;
t95 = t164 + t243 + (-rSges(4,3) + t237) * t173;
t92 = (-t123 - t155) * t176;
t91 = t173 * t123 + t146;
t86 = rSges(5,3) * t211 - t187;
t85 = -rSges(5,3) * t213 + t204;
t62 = -t173 * t195 + (t173 * rSges(4,3) - t243) * t176;
t61 = (-t125 + t203) * t176;
t60 = t173 * t125 + t207;
t59 = t173 * t192 + t189 + t204;
t58 = t176 * t192 + t181 + t187;
t57 = (t203 - t208) * t176;
t56 = t173 * t208 + t207;
t55 = t123 * t211 - t172 * t86;
t54 = t123 * t213 + t172 * t85;
t50 = (-t173 * t86 - t176 * t85) * t175;
t43 = t173 * t193 + t179 + t206;
t42 = t227 + t176 * t193 + (-rSges(6,2) + pkin(4)) * t142 + t180;
t41 = t176 * t86 + t129 + (-t145 - t85) * t173;
t40 = t173 * t190 + t179 + t242;
t39 = t228 + t176 * t190 + (pkin(4) + t249) * t142 + t180;
t38 = t125 * t211 + t172 * t229 + t107;
t37 = t125 * t213 + t172 * t82 + t222;
t18 = t172 * t197 + t176 * t191 + t107;
t17 = t172 * t226 + t173 * t191 + t222;
t16 = (t173 * t229 + t176 * t230) * t175;
t15 = t176 * t84 + (-t145 + t230) * t173 + t221;
t14 = (t173 * t197 + t176 * t198) * t175;
t13 = t225 * t176 + (-t145 + t198) * t173 + t221;
t12 = t30 * t173 + t29 * t176;
t11 = t28 * t173 + t27 * t176;
t10 = t26 * t173 + t25 * t176;
t9 = t24 * t173 + t23 * t176;
t8 = t22 * t173 + t21 * t176;
t7 = t20 * t173 + t19 * t176;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t175 - t171 * t111 - t174 * t120) * t175 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t95 ^ 2 + t96 ^ 2) + m(2) * (t152 ^ 2 + t154 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + t245 + (-0.2e1 * Icges(4,4) * t175 + Icges(4,2) * t172) * t172; m(7) * (t173 * t39 - t176 * t40) + m(5) * (t173 * t58 - t176 * t59) + m(6) * (t173 * t42 - t176 * t43) + m(4) * (t173 * t95 - t176 * t96) + m(3) * (t173 * t126 - t176 * t127); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t200) * t201; m(7) * (t39 * t56 + t40 * t57) + m(5) * (t91 * t58 + t92 * t59) + m(6) * (t42 * t60 + t43 * t61) + (t149 * t232 + t176 * t248 - t231 * t96 - t177) * t176 + (t149 * t234 + t173 * t248 + t231 * t95 + t178) * t173; m(5) * (t91 * t173 - t92 * t176) + m(6) * (t60 * t173 - t61 * t176) + m(7) * (t56 * t173 - t57 * t176) + t201 * t231; m(7) * (t13 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t41 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t153 ^ 2 * t201 + t62 ^ 2) + (t170 * t109 + t11 + t7 + t8) * t176 + (t168 * t110 + t10 + t12 + t9 + (t173 * t109 + t176 * t110) * t176) * t173; m(7) * (t17 * t40 + t18 * t39) + m(5) * (t54 * t59 + t55 * t58) + m(6) * (t37 * t43 + t38 * t42) + (t173 * t177 + t176 * t178) * t175 + t199; m(5) * (t55 * t173 - t54 * t176) + m(6) * (t38 * t173 - t37 * t176) + m(7) * (-t17 * t176 + t18 * t173); m(7) * (t13 * t14 + t17 * t57 + t18 * t56) + m(6) * (t15 * t16 + t37 * t61 + t38 * t60) + m(5) * (t50 * t41 + t54 * t92 + t55 * t91) + ((t12 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t176 + (-t8 / 0.2e1 - t7 / 0.2e1 - t11 / 0.2e1) * t173) * t175 + (t238 * t173 + t239 * t176) * t172 / 0.2e1 + t240 * t234 + t241 * t232; t199 * t172 + m(7) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t16 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t50 ^ 2 + t54 ^ 2 + t55 ^ 2) + (t240 * t176 - t241 * t173 + (-t239 * t173 + t238 * t176) * t172) * t175; m(7) * (t138 * t39 - t140 * t40) + m(6) * (t138 * t42 - t140 * t43); (t138 * t173 + t140 * t176) * t244; m(7) * (t13 * t217 + t138 * t56 - t140 * t57) + m(6) * (t138 * t60 - t140 * t61 + t15 * t217); m(7) * (t138 * t18 + t14 * t217 - t140 * t17) + m(6) * (t138 * t38 - t140 * t37 + t16 * t217); (t169 * t171 ^ 2 + t138 ^ 2 + t140 ^ 2) * t244; m(7) * (t139 * t39 - t142 * t40); m(7) * (t139 * t173 + t142 * t176); m(7) * (t13 * t212 + t139 * t56 - t142 * t57); m(7) * (t139 * t18 + t14 * t212 - t142 * t17); m(7) * (t169 * t174 * t171 + t139 * t138 + t142 * t140); m(7) * (t169 * t174 ^ 2 + t139 ^ 2 + t142 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
