% Calculate joint inertia matrix for
% S6RPRRPP7
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:13
% EndTime: 2019-03-09 04:50:20
% DurationCPUTime: 3.00s
% Computational Cost: add. (3117->395), mult. (7566->563), div. (0->0), fcn. (8196->6), ass. (0->187)
t176 = cos(qJ(3));
t253 = Icges(4,5) * t176;
t173 = sin(qJ(3));
t252 = Icges(4,6) * t173;
t251 = t253 / 0.2e1 - t252 / 0.2e1;
t250 = rSges(7,1) + pkin(5);
t249 = rSges(7,3) + qJ(6);
t172 = sin(qJ(4));
t175 = cos(qJ(4));
t109 = Icges(6,6) * t173 + (Icges(6,5) * t175 + Icges(6,3) * t172) * t176;
t110 = Icges(5,3) * t173 + (Icges(5,5) * t175 - Icges(5,6) * t172) * t176;
t113 = -Icges(7,6) * t173 + (Icges(7,4) * t175 + Icges(7,2) * t172) * t176;
t114 = Icges(6,2) * t173 + (Icges(6,4) * t175 + Icges(6,6) * t172) * t176;
t118 = -Icges(7,5) * t173 + (Icges(7,1) * t175 + Icges(7,4) * t172) * t176;
t119 = Icges(6,4) * t173 + (Icges(6,1) * t175 + Icges(6,5) * t172) * t176;
t120 = Icges(5,5) * t173 + (Icges(5,1) * t175 - Icges(5,4) * t172) * t176;
t217 = t175 * t176;
t222 = t172 * t176;
t248 = (t109 + t113) * t222 + (t110 + t114) * t173 + (t118 + t119 + t120) * t217;
t205 = m(6) / 0.2e1 + m(7) / 0.2e1;
t247 = 0.2e1 * t205;
t177 = cos(qJ(1));
t246 = (rSges(4,1) * t173 + rSges(4,2) * t176) * t177;
t174 = sin(qJ(1));
t214 = t177 * t175;
t220 = t174 * t172;
t139 = t173 * t220 - t214;
t215 = t177 * t172;
t219 = t174 * t175;
t140 = t173 * t219 + t215;
t218 = t174 * t176;
t63 = Icges(7,5) * t140 + Icges(7,6) * t139 + Icges(7,3) * t218;
t69 = Icges(7,4) * t140 + Icges(7,2) * t139 + Icges(7,6) * t218;
t75 = Icges(7,1) * t140 + Icges(7,4) * t139 + Icges(7,5) * t218;
t19 = t139 * t69 + t140 * t75 + t218 * t63;
t141 = t173 * t215 + t219;
t143 = -t173 * t214 + t220;
t216 = t176 * t177;
t64 = Icges(7,5) * t143 - Icges(7,6) * t141 - Icges(7,3) * t216;
t70 = Icges(7,4) * t143 - Icges(7,2) * t141 - Icges(7,6) * t216;
t76 = Icges(7,1) * t143 - Icges(7,4) * t141 - Icges(7,5) * t216;
t20 = t139 * t70 + t140 * t76 + t218 * t64;
t65 = Icges(6,5) * t140 - Icges(6,6) * t218 + Icges(6,3) * t139;
t71 = Icges(6,4) * t140 - Icges(6,2) * t218 + Icges(6,6) * t139;
t77 = Icges(6,1) * t140 - Icges(6,4) * t218 + Icges(6,5) * t139;
t21 = t139 * t65 + t140 * t77 - t218 * t71;
t66 = Icges(6,5) * t143 + Icges(6,6) * t216 - Icges(6,3) * t141;
t72 = Icges(6,4) * t143 + Icges(6,2) * t216 - Icges(6,6) * t141;
t78 = Icges(6,1) * t143 + Icges(6,4) * t216 - Icges(6,5) * t141;
t22 = t139 * t66 + t140 * t78 - t218 * t72;
t67 = Icges(5,5) * t140 - Icges(5,6) * t139 - Icges(5,3) * t218;
t73 = Icges(5,4) * t140 - Icges(5,2) * t139 - Icges(5,6) * t218;
t79 = Icges(5,1) * t140 - Icges(5,4) * t139 - Icges(5,5) * t218;
t23 = -t139 * t73 + t140 * t79 - t218 * t67;
t68 = Icges(5,5) * t143 + Icges(5,6) * t141 + Icges(5,3) * t216;
t74 = Icges(5,4) * t143 + Icges(5,2) * t141 + Icges(5,6) * t216;
t80 = Icges(5,1) * t143 + Icges(5,4) * t141 + Icges(5,5) * t216;
t24 = -t139 * t74 + t140 * t80 - t218 * t68;
t108 = -Icges(7,3) * t173 + (Icges(7,5) * t175 + Icges(7,6) * t172) * t176;
t44 = t108 * t218 + t139 * t113 + t140 * t118;
t45 = t139 * t109 - t114 * t218 + t140 * t119;
t115 = Icges(5,6) * t173 + (Icges(5,4) * t175 - Icges(5,2) * t172) * t176;
t46 = -t110 * t218 - t139 * t115 + t140 * t120;
t245 = ((t20 + t22 + t24) * t177 + (-t19 - t21 - t23) * t174) * t176 + (t44 + t45 + t46) * t173;
t25 = -t141 * t69 + t143 * t75 - t216 * t63;
t26 = -t141 * t70 + t143 * t76 - t216 * t64;
t27 = -t141 * t65 + t143 * t77 + t216 * t71;
t28 = -t141 * t66 + t143 * t78 + t216 * t72;
t29 = t141 * t73 + t143 * t79 + t216 * t67;
t30 = t141 * t74 + t143 * t80 + t216 * t68;
t47 = -t108 * t216 - t141 * t113 + t143 * t118;
t48 = -t141 * t109 + t114 * t216 + t143 * t119;
t49 = t110 * t216 + t141 * t115 + t143 * t120;
t244 = ((t26 + t28 + t30) * t177 + (-t25 - t27 - t29) * t174) * t176 + (t47 + t48 + t49) * t173;
t31 = -t173 * t63 + (t172 * t69 + t175 * t75) * t176;
t33 = t173 * t71 + (t172 * t65 + t175 * t77) * t176;
t35 = t173 * t67 + (-t172 * t73 + t175 * t79) * t176;
t243 = t31 + t33 + t35;
t32 = -t173 * t64 + (t172 * t70 + t175 * t76) * t176;
t34 = t173 * t72 + (t172 * t66 + t175 * t78) * t176;
t36 = t173 * t68 + (-t172 * t74 + t175 * t80) * t176;
t242 = t32 + t34 + t36;
t169 = t174 ^ 2;
t171 = t177 ^ 2;
t241 = -pkin(1) - pkin(7);
t238 = t174 / 0.2e1;
t236 = t177 / 0.2e1;
t153 = t176 * rSges(4,1) - t173 * rSges(4,2);
t235 = m(4) * t153;
t234 = m(7) * t176;
t210 = t140 * rSges(6,1) + t139 * rSges(6,3);
t82 = -rSges(6,2) * t218 + t210;
t93 = t140 * pkin(4) + t139 * qJ(5);
t233 = -t82 - t93;
t229 = t141 * rSges(6,3);
t85 = t143 * rSges(6,1) + rSges(6,2) * t216 - t229;
t130 = t141 * qJ(5);
t94 = t143 * pkin(4) - t130;
t232 = -t85 - t94;
t230 = t141 * rSges(7,2);
t228 = t139 * rSges(7,2) + t250 * t140 + t249 * t218;
t227 = t250 * t143 - t249 * t216 - t230;
t144 = (pkin(4) * t175 + qJ(5) * t172) * t176;
t226 = t144 * t218 + t173 * t93;
t162 = t177 * t173 * pkin(3);
t129 = t177 * (pkin(8) * t216 - t162);
t225 = t177 * t94 + t129;
t221 = t173 * t174;
t212 = (rSges(7,1) * t175 + rSges(7,2) * t172) * t176 + pkin(5) * t217 - t249 * t173;
t155 = t176 * pkin(3) + t173 * pkin(8);
t146 = t174 * t155;
t211 = t174 * t144 + t146;
t209 = t140 * rSges(5,1) - t139 * rSges(5,2);
t208 = -t144 - t155;
t207 = t177 * pkin(1) + t174 * qJ(2);
t206 = t169 + t171;
t204 = (-t173 * t108 - t115 * t222 + t248) * t173;
t203 = pkin(8) * t218;
t202 = -t93 - t228;
t201 = -t94 - t227;
t199 = rSges(4,1) * t221 + rSges(4,2) * t218 + t177 * rSges(4,3);
t198 = t177 * pkin(7) + t207;
t197 = (-rSges(6,2) - pkin(8)) * t176;
t196 = (-rSges(5,3) - pkin(8)) * t176;
t195 = t212 * t176;
t161 = pkin(3) * t221;
t194 = t161 + t198;
t192 = -t143 * rSges(5,1) - t141 * rSges(5,2);
t17 = t173 * t228 + t174 * t195 + t226;
t107 = t144 * t216;
t18 = t173 * t201 + t177 * t195 + t107;
t191 = -t17 * t177 + t18 * t174;
t56 = t174 * t212 + t211;
t57 = (t208 - t212) * t177;
t190 = t56 * t174 - t57 * t177;
t187 = Icges(4,5) * t173 + Icges(4,6) * t176;
t184 = t139 * t174 + t141 * t177;
t165 = t177 * qJ(2);
t183 = t174 * t241 + t162 + t165;
t181 = t130 + t183;
t39 = t230 + (-pkin(8) + t249) * t216 + (-pkin(4) - t250) * t143 + t181;
t180 = t194 + t93;
t40 = t180 - t203 + t228;
t182 = m(7) * (t174 * t39 - t177 * t40);
t179 = -t35 / 0.2e1 - t33 / 0.2e1 - t31 / 0.2e1 - t46 / 0.2e1 - t45 / 0.2e1 - t44 / 0.2e1;
t178 = t49 / 0.2e1 + t48 / 0.2e1 + t36 / 0.2e1 + t47 / 0.2e1 + t34 / 0.2e1 + t32 / 0.2e1;
t170 = t176 ^ 2;
t154 = t177 * rSges(2,1) - t174 * rSges(2,2);
t152 = -t174 * rSges(2,1) - t177 * rSges(2,2);
t149 = -t252 + t253;
t145 = t161 - t203;
t127 = -t177 * rSges(3,2) + t174 * rSges(3,3) + t207;
t126 = t177 * rSges(3,3) + t165 + (rSges(3,2) - pkin(1)) * t174;
t125 = t173 * rSges(5,3) + (rSges(5,1) * t175 - rSges(5,2) * t172) * t176;
t124 = t173 * rSges(6,2) + (rSges(6,1) * t175 + rSges(6,3) * t172) * t176;
t112 = Icges(4,3) * t174 - t177 * t187;
t111 = Icges(4,3) * t177 + t174 * t187;
t96 = t198 + t199;
t95 = t165 + t246 + (-rSges(4,3) + t241) * t174;
t92 = (-t125 - t155) * t177;
t91 = t174 * t125 + t146;
t86 = rSges(5,3) * t216 - t192;
t83 = -rSges(5,3) * t218 + t209;
t62 = -t174 * t199 + (t174 * rSges(4,3) - t246) * t177;
t61 = (-t124 + t208) * t177;
t60 = t174 * t124 + t211;
t59 = t174 * t196 + t194 + t209;
t58 = t177 * t196 + t183 + t192;
t55 = t125 * t216 - t173 * t86;
t54 = t125 * t218 + t173 * t83;
t50 = (-t174 * t86 - t177 * t83) * t176;
t43 = t174 * t197 + t180 + t210;
t42 = t229 + t177 * t197 + (-rSges(6,1) - pkin(4)) * t143 + t181;
t41 = t177 * t86 + t129 + (-t145 - t83) * t174;
t38 = t124 * t216 + t173 * t232 + t107;
t37 = t124 * t218 + t173 * t82 + t226;
t16 = (t174 * t232 + t177 * t233) * t176;
t15 = t177 * t85 + (-t145 + t233) * t174 + t225;
t14 = (t174 * t201 + t177 * t202) * t176;
t13 = t227 * t177 + (-t145 + t202) * t174 + t225;
t12 = t30 * t174 + t29 * t177;
t11 = t28 * t174 + t27 * t177;
t10 = t26 * t174 + t25 * t177;
t9 = t24 * t174 + t23 * t177;
t8 = t22 * t174 + t21 * t177;
t7 = t20 * t174 + t19 * t177;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t176 - t172 * t115) * t176 + m(6) * (t42 ^ 2 + t43 ^ 2) + m(7) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t95 ^ 2 + t96 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + m(2) * (t152 ^ 2 + t154 ^ 2) + t248 + (-0.2e1 * Icges(4,4) * t176 + Icges(4,2) * t173 - t108) * t173; m(6) * (t174 * t42 - t177 * t43) + t182 + m(5) * (t174 * t58 - t177 * t59) + m(4) * (t174 * t95 - t177 * t96) + m(3) * (t174 * t126 - t177 * t127); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t205) * t206; m(6) * (t42 * t60 + t43 * t61) + m(7) * (t39 * t56 + t40 * t57) + m(5) * (t58 * t91 + t59 * t92) + (t149 * t236 + t251 * t177 - t235 * t96 - t179) * t177 + (t149 * t238 + t251 * t174 + t235 * t95 + t178) * t174; m(5) * (t91 * t174 - t92 * t177) + m(6) * (t60 * t174 - t61 * t177) + m(7) * t190 + t206 * t235; m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(7) * (t13 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t41 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t153 ^ 2 * t206 + t62 ^ 2) + (t171 * t111 + t7 + t8 + t9) * t177 + (t169 * t112 + t10 + t11 + t12 + (t174 * t111 + t177 * t112) * t177) * t174; m(6) * (t37 * t43 + t38 * t42) + m(7) * (t17 * t40 + t18 * t39) + m(5) * (t54 * t59 + t55 * t58) + (t174 * t179 + t177 * t178) * t176 + t204; m(5) * (t55 * t174 - t54 * t177) + m(6) * (t38 * t174 - t37 * t177) + m(7) * t191; m(6) * (t16 * t15 + t37 * t61 + t38 * t60) + m(7) * (t14 * t13 + t17 * t57 + t18 * t56) + m(5) * (t41 * t50 + t54 * t92 + t55 * t91) + ((t11 / 0.2e1 + t10 / 0.2e1 + t12 / 0.2e1) * t177 + (-t7 / 0.2e1 - t9 / 0.2e1 - t8 / 0.2e1) * t174) * t176 + (t242 * t174 + t243 * t177) * t173 / 0.2e1 + t244 * t238 + t245 * t236; t204 * t173 + m(6) * (t16 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(7) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(5) * (t50 ^ 2 + t54 ^ 2 + t55 ^ 2) + (t244 * t177 - t245 * t174 + (-t243 * t174 + t242 * t177) * t173) * t176; m(6) * (t139 * t42 - t141 * t43) + m(7) * (t139 * t39 - t141 * t40); t184 * t247; m(6) * (t139 * t60 - t141 * t61 + t15 * t222) + m(7) * (t13 * t222 + t139 * t56 - t141 * t57); m(6) * (t139 * t38 - t141 * t37 + t16 * t222) + m(7) * (t139 * t18 + t14 * t222 - t141 * t17); (t170 * t172 ^ 2 + t139 ^ 2 + t141 ^ 2) * t247; t176 * t182; t206 * t234; m(7) * (-t173 * t13 + t176 * t190); m(7) * (-t173 * t14 + t176 * t191); (-t172 * t173 + t184) * t234; m(7) * (t170 * t206 + t173 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
