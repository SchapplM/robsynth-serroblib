% Calculate joint inertia matrix for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:44
% EndTime: 2019-03-09 04:42:51
% DurationCPUTime: 3.12s
% Computational Cost: add. (5778->400), mult. (7564->582), div. (0->0), fcn. (8204->8), ass. (0->185)
t166 = pkin(9) + qJ(3);
t164 = cos(t166);
t172 = sin(qJ(4));
t175 = cos(qJ(1));
t207 = t175 * t172;
t173 = sin(qJ(1));
t174 = cos(qJ(4));
t209 = t173 * t174;
t142 = t164 * t207 - t209;
t206 = t175 * t174;
t210 = t173 * t172;
t143 = t164 * t206 + t210;
t163 = sin(t166);
t213 = t163 * t175;
t243 = rSges(7,3) + qJ(6);
t244 = rSges(7,1) + pkin(5);
t220 = t142 * rSges(7,2) + t143 * t244 - t243 * t213;
t246 = Icges(4,5) * t163;
t245 = t246 / 0.2e1;
t112 = -Icges(5,3) * t164 + (Icges(5,5) * t174 - Icges(5,6) * t172) * t163;
t114 = -Icges(6,2) * t164 + (Icges(6,4) * t174 + Icges(6,6) * t172) * t163;
t242 = t112 + t114;
t198 = m(6) / 0.2e1 + m(7) / 0.2e1;
t241 = 0.2e1 * t198;
t140 = t164 * t210 + t206;
t141 = t164 * t209 - t207;
t215 = t163 * t173;
t63 = Icges(7,5) * t141 + Icges(7,6) * t140 - Icges(7,3) * t215;
t69 = Icges(7,4) * t141 + Icges(7,2) * t140 - Icges(7,6) * t215;
t75 = Icges(7,1) * t141 + Icges(7,4) * t140 - Icges(7,5) * t215;
t19 = t140 * t69 + t141 * t75 - t215 * t63;
t64 = Icges(7,5) * t143 + Icges(7,6) * t142 - Icges(7,3) * t213;
t70 = Icges(7,4) * t143 + Icges(7,2) * t142 - Icges(7,6) * t213;
t76 = Icges(7,1) * t143 + Icges(7,4) * t142 - Icges(7,5) * t213;
t20 = t140 * t70 + t141 * t76 - t215 * t64;
t65 = Icges(6,5) * t141 + Icges(6,6) * t215 + Icges(6,3) * t140;
t71 = Icges(6,4) * t141 + Icges(6,2) * t215 + Icges(6,6) * t140;
t77 = Icges(6,1) * t141 + Icges(6,4) * t215 + Icges(6,5) * t140;
t21 = t140 * t65 + t141 * t77 + t215 * t71;
t66 = Icges(6,5) * t143 + Icges(6,6) * t213 + Icges(6,3) * t142;
t72 = Icges(6,4) * t143 + Icges(6,2) * t213 + Icges(6,6) * t142;
t78 = Icges(6,1) * t143 + Icges(6,4) * t213 + Icges(6,5) * t142;
t22 = t140 * t66 + t141 * t78 + t215 * t72;
t67 = Icges(5,5) * t141 - Icges(5,6) * t140 + Icges(5,3) * t215;
t73 = Icges(5,4) * t141 - Icges(5,2) * t140 + Icges(5,6) * t215;
t79 = Icges(5,1) * t141 - Icges(5,4) * t140 + Icges(5,5) * t215;
t23 = -t140 * t73 + t141 * t79 + t215 * t67;
t68 = Icges(5,5) * t143 - Icges(5,6) * t142 + Icges(5,3) * t213;
t74 = Icges(5,4) * t143 - Icges(5,2) * t142 + Icges(5,6) * t213;
t80 = Icges(5,1) * t143 - Icges(5,4) * t142 + Icges(5,5) * t213;
t24 = -t140 * t74 + t141 * t80 + t215 * t68;
t110 = Icges(7,3) * t164 + (Icges(7,5) * t174 + Icges(7,6) * t172) * t163;
t113 = Icges(7,6) * t164 + (Icges(7,4) * t174 + Icges(7,2) * t172) * t163;
t116 = Icges(7,5) * t164 + (Icges(7,1) * t174 + Icges(7,4) * t172) * t163;
t42 = -t110 * t215 + t140 * t113 + t141 * t116;
t111 = -Icges(6,6) * t164 + (Icges(6,5) * t174 + Icges(6,3) * t172) * t163;
t117 = -Icges(6,4) * t164 + (Icges(6,1) * t174 + Icges(6,5) * t172) * t163;
t43 = t140 * t111 + t114 * t215 + t141 * t117;
t115 = -Icges(5,6) * t164 + (Icges(5,4) * t174 - Icges(5,2) * t172) * t163;
t118 = -Icges(5,5) * t164 + (Icges(5,1) * t174 - Icges(5,4) * t172) * t163;
t44 = t112 * t215 - t140 * t115 + t141 * t118;
t239 = (-t42 - t43 - t44) * t164 + ((t20 + t22 + t24) * t175 + (t19 + t21 + t23) * t173) * t163;
t25 = t142 * t69 + t143 * t75 - t213 * t63;
t26 = t142 * t70 + t143 * t76 - t213 * t64;
t27 = t142 * t65 + t143 * t77 + t213 * t71;
t28 = t142 * t66 + t143 * t78 + t213 * t72;
t29 = -t142 * t73 + t143 * t79 + t213 * t67;
t30 = -t142 * t74 + t143 * t80 + t213 * t68;
t45 = -t110 * t213 + t142 * t113 + t143 * t116;
t46 = t142 * t111 + t114 * t213 + t143 * t117;
t47 = t112 * t213 - t142 * t115 + t143 * t118;
t238 = (-t45 - t46 - t47) * t164 + ((t26 + t28 + t30) * t175 + (t25 + t27 + t29) * t173) * t163;
t31 = t164 * t63 + (t172 * t69 + t174 * t75) * t163;
t33 = -t164 * t71 + (t172 * t65 + t174 * t77) * t163;
t35 = -t164 * t67 + (-t172 * t73 + t174 * t79) * t163;
t237 = -t31 - t33 - t35;
t32 = t164 * t64 + (t172 * t70 + t174 * t76) * t163;
t34 = -t164 * t72 + (t172 * t66 + t174 * t78) * t163;
t36 = -t164 * t68 + (-t172 * t74 + t174 * t80) * t163;
t236 = t32 + t34 + t36;
t214 = t163 * t174;
t216 = t163 * t172;
t235 = t164 * t110 + (t111 + t113) * t216 + (t116 + t117 + t118) * t214;
t234 = t164 ^ 2;
t167 = t173 ^ 2;
t168 = t175 ^ 2;
t233 = -t164 / 0.2e1;
t232 = t173 / 0.2e1;
t149 = t163 * rSges(4,1) + t164 * rSges(4,2);
t230 = m(4) * t149;
t229 = m(7) * t163;
t228 = pkin(3) * t164;
t85 = t143 * rSges(6,1) + rSges(6,2) * t213 + t142 * rSges(6,3);
t98 = t143 * pkin(4) + t142 * qJ(5);
t227 = -t85 - t98;
t226 = t140 * rSges(7,2);
t225 = t140 * rSges(6,3);
t224 = rSges(3,3) + qJ(2);
t221 = t244 * t141 - t243 * t215 + t226;
t131 = (pkin(4) * t174 + qJ(5) * t172) * t163;
t130 = t140 * qJ(5);
t97 = t141 * pkin(4) + t130;
t219 = t131 * t215 + t164 * t97;
t217 = Icges(4,4) * t164;
t212 = t164 * t175;
t211 = t172 * t115;
t171 = -pkin(7) - qJ(2);
t208 = t175 * t171;
t205 = (rSges(7,1) * t174 + rSges(7,2) * t172) * t163 + pkin(5) * t214 + t243 * t164;
t120 = -t164 * rSges(6,2) + (rSges(6,1) * t174 + rSges(6,3) * t172) * t163;
t204 = -t120 - t131;
t121 = -t164 * rSges(5,3) + (rSges(5,1) * t174 - rSges(5,2) * t172) * t163;
t150 = t163 * pkin(3) - t164 * pkin(8);
t203 = -t121 - t150;
t200 = pkin(3) * t212 + pkin(8) * t213;
t202 = t167 * (pkin(8) * t163 + t228) + t175 * t200;
t199 = t167 + t168;
t197 = t163 * t211 + t242 * t164 - t235;
t196 = -t98 - t220;
t195 = -t131 - t205;
t194 = -t150 + t204;
t86 = t143 * rSges(5,1) - t142 * rSges(5,2) + rSges(5,3) * t213;
t170 = cos(pkin(9));
t160 = t170 * pkin(2) + pkin(1);
t193 = -t160 - t228;
t192 = -t130 - t208;
t191 = t175 * t160 - t173 * t171;
t190 = t173 * t97 + t175 * t98 + t202;
t189 = -t150 + t195;
t188 = rSges(4,1) * t164 - rSges(4,2) * t163;
t187 = -t141 * rSges(5,1) + t140 * rSges(5,2);
t185 = -Icges(4,2) * t163 + t217;
t184 = Icges(4,5) * t164 - Icges(4,6) * t163;
t181 = rSges(4,1) * t212 - rSges(4,2) * t213 + t173 * rSges(4,3);
t169 = sin(pkin(9));
t180 = rSges(3,1) * t170 - rSges(3,2) * t169 + pkin(1);
t179 = t191 + t200;
t178 = t34 / 0.2e1 + t47 / 0.2e1 + t46 / 0.2e1 + t45 / 0.2e1 + t32 / 0.2e1 + t36 / 0.2e1;
t177 = t35 / 0.2e1 + t33 / 0.2e1 + t44 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1 + t31 / 0.2e1;
t176 = t179 + t98;
t162 = t163 ^ 2;
t152 = t175 * rSges(2,1) - t173 * rSges(2,2);
t151 = -t173 * rSges(2,1) - t175 * rSges(2,2);
t146 = Icges(4,6) * t164 + t246;
t123 = Icges(4,3) * t173 + t175 * t184;
t122 = -Icges(4,3) * t175 + t173 * t184;
t108 = t173 * t224 + t175 * t180;
t107 = -t173 * t180 + t175 * t224;
t96 = t181 + t191;
t95 = (rSges(4,3) - t171) * t175 + (-t160 - t188) * t173;
t89 = t203 * t175;
t88 = t203 * t173;
t87 = t97 * t213;
t83 = rSges(5,3) * t215 - t187;
t82 = t141 * rSges(6,1) + rSges(6,2) * t215 + t225;
t62 = t175 * t181 + (-t175 * rSges(4,3) + t173 * t188) * t173;
t61 = t194 * t175;
t60 = t194 * t173;
t59 = t179 + t86;
t58 = -t208 + ((-rSges(5,3) - pkin(8)) * t163 + t193) * t173 + t187;
t57 = -t121 * t213 - t164 * t86;
t56 = t121 * t215 + t164 * t83;
t55 = t189 * t175;
t54 = t189 * t173;
t50 = (-t173 * t86 + t175 * t83) * t163;
t49 = t176 + t85;
t48 = -t225 + (-rSges(6,1) - pkin(4)) * t141 + ((-rSges(6,2) - pkin(8)) * t163 + t193) * t173 + t192;
t41 = t173 * t83 + t175 * t86 + t202;
t40 = t176 + t220;
t39 = -t226 + (-pkin(4) - t244) * t141 + ((-pkin(8) + t243) * t163 + t193) * t173 + t192;
t38 = t164 * t227 + t204 * t213;
t37 = t120 * t215 + t164 * t82 + t219;
t18 = t87 + (t173 * t227 + t175 * t82) * t163;
t17 = t164 * t196 + t195 * t213;
t16 = t164 * t221 + t205 * t215 + t219;
t15 = t173 * t82 + t175 * t85 + t190;
t14 = t87 + (t173 * t196 + t175 * t221) * t163;
t13 = t173 * t221 + t175 * t220 + t190;
t12 = t30 * t173 - t29 * t175;
t11 = t28 * t173 - t27 * t175;
t10 = t26 * t173 - t25 * t175;
t9 = t24 * t173 - t23 * t175;
t8 = t22 * t173 - t21 * t175;
t7 = t20 * t173 - t19 * t175;
t1 = [Icges(3,2) * t170 ^ 2 + Icges(2,3) + (Icges(3,1) * t169 + 0.2e1 * Icges(3,4) * t170) * t169 + (Icges(4,1) * t163 - t211 + t217) * t163 + (Icges(4,4) * t163 + Icges(4,2) * t164 - t242) * t164 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t95 ^ 2 + t96 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t151 ^ 2 + t152 ^ 2) + t235; m(7) * (t173 * t39 - t175 * t40) + m(6) * (t173 * t48 - t175 * t49) + m(5) * (t173 * t58 - t175 * t59) + m(4) * (t173 * t95 - t175 * t96) + m(3) * (t173 * t107 - t175 * t108); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t198) * t199; m(7) * (t55 * t39 + t54 * t40) + m(6) * (t61 * t48 + t60 * t49) + m(5) * (t89 * t58 + t88 * t59) + (t173 * t185 * t233 - t95 * t230 - t177 + (-Icges(4,6) * t233 + t245 + t146 / 0.2e1) * t175) * t175 + (t164 * (Icges(4,6) * t173 + t175 * t185) / 0.2e1 + t173 * t245 - t96 * t230 + t146 * t232 + t178) * t173; m(5) * (t89 * t173 - t88 * t175) + m(6) * (t61 * t173 - t60 * t175) + m(7) * (t55 * t173 - t54 * t175); m(7) * (t13 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t41 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t149 ^ 2 * t199 + t62 ^ 2) + (-t168 * t122 - t7 - t8 - t9) * t175 + (t167 * t123 + t10 + t11 + t12 + (-t173 * t122 + t175 * t123) * t175) * t173; t197 * t164 + m(7) * (t16 * t39 + t17 * t40) + m(6) * (t37 * t48 + t38 * t49) + m(5) * (t56 * t58 + t57 * t59) + (t173 * t177 + t175 * t178) * t163; m(5) * (t56 * t173 - t57 * t175) + m(6) * (t37 * t173 - t38 * t175) + m(7) * (t16 * t173 - t17 * t175); m(7) * (t14 * t13 + t16 * t55 + t17 * t54) + m(6) * (t18 * t15 + t37 * t61 + t38 * t60) + m(5) * (t50 * t41 + t56 * t89 + t57 * t88) + ((t10 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t175 + (t9 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t173) * t163 + (t236 * t173 + t237 * t175) * t233 + t238 * t232 - t239 * t175 / 0.2e1; m(6) * (t18 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(5) * (t50 ^ 2 + t56 ^ 2 + t57 ^ 2) - t197 * t234 + (t238 * t175 + t239 * t173 + (t237 * t173 - t236 * t175) * t164) * t163; m(7) * (t140 * t40 + t142 * t39) + m(6) * (t140 * t49 + t142 * t48); (-t140 * t175 + t142 * t173) * t241; m(7) * (t13 * t216 + t140 * t54 + t142 * t55) + m(6) * (t140 * t60 + t142 * t61 + t15 * t216); m(6) * (t140 * t38 + t142 * t37 + t18 * t216) + m(7) * (t14 * t216 + t140 * t17 + t142 * t16); (t162 * t172 ^ 2 + t140 ^ 2 + t142 ^ 2) * t241; (-t173 * t40 - t175 * t39) * t229; 0; m(7) * (t164 * t13 + (-t173 * t54 - t175 * t55) * t163); m(7) * (t164 * t14 + (-t16 * t175 - t17 * t173) * t163); (-t140 * t173 - t142 * t175 + t164 * t172) * t229; m(7) * (t162 * t199 + t234);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
