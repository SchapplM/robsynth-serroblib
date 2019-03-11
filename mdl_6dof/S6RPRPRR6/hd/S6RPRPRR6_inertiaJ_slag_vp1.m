% Calculate joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:15
% EndTime: 2019-03-09 03:51:20
% DurationCPUTime: 2.88s
% Computational Cost: add. (9667->426), mult. (7747->629), div. (0->0), fcn. (8241->12), ass. (0->202)
t187 = sin(qJ(1));
t250 = t187 / 0.2e1;
t188 = cos(qJ(1));
t248 = t188 / 0.2e1;
t183 = cos(pkin(11));
t166 = t183 * pkin(4) + pkin(3);
t177 = pkin(11) + qJ(5);
t171 = cos(t177);
t149 = pkin(5) * t171 + t166;
t169 = sin(t177);
t181 = sin(pkin(11));
t150 = pkin(4) * t181 + pkin(5) * t169;
t178 = pkin(10) + qJ(3);
t172 = cos(t178);
t232 = t172 * t188;
t173 = qJ(6) + t177;
t164 = sin(t173);
t165 = cos(t173);
t227 = t187 * t165;
t117 = -t164 * t232 + t227;
t228 = t187 * t164;
t118 = t165 * t232 + t228;
t170 = sin(t178);
t233 = t170 * t188;
t72 = t118 * rSges(7,1) + t117 * rSges(7,2) + rSges(7,3) * t233;
t256 = t149 * t232 + t187 * t150 + t72;
t179 = t187 ^ 2;
t180 = t188 ^ 2;
t255 = m(5) / 0.2e1;
t254 = m(6) / 0.2e1;
t253 = m(7) / 0.2e1;
t234 = t170 * t187;
t115 = -t165 * t188 - t172 * t228;
t116 = -t164 * t188 + t172 * t227;
t65 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t234;
t67 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t234;
t69 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t234;
t22 = t115 * t67 + t116 * t69 + t234 * t65;
t66 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t233;
t68 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t233;
t70 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t233;
t23 = t115 * t68 + t116 * t70 + t234 * t66;
t100 = -Icges(7,6) * t172 + (Icges(7,4) * t165 - Icges(7,2) * t164) * t170;
t101 = -Icges(7,5) * t172 + (Icges(7,1) * t165 - Icges(7,4) * t164) * t170;
t99 = -Icges(7,3) * t172 + (Icges(7,5) * t165 - Icges(7,6) * t164) * t170;
t38 = t100 * t115 + t101 * t116 + t234 * t99;
t5 = -t38 * t172 + (t187 * t22 + t188 * t23) * t170;
t24 = t117 * t67 + t118 * t69 + t233 * t65;
t25 = t117 * t68 + t118 * t70 + t233 * t66;
t39 = t117 * t100 + t118 * t101 + t233 * t99;
t6 = -t39 * t172 + (t187 * t24 + t188 * t25) * t170;
t252 = t6 * t233 + t5 * t234;
t251 = -t172 / 0.2e1;
t249 = -t188 / 0.2e1;
t148 = rSges(4,1) * t170 + rSges(4,2) * t172;
t247 = m(4) * t148;
t246 = pkin(3) * t172;
t245 = -pkin(3) + t166;
t185 = -pkin(8) - qJ(4);
t176 = -pkin(9) + t185;
t215 = t176 - t185;
t224 = t187 * t181;
t219 = -pkin(4) * t224 - t166 * t232;
t244 = -t215 * t233 + t219 + t256;
t102 = -t172 * rSges(7,3) + (rSges(7,1) * t165 - rSges(7,2) * t164) * t170;
t197 = -t116 * rSges(7,1) - t115 * rSges(7,2);
t71 = rSges(7,3) * t234 - t197;
t51 = t102 * t234 + t172 * t71;
t236 = t100 * t164;
t93 = t170 * t165 * t101;
t47 = -t170 * t236 - t172 * t99 + t93;
t243 = t47 * t172;
t242 = rSges(3,3) + qJ(2);
t218 = t149 - t166;
t90 = t170 * t218 + t172 * t215;
t241 = -t102 - t90;
t147 = t170 * pkin(3) - t172 * qJ(4);
t240 = -t147 - (qJ(4) + t185) * t172 - t245 * t170;
t239 = Icges(4,4) * t170;
t238 = Icges(4,4) * t172;
t237 = qJ(4) * t170;
t104 = -Icges(6,6) * t172 + (Icges(6,4) * t171 - Icges(6,2) * t169) * t170;
t235 = t104 * t169;
t231 = t181 * t188;
t230 = t183 * t188;
t186 = -pkin(7) - qJ(2);
t229 = t186 * t188;
t226 = t187 * t169;
t225 = t187 * t171;
t223 = t187 * t183;
t222 = t172 * rSges(5,3) - (rSges(5,1) * t183 - rSges(5,2) * t181) * t170 - t147;
t216 = pkin(3) * t232 + qJ(4) * t233;
t221 = t179 * (t237 + t246) + t188 * t216;
t217 = pkin(4) * t231 + t185 * t234;
t214 = t179 + t180;
t130 = -t171 * t188 - t172 * t226;
t131 = -t169 * t188 + t172 * t225;
t73 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t234;
t75 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t234;
t77 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t234;
t34 = -t172 * t73 + (-t169 * t75 + t171 * t77) * t170;
t103 = -Icges(6,3) * t172 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t170;
t105 = -Icges(6,5) * t172 + (Icges(6,1) * t171 - Icges(6,4) * t169) * t170;
t41 = t103 * t234 + t104 * t130 + t105 * t131;
t213 = t34 / 0.2e1 + t41 / 0.2e1;
t132 = -t169 * t232 + t225;
t133 = t171 * t232 + t226;
t74 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t233;
t76 = Icges(6,4) * t133 + Icges(6,2) * t132 + Icges(6,6) * t233;
t78 = Icges(6,1) * t133 + Icges(6,4) * t132 + Icges(6,5) * t233;
t35 = -t172 * t74 + (-t169 * t76 + t171 * t78) * t170;
t42 = t103 * t233 + t132 * t104 + t133 * t105;
t212 = t35 / 0.2e1 + t42 / 0.2e1;
t106 = -t172 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t169) * t170;
t211 = -t106 + t240;
t80 = t133 * rSges(6,1) + t132 * rSges(6,2) + rSges(6,3) * t233;
t139 = -t172 * t231 + t223;
t140 = t172 * t230 + t224;
t210 = t140 * rSges(5,1) + t139 * rSges(5,2) + rSges(5,3) * t233;
t209 = t234 / 0.2e1;
t208 = t233 / 0.2e1;
t32 = -t172 * t65 + (-t164 * t67 + t165 * t69) * t170;
t33 = -t172 * t66 + (-t164 * t68 + t165 * t70) * t170;
t207 = (t32 + t38) * t209 + (t33 + t39) * t208;
t184 = cos(pkin(10));
t167 = pkin(2) * t184 + pkin(1);
t206 = t188 * t167 - t187 * t186;
t205 = t240 + t241;
t189 = -t185 * t233 - t219;
t204 = t187 * ((t172 * t245 - t237) * t187 - t217) + t188 * (t189 - t216) + t221;
t9 = -t243 + (t187 * t32 + t188 * t33) * t170;
t203 = -t172 * t9 + t252;
t202 = t255 + t254 + t253;
t12 = t23 * t187 - t188 * t22;
t13 = t25 * t187 - t188 * t24;
t201 = t12 * t209 + t13 * t208 + t5 * t249 + t6 * t250 + (t33 * t187 - t32 * t188) * t251;
t200 = rSges(4,1) * t172 - rSges(4,2) * t170;
t137 = -t172 * t224 - t230;
t138 = t172 * t223 - t231;
t199 = -t138 * rSges(5,1) - t137 * rSges(5,2);
t198 = -t131 * rSges(6,1) - t130 * rSges(6,2);
t196 = Icges(4,1) * t172 - t239;
t195 = -Icges(4,2) * t170 + t238;
t194 = Icges(4,5) * t172 - Icges(4,6) * t170;
t191 = rSges(4,1) * t232 - rSges(4,2) * t233 + t187 * rSges(4,3);
t182 = sin(pkin(10));
t190 = rSges(3,1) * t184 - rSges(3,2) * t182 + pkin(1);
t153 = rSges(2,1) * t188 - t187 * rSges(2,2);
t152 = -t187 * rSges(2,1) - rSges(2,2) * t188;
t143 = Icges(4,5) * t170 + Icges(4,6) * t172;
t120 = Icges(4,3) * t187 + t188 * t194;
t119 = -Icges(4,3) * t188 + t187 * t194;
t111 = -Icges(5,5) * t172 + (Icges(5,1) * t183 - Icges(5,4) * t181) * t170;
t110 = -Icges(5,6) * t172 + (Icges(5,4) * t183 - Icges(5,2) * t181) * t170;
t108 = t187 * t242 + t188 * t190;
t107 = -t187 * t190 + t188 * t242;
t97 = t191 + t206;
t96 = (rSges(4,3) - t186) * t188 + (-t167 - t200) * t187;
t94 = t170 * t171 * t105;
t92 = t222 * t188;
t91 = t222 * t187;
t89 = Icges(5,1) * t140 + Icges(5,4) * t139 + Icges(5,5) * t233;
t88 = Icges(5,1) * t138 + Icges(5,4) * t137 + Icges(5,5) * t234;
t87 = Icges(5,4) * t140 + Icges(5,2) * t139 + Icges(5,6) * t233;
t86 = Icges(5,4) * t138 + Icges(5,2) * t137 + Icges(5,6) * t234;
t85 = Icges(5,5) * t140 + Icges(5,6) * t139 + Icges(5,3) * t233;
t84 = Icges(5,5) * t138 + Icges(5,6) * t137 + Icges(5,3) * t234;
t83 = t188 * t191 + (-t188 * rSges(4,3) + t187 * t200) * t187;
t79 = rSges(6,3) * t234 - t198;
t63 = t71 * t233;
t61 = -t150 * t188 + (-t170 * t176 + t172 * t218) * t187 + t217;
t60 = t206 + t210 + t216;
t59 = -t229 + (-t246 - t167 + (-rSges(5,3) - qJ(4)) * t170) * t187 + t199;
t58 = t211 * t188;
t57 = t211 * t187;
t56 = -t106 * t233 - t172 * t80;
t55 = t106 * t234 + t172 * t79;
t54 = t189 + t206 + t80;
t53 = -t229 + (-rSges(6,3) * t170 - t166 * t172 - t167) * t187 + t198 + t217;
t52 = -t102 * t233 - t172 * t72;
t50 = -t172 * t103 - t170 * t235 + t94;
t49 = -t176 * t233 + t206 + t256;
t48 = (t150 - t186) * t188 + (-t149 * t172 - t167 + (-rSges(7,3) + t176) * t170) * t187 + t197;
t46 = (-t187 * t80 + t188 * t79) * t170;
t45 = -t234 * t72 + t63;
t44 = t205 * t188;
t43 = t205 * t187;
t40 = t187 * (rSges(5,3) * t234 - t199) + t188 * t210 + t221;
t29 = t132 * t76 + t133 * t78 + t233 * t74;
t28 = t132 * t75 + t133 * t77 + t233 * t73;
t27 = t130 * t76 + t131 * t78 + t234 * t74;
t26 = t130 * t75 + t131 * t77 + t234 * t73;
t21 = -t172 * t244 + t233 * t241;
t20 = t172 * t61 + t234 * t90 + t51;
t19 = t187 * t79 + t188 * t80 + t204;
t18 = t63 + (-t187 * t244 + t188 * t61) * t170;
t16 = t244 * t188 + (t61 + t71) * t187 + t204;
t15 = t29 * t187 - t188 * t28;
t14 = t27 * t187 - t188 * t26;
t8 = -t42 * t172 + (t187 * t28 + t188 * t29) * t170;
t7 = -t41 * t172 + (t187 * t26 + t188 * t27) * t170;
t1 = [Icges(3,2) * t184 ^ 2 + Icges(2,3) + t93 + t94 + (Icges(3,1) * t182 + 0.2e1 * Icges(3,4) * t184) * t182 + (-t99 - t103 + t239 - (Icges(5,5) * t183 - Icges(5,6) * t181) * t170 + (Icges(4,2) + Icges(5,3)) * t172) * t172 + (Icges(4,1) * t170 - t110 * t181 + t111 * t183 - t235 - t236 + t238) * t170 + m(7) * (t48 ^ 2 + t49 ^ 2) + m(6) * (t53 ^ 2 + t54 ^ 2) + m(5) * (t59 ^ 2 + t60 ^ 2) + m(4) * (t96 ^ 2 + t97 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t152 ^ 2 + t153 ^ 2); m(7) * (t187 * t48 - t188 * t49) + m(6) * (t187 * t53 - t188 * t54) + m(5) * (t187 * t59 - t188 * t60) + m(4) * (t187 * t96 - t188 * t97) + m(3) * (t187 * t107 - t108 * t188); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t202) * t214; m(7) * (t43 * t49 + t44 * t48) + m(6) * (t53 * t58 + t54 * t57) + m(5) * (t59 * t92 + t60 * t91) + (-t38 / 0.2e1 - t32 / 0.2e1 - t96 * t247 - t137 * t110 / 0.2e1 - t138 * t111 / 0.2e1 + t143 * t248 + (Icges(4,6) * t248 - t187 * t195 / 0.2e1 + t84 / 0.2e1) * t172 - t213) * t188 + (t33 / 0.2e1 + t39 / 0.2e1 - t97 * t247 + t139 * t110 / 0.2e1 + t140 * t111 / 0.2e1 + t143 * t250 + (Icges(4,6) * t250 + t195 * t248 - t85 / 0.2e1) * t172 + t212) * t187 + ((Icges(4,5) * t187 - t181 * t87 + t183 * t89 + t188 * t196) * t250 + (-Icges(4,5) * t188 - t181 * t86 + t183 * t88 + t187 * t196) * t249) * t170; m(5) * (t92 * t187 - t188 * t91) + m(6) * (t58 * t187 - t188 * t57) + m(7) * (t44 * t187 - t188 * t43); m(7) * (t16 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t19 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t40 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t148 ^ 2 * t214 + t83 ^ 2) + (-t180 * t119 - t12 - t14 + (t137 * t86 + t138 * t88 + t84 * t234) * t188) * t188 + (t13 + t15 + t179 * t120 + (t139 * t87 + t140 * t89 + t85 * t233) * t187 + (-t187 * t119 + t188 * t120 - t137 * t87 - t138 * t89 - t139 * t86 - t140 * t88 - t233 * t84 - t234 * t85) * t188) * t187; 0.2e1 * ((t187 * t49 + t188 * t48) * t253 + (t187 * t54 + t188 * t53) * t254 + (t187 * t60 + t188 * t59) * t255) * t170; 0; m(7) * (-t172 * t16 + (t187 * t43 + t188 * t44) * t170) + m(6) * (-t172 * t19 + (t187 * t57 + t188 * t58) * t170) + m(5) * (-t172 * t40 + (t187 * t91 + t188 * t92) * t170); 0.2e1 * t202 * (t170 ^ 2 * t214 + t172 ^ 2); (-t47 - t50) * t172 + m(7) * (t20 * t48 + t21 * t49) + m(6) * (t53 * t55 + t54 * t56) + (t187 * t213 + t188 * t212) * t170 + t207; m(6) * (t55 * t187 - t188 * t56) + m(7) * (t20 * t187 - t188 * t21); t8 * t250 + t7 * t249 + (t35 * t187 - t34 * t188) * t251 + (t14 * t250 + t15 * t248) * t170 + m(7) * (t16 * t18 + t20 * t44 + t21 * t43) + m(6) * (t19 * t46 + t55 * t58 + t56 * t57) + t201; m(6) * (-t46 * t172 + (t187 * t56 + t188 * t55) * t170) + m(7) * (-t18 * t172 + (t187 * t21 + t188 * t20) * t170); (t50 * t172 - t9) * t172 + m(7) * (t18 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t46 ^ 2 + t55 ^ 2 + t56 ^ 2) + (t188 * t8 + t187 * t7 - t172 * (t187 * t34 + t188 * t35)) * t170 + t252; -t243 + m(7) * (t48 * t51 + t49 * t52) + t207; m(7) * (t51 * t187 - t188 * t52); m(7) * (t16 * t45 + t43 * t52 + t44 * t51) + t201; m(7) * (-t45 * t172 + (t187 * t52 + t188 * t51) * t170); m(7) * (t18 * t45 + t20 * t51 + t21 * t52) + t203; m(7) * (t45 ^ 2 + t51 ^ 2 + t52 ^ 2) + t203;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
