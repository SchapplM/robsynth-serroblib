% Calculate joint inertia matrix for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:47
% EndTime: 2019-03-09 04:38:54
% DurationCPUTime: 3.05s
% Computational Cost: add. (7808->406), mult. (7722->597), div. (0->0), fcn. (8249->10), ass. (0->198)
t173 = pkin(9) + qJ(3);
t168 = sin(t173);
t258 = Icges(4,5) * t168;
t257 = t258 / 0.2e1;
t170 = cos(t173);
t174 = qJ(4) + pkin(10);
t169 = sin(t174);
t182 = sin(qJ(1));
t221 = t182 * t169;
t171 = cos(t174);
t184 = cos(qJ(1));
t223 = t171 * t184;
t134 = t170 * t221 + t223;
t220 = t182 * t171;
t135 = -t169 * t184 + t170 * t220;
t254 = rSges(7,3) + qJ(6);
t255 = rSges(7,1) + pkin(5);
t256 = -t254 * t134 - t255 * t135;
t106 = -Icges(6,3) * t170 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t168;
t107 = -Icges(7,2) * t170 + (Icges(7,4) * t171 + Icges(7,6) * t169) * t168;
t181 = sin(qJ(4));
t183 = cos(qJ(4));
t115 = -Icges(5,3) * t170 + (Icges(5,5) * t183 - Icges(5,6) * t181) * t168;
t253 = -t106 - t107 - t115;
t226 = t168 * t182;
t63 = Icges(7,5) * t135 + Icges(7,6) * t226 + Icges(7,3) * t134;
t67 = Icges(7,4) * t135 + Icges(7,2) * t226 + Icges(7,6) * t134;
t71 = Icges(7,1) * t135 + Icges(7,4) * t226 + Icges(7,5) * t134;
t19 = t134 * t63 + t135 * t71 + t67 * t226;
t224 = t170 * t184;
t136 = t169 * t224 - t220;
t137 = t170 * t223 + t221;
t225 = t168 * t184;
t64 = Icges(7,5) * t137 + Icges(7,6) * t225 + Icges(7,3) * t136;
t68 = Icges(7,4) * t137 + Icges(7,2) * t225 + Icges(7,6) * t136;
t72 = Icges(7,1) * t137 + Icges(7,4) * t225 + Icges(7,5) * t136;
t20 = t134 * t64 + t135 * t72 + t68 * t226;
t65 = Icges(6,5) * t135 - Icges(6,6) * t134 + Icges(6,3) * t226;
t69 = Icges(6,4) * t135 - Icges(6,2) * t134 + Icges(6,6) * t226;
t73 = Icges(6,1) * t135 - Icges(6,4) * t134 + Icges(6,5) * t226;
t21 = -t134 * t69 + t135 * t73 + t65 * t226;
t66 = Icges(6,5) * t137 - Icges(6,6) * t136 + Icges(6,3) * t225;
t70 = Icges(6,4) * t137 - Icges(6,2) * t136 + Icges(6,6) * t225;
t74 = Icges(6,1) * t137 - Icges(6,4) * t136 + Icges(6,5) * t225;
t22 = -t134 * t70 + t135 * t74 + t66 * t226;
t217 = t183 * t184;
t219 = t182 * t181;
t140 = -t170 * t219 - t217;
t218 = t182 * t183;
t222 = t181 * t184;
t141 = t170 * t218 - t222;
t85 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t226;
t87 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t226;
t89 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t226;
t33 = t140 * t87 + t141 * t89 + t85 * t226;
t142 = -t170 * t222 + t218;
t143 = t170 * t217 + t219;
t86 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t225;
t88 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t225;
t90 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t225;
t34 = t140 * t88 + t141 * t90 + t86 * t226;
t105 = -Icges(7,6) * t170 + (Icges(7,5) * t171 + Icges(7,3) * t169) * t168;
t109 = -Icges(7,4) * t170 + (Icges(7,1) * t171 + Icges(7,5) * t169) * t168;
t42 = t105 * t134 + t107 * t226 + t109 * t135;
t108 = -Icges(6,6) * t170 + (Icges(6,4) * t171 - Icges(6,2) * t169) * t168;
t110 = -Icges(6,5) * t170 + (Icges(6,1) * t171 - Icges(6,4) * t169) * t168;
t43 = t106 * t226 - t108 * t134 + t110 * t135;
t116 = -Icges(5,6) * t170 + (Icges(5,4) * t183 - Icges(5,2) * t181) * t168;
t117 = -Icges(5,5) * t170 + (Icges(5,1) * t183 - Icges(5,4) * t181) * t168;
t46 = t115 * t226 + t116 * t140 + t117 * t141;
t252 = (-t42 - t43 - t46) * t170 + ((t20 + t22 + t34) * t184 + (t19 + t21 + t33) * t182) * t168;
t23 = t136 * t63 + t137 * t71 + t67 * t225;
t24 = t136 * t64 + t137 * t72 + t68 * t225;
t25 = -t136 * t69 + t137 * t73 + t65 * t225;
t26 = -t136 * t70 + t137 * t74 + t66 * t225;
t35 = t142 * t87 + t143 * t89 + t85 * t225;
t36 = t142 * t88 + t143 * t90 + t86 * t225;
t44 = t136 * t105 + t107 * t225 + t137 * t109;
t45 = t106 * t225 - t136 * t108 + t137 * t110;
t47 = t115 * t225 + t142 * t116 + t143 * t117;
t251 = (-t44 - t45 - t47) * t170 + ((t24 + t26 + t36) * t184 + (t23 + t25 + t35) * t182) * t168;
t27 = -t170 * t67 + (t169 * t63 + t171 * t71) * t168;
t29 = -t170 * t65 + (-t169 * t69 + t171 * t73) * t168;
t37 = -t170 * t85 + (-t181 * t87 + t183 * t89) * t168;
t250 = -t27 - t29 - t37;
t28 = -t170 * t68 + (t169 * t64 + t171 * t72) * t168;
t30 = -t170 * t66 + (-t169 * t70 + t171 * t74) * t168;
t38 = -t170 * t86 + (-t181 * t88 + t183 * t90) * t168;
t249 = t28 + t30 + t38;
t228 = t168 * t169;
t248 = t105 * t228 + (t183 * t117 + (t109 + t110) * t171) * t168;
t247 = t170 ^ 2;
t175 = t182 ^ 2;
t176 = t184 ^ 2;
t246 = m(6) / 0.2e1;
t245 = m(7) / 0.2e1;
t244 = -t170 / 0.2e1;
t243 = t182 / 0.2e1;
t148 = rSges(4,1) * t168 + rSges(4,2) * t170;
t241 = m(4) * t148;
t240 = pkin(3) * t170;
t239 = pkin(8) * t168;
t166 = pkin(4) * t183 + pkin(3);
t238 = -pkin(3) + t166;
t237 = rSges(7,2) * t226 - t256;
t236 = rSges(7,2) * t225 + t254 * t136 + t255 * t137;
t78 = t137 * rSges(6,1) - t136 * rSges(6,2) + rSges(6,3) * t225;
t179 = -qJ(5) - pkin(8);
t188 = pkin(4) * t219 + t166 * t224 - t179 * t225;
t210 = pkin(3) * t224 + pkin(8) * t225;
t84 = t188 - t210;
t235 = -t78 - t84;
t104 = (pkin(8) + t179) * t170 + t238 * t168;
t211 = -pkin(4) * t222 - t179 * t226;
t83 = (t238 * t170 - t239) * t182 + t211;
t234 = t104 * t226 + t170 * t83;
t232 = rSges(3,3) + qJ(2);
t230 = Icges(4,4) * t170;
t229 = t116 * t181;
t180 = -pkin(7) - qJ(2);
t216 = t184 * t180;
t112 = -t170 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t169) * t168;
t215 = -t104 - t112;
t214 = -t170 * rSges(7,2) + (t254 * t169 + t255 * t171) * t168;
t118 = -t170 * rSges(5,3) + (rSges(5,1) * t183 - rSges(5,2) * t181) * t168;
t150 = t168 * pkin(3) - t170 * pkin(8);
t213 = -t118 - t150;
t212 = t175 * (t239 + t240) + t184 * t210;
t209 = t175 + t176;
t208 = t246 + t245;
t207 = -t108 * t228 - t168 * t229 + t253 * t170 + t248;
t206 = -t84 - t236;
t205 = -t104 - t214;
t204 = -t150 + t215;
t92 = t143 * rSges(5,1) + t142 * rSges(5,2) + rSges(5,3) * t225;
t178 = cos(pkin(9));
t165 = pkin(2) * t178 + pkin(1);
t203 = t184 * t165 - t182 * t180;
t202 = -t166 * t170 - t165;
t201 = t182 * t83 + t184 * t84 + t212;
t200 = -t150 + t205;
t199 = -t211 - t216;
t198 = rSges(4,1) * t170 - rSges(4,2) * t168;
t197 = -t141 * rSges(5,1) - t140 * rSges(5,2);
t196 = -t135 * rSges(6,1) + t134 * rSges(6,2);
t194 = -Icges(4,2) * t168 + t230;
t193 = Icges(4,5) * t170 - Icges(4,6) * t168;
t190 = rSges(4,1) * t224 - rSges(4,2) * t225 + t182 * rSges(4,3);
t177 = sin(pkin(9));
t189 = rSges(3,1) * t178 - rSges(3,2) * t177 + pkin(1);
t187 = t30 / 0.2e1 + t47 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t28 / 0.2e1 + t38 / 0.2e1;
t186 = t37 / 0.2e1 + t29 / 0.2e1 + t46 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1 + t27 / 0.2e1;
t185 = t188 + t203;
t167 = t168 ^ 2;
t153 = rSges(2,1) * t184 - t182 * rSges(2,2);
t152 = -t182 * rSges(2,1) - rSges(2,2) * t184;
t145 = Icges(4,6) * t170 + t258;
t120 = Icges(4,3) * t182 + t193 * t184;
t119 = -Icges(4,3) * t184 + t193 * t182;
t114 = t232 * t182 + t189 * t184;
t113 = -t189 * t182 + t232 * t184;
t102 = t190 + t203;
t101 = (rSges(4,3) - t180) * t184 + (-t165 - t198) * t182;
t96 = t213 * t184;
t95 = t213 * t182;
t91 = rSges(5,3) * t226 - t197;
t80 = t184 * t190 + (-t184 * rSges(4,3) + t198 * t182) * t182;
t76 = rSges(6,3) * t226 - t196;
t62 = t83 * t225;
t61 = t203 + t92 + t210;
t60 = -t216 + (-t240 - t165 + (-rSges(5,3) - pkin(8)) * t168) * t182 + t197;
t59 = t204 * t184;
t58 = t204 * t182;
t57 = -t118 * t225 - t170 * t92;
t56 = t118 * t226 + t170 * t91;
t54 = t185 + t78;
t53 = (-rSges(6,3) * t168 + t202) * t182 + t196 + t199;
t52 = (-t182 * t92 + t184 * t91) * t168;
t51 = t200 * t184;
t50 = t200 * t182;
t41 = t182 * t91 + t184 * t92 + t212;
t40 = t185 + t236;
t39 = (-rSges(7,2) * t168 + t202) * t182 + t199 + t256;
t32 = t235 * t170 + t215 * t225;
t31 = t112 * t226 + t170 * t76 + t234;
t18 = t62 + (t235 * t182 + t184 * t76) * t168;
t17 = t182 * t76 + t184 * t78 + t201;
t16 = t206 * t170 + t205 * t225;
t15 = t237 * t170 + t214 * t226 + t234;
t14 = t62 + (t206 * t182 + t237 * t184) * t168;
t13 = t237 * t182 + t236 * t184 + t201;
t12 = t36 * t182 - t184 * t35;
t11 = t34 * t182 - t184 * t33;
t10 = t26 * t182 - t184 * t25;
t9 = t24 * t182 - t184 * t23;
t8 = t22 * t182 - t184 * t21;
t7 = t20 * t182 - t184 * t19;
t1 = [Icges(3,2) * t178 ^ 2 + Icges(2,3) + (Icges(3,1) * t177 + 0.2e1 * Icges(3,4) * t178) * t177 + (Icges(4,1) * t168 - t108 * t169 - t229 + t230) * t168 + (Icges(4,4) * t168 + Icges(4,2) * t170 + t253) * t170 + m(6) * (t53 ^ 2 + t54 ^ 2) + m(7) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t152 ^ 2 + t153 ^ 2) + t248; m(6) * (t182 * t53 - t184 * t54) + m(7) * (t182 * t39 - t184 * t40) + m(5) * (t182 * t60 - t184 * t61) + m(4) * (t182 * t101 - t102 * t184) + m(3) * (t182 * t113 - t114 * t184); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t208) * t209; m(6) * (t53 * t59 + t54 * t58) + m(7) * (t39 * t51 + t40 * t50) + m(5) * (t60 * t96 + t61 * t95) + (t194 * t182 * t244 - t101 * t241 - t186 + (-Icges(4,6) * t244 + t257 + t145 / 0.2e1) * t184) * t184 + (t170 * (Icges(4,6) * t182 + t194 * t184) / 0.2e1 + t182 * t257 - t102 * t241 + t145 * t243 + t187) * t182; m(5) * (t96 * t182 - t184 * t95) + m(6) * (t59 * t182 - t184 * t58) + m(7) * (t51 * t182 - t184 * t50); m(7) * (t13 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t17 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t41 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t209 * t148 ^ 2 + t80 ^ 2) + (-t176 * t119 - t11 - t7 - t8) * t184 + (t175 * t120 + t10 + t12 + t9 + (-t182 * t119 + t184 * t120) * t184) * t182; -t207 * t170 + m(6) * (t31 * t53 + t32 * t54) + m(7) * (t15 * t39 + t16 * t40) + m(5) * (t56 * t60 + t57 * t61) + (t186 * t182 + t187 * t184) * t168; m(5) * (t56 * t182 - t184 * t57) + m(6) * (t31 * t182 - t184 * t32) + m(7) * (t15 * t182 - t16 * t184); m(7) * (t13 * t14 + t15 * t51 + t16 * t50) + m(6) * (t17 * t18 + t31 * t59 + t32 * t58) + m(5) * (t41 * t52 + t56 * t96 + t57 * t95) + ((t9 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t184 + (t11 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t182) * t168 + (t182 * t249 + t184 * t250) * t244 + t251 * t243 - t252 * t184 / 0.2e1; m(7) * (t14 ^ 2 + t15 ^ 2 + t16 ^ 2) + m(6) * (t18 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t52 ^ 2 + t56 ^ 2 + t57 ^ 2) + t207 * t247 + (t251 * t184 + t252 * t182 + (t182 * t250 - t184 * t249) * t170) * t168; 0.2e1 * ((t182 * t54 + t184 * t53) * t246 + (t182 * t40 + t184 * t39) * t245) * t168; 0; m(7) * (-t170 * t13 + (t182 * t50 + t184 * t51) * t168) + m(6) * (-t170 * t17 + (t182 * t58 + t184 * t59) * t168); m(7) * (-t170 * t14 + (t15 * t184 + t16 * t182) * t168) + m(6) * (-t170 * t18 + (t182 * t32 + t184 * t31) * t168); 0.2e1 * t208 * (t209 * t167 + t247); m(7) * (t134 * t40 + t136 * t39); m(7) * (-t134 * t184 + t136 * t182); m(7) * (t13 * t228 + t134 * t50 + t136 * t51); m(7) * (t134 * t16 + t136 * t15 + t14 * t228); m(7) * (t134 * t182 + t136 * t184 - t169 * t170) * t168; m(7) * (t167 * t169 ^ 2 + t134 ^ 2 + t136 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
