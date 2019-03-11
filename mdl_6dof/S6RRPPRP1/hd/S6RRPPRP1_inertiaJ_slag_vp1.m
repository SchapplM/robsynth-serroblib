% Calculate joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:27
% EndTime: 2019-03-09 08:25:35
% DurationCPUTime: 3.37s
% Computational Cost: add. (6504->420), mult. (6361->624), div. (0->0), fcn. (6725->10), ass. (0->192)
t268 = Icges(3,3) + Icges(4,3);
t168 = qJ(2) + pkin(9);
t161 = sin(t168);
t163 = cos(t168);
t175 = sin(qJ(2));
t177 = cos(qJ(2));
t267 = Icges(3,5) * t177 + Icges(4,5) * t163 - Icges(3,6) * t175 - Icges(4,6) * t161;
t167 = pkin(10) + qJ(5);
t162 = cos(t167);
t178 = cos(qJ(1));
t160 = sin(t167);
t176 = sin(qJ(1));
t220 = t176 * t160;
t115 = t162 * t178 + t163 * t220;
t219 = t176 * t162;
t116 = -t160 * t178 + t163 * t219;
t262 = rSges(7,3) + qJ(6);
t264 = rSges(7,1) + pkin(5);
t266 = -t262 * t115 - t264 * t116;
t247 = t176 / 0.2e1;
t265 = t178 / 0.2e1;
t83 = -Icges(6,3) * t163 + (Icges(6,5) * t162 - Icges(6,6) * t160) * t161;
t84 = -Icges(7,2) * t163 + (Icges(7,4) * t162 + Icges(7,6) * t160) * t161;
t263 = -t83 - t84;
t225 = t161 * t176;
t47 = Icges(7,5) * t116 + Icges(7,6) * t225 + Icges(7,3) * t115;
t51 = Icges(7,4) * t116 + Icges(7,2) * t225 + Icges(7,6) * t115;
t55 = Icges(7,1) * t116 + Icges(7,4) * t225 + Icges(7,5) * t115;
t12 = t115 * t47 + t116 * t55 + t225 * t51;
t223 = t163 * t178;
t117 = t160 * t223 - t219;
t118 = t162 * t223 + t220;
t224 = t161 * t178;
t48 = Icges(7,5) * t118 + Icges(7,6) * t224 + Icges(7,3) * t117;
t52 = Icges(7,4) * t118 + Icges(7,2) * t224 + Icges(7,6) * t117;
t56 = Icges(7,1) * t118 + Icges(7,4) * t224 + Icges(7,5) * t117;
t13 = t115 * t48 + t116 * t56 + t225 * t52;
t49 = Icges(6,5) * t116 - Icges(6,6) * t115 + Icges(6,3) * t225;
t53 = Icges(6,4) * t116 - Icges(6,2) * t115 + Icges(6,6) * t225;
t57 = Icges(6,1) * t116 - Icges(6,4) * t115 + Icges(6,5) * t225;
t14 = -t115 * t53 + t116 * t57 + t225 * t49;
t50 = Icges(6,5) * t118 - Icges(6,6) * t117 + Icges(6,3) * t224;
t54 = Icges(6,4) * t118 - Icges(6,2) * t117 + Icges(6,6) * t224;
t58 = Icges(6,1) * t118 - Icges(6,4) * t117 + Icges(6,5) * t224;
t15 = -t115 * t54 + t116 * t58 + t225 * t50;
t82 = -Icges(7,6) * t163 + (Icges(7,5) * t162 + Icges(7,3) * t160) * t161;
t86 = -Icges(7,4) * t163 + (Icges(7,1) * t162 + Icges(7,5) * t160) * t161;
t29 = t115 * t82 + t116 * t86 + t225 * t84;
t85 = -Icges(6,6) * t163 + (Icges(6,4) * t162 - Icges(6,2) * t160) * t161;
t87 = -Icges(6,5) * t163 + (Icges(6,1) * t162 - Icges(6,4) * t160) * t161;
t30 = -t115 * t85 + t116 * t87 + t225 * t83;
t261 = (-t29 - t30) * t163 + ((t13 + t15) * t178 + (t12 + t14) * t176) * t161;
t16 = t117 * t47 + t118 * t55 + t224 * t51;
t17 = t117 * t48 + t118 * t56 + t224 * t52;
t18 = -t117 * t53 + t118 * t57 + t224 * t49;
t19 = -t117 * t54 + t118 * t58 + t224 * t50;
t31 = t117 * t82 + t118 * t86 + t224 * t84;
t32 = -t117 * t85 + t118 * t87 + t224 * t83;
t260 = (-t31 - t32) * t163 + ((t17 + t19) * t178 + (t16 + t18) * t176) * t161;
t259 = t175 / 0.2e1;
t258 = t177 / 0.2e1;
t21 = -t163 * t51 + (t160 * t47 + t162 * t55) * t161;
t23 = -t163 * t49 + (-t160 * t53 + t162 * t57) * t161;
t257 = -t21 - t23;
t22 = -t163 * t52 + (t160 * t48 + t162 * t56) * t161;
t24 = -t163 * t50 + (-t160 * t54 + t162 * t58) * t161;
t256 = t22 + t24;
t227 = t160 * t161;
t255 = t82 * t227 + (t86 + t87) * t161 * t162;
t254 = -t267 * t176 + t268 * t178;
t253 = t268 * t176 + t267 * t178;
t252 = t163 ^ 2;
t169 = t176 ^ 2;
t170 = t178 ^ 2;
t251 = m(5) / 0.2e1;
t250 = m(6) / 0.2e1;
t249 = m(7) / 0.2e1;
t246 = -t178 / 0.2e1;
t142 = rSges(3,1) * t175 + rSges(3,2) * t177;
t245 = m(3) * t142;
t244 = pkin(2) * t175;
t243 = pkin(3) * t163;
t172 = cos(pkin(10));
t157 = pkin(4) * t172 + pkin(3);
t242 = -pkin(3) + t157;
t241 = t163 * t263 - t85 * t227 + t255;
t240 = rSges(7,2) * t225 - t266;
t239 = rSges(7,2) * t224 + t117 * t262 + t118 * t264;
t158 = pkin(2) * t177 + pkin(1);
t151 = t178 * t158;
t166 = t178 * pkin(7);
t173 = -qJ(3) - pkin(7);
t216 = t178 * t173;
t237 = t176 * (t216 + t166 + (-pkin(1) + t158) * t176) + t178 * (-t178 * pkin(1) + t151 + (-pkin(7) - t173) * t176);
t236 = rSges(3,1) * t177;
t235 = rSges(3,2) * t175;
t234 = t178 * rSges(3,3);
t233 = -t163 * rSges(7,2) + (t160 * t262 + t162 * t264) * t161;
t232 = Icges(3,4) * t175;
t231 = Icges(3,4) * t177;
t230 = Icges(4,4) * t161;
t229 = Icges(4,4) * t163;
t228 = qJ(4) * t161;
t171 = sin(pkin(10));
t222 = t171 * t178;
t221 = t172 * t178;
t218 = t176 * t171;
t217 = t176 * t172;
t174 = -pkin(8) - qJ(4);
t215 = -pkin(4) * t222 - t174 * t225;
t214 = pkin(3) * t223 + qJ(4) * t224;
t213 = t176 * rSges(3,3) + t178 * t236;
t212 = t169 + t170;
t62 = t118 * rSges(6,1) - t117 * rSges(6,2) + rSges(6,3) * t224;
t129 = -t163 * t222 + t217;
t130 = t163 * t221 + t218;
t211 = t130 * rSges(5,1) + t129 * rSges(5,2) + rSges(5,3) * t224;
t210 = Icges(4,5) * t161 / 0.2e1 + Icges(4,6) * t163 / 0.2e1 + Icges(3,5) * t259 + Icges(3,6) * t258;
t209 = -t161 * pkin(3) + t163 * qJ(4) - t244;
t208 = -rSges(4,1) * t161 - rSges(4,2) * t163 - t244;
t207 = -t176 * t173 + t151;
t206 = -t157 * t163 - t158;
t205 = t169 * (t228 + t243) + t178 * t214 + t237;
t204 = t251 + t250 + t249;
t203 = t209 - (qJ(4) + t174) * t163 - t242 * t161;
t202 = t209 + t163 * rSges(5,3) - (rSges(5,1) * t172 - rSges(5,2) * t171) * t161;
t201 = -t215 - t216;
t200 = -t235 + t236;
t199 = rSges(4,1) * t163 - rSges(4,2) * t161;
t127 = -t163 * t218 - t221;
t128 = t163 * t217 - t222;
t198 = -t128 * rSges(5,1) - t127 * rSges(5,2);
t197 = -t116 * rSges(6,1) + t115 * rSges(6,2);
t89 = -t163 * rSges(6,3) + (rSges(6,1) * t162 - rSges(6,2) * t160) * t161;
t196 = t203 - t89;
t195 = Icges(3,1) * t177 - t232;
t194 = Icges(4,1) * t163 - t230;
t193 = -Icges(3,2) * t175 + t231;
t192 = -Icges(4,2) * t161 + t229;
t183 = pkin(4) * t218 + t157 * t223 - t174 * t224;
t185 = t176 * ((t242 * t163 - t228) * t176 + t215) + t178 * (t183 - t214) + t205;
t184 = rSges(4,1) * t223 - rSges(4,2) * t224 + t176 * rSges(4,3);
t182 = t24 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1 + t22 / 0.2e1;
t181 = t30 / 0.2e1 + t29 / 0.2e1 + t21 / 0.2e1 + t23 / 0.2e1;
t180 = t203 - t233;
t179 = t183 + t207;
t159 = t161 ^ 2;
t144 = rSges(2,1) * t178 - t176 * rSges(2,2);
t143 = -t176 * rSges(2,1) - rSges(2,2) * t178;
t99 = t208 * t178;
t98 = t208 * t176;
t96 = -Icges(5,5) * t163 + (Icges(5,1) * t172 - Icges(5,4) * t171) * t161;
t95 = -Icges(5,6) * t163 + (Icges(5,4) * t172 - Icges(5,2) * t171) * t161;
t93 = t176 * pkin(7) + (pkin(1) - t235) * t178 + t213;
t92 = t234 + t166 + (-pkin(1) - t200) * t176;
t80 = t184 + t207;
t79 = (rSges(4,3) - t173) * t178 + (-t158 - t199) * t176;
t75 = t178 * (-t178 * t235 + t213) + (t176 * t200 - t234) * t176;
t72 = Icges(5,1) * t130 + Icges(5,4) * t129 + Icges(5,5) * t224;
t71 = Icges(5,1) * t128 + Icges(5,4) * t127 + Icges(5,5) * t225;
t70 = Icges(5,4) * t130 + Icges(5,2) * t129 + Icges(5,6) * t224;
t69 = Icges(5,4) * t128 + Icges(5,2) * t127 + Icges(5,6) * t225;
t68 = Icges(5,5) * t130 + Icges(5,6) * t129 + Icges(5,3) * t224;
t67 = Icges(5,5) * t128 + Icges(5,6) * t127 + Icges(5,3) * t225;
t66 = t202 * t178;
t65 = t202 * t176;
t60 = rSges(6,3) * t225 - t197;
t46 = t207 + t211 + t214;
t45 = -t216 + (-t243 - t158 + (-rSges(5,3) - qJ(4)) * t161) * t176 + t198;
t44 = t196 * t178;
t43 = t196 * t176;
t42 = -t163 * t62 - t224 * t89;
t41 = t163 * t60 + t225 * t89;
t40 = t179 + t62;
t39 = (-rSges(6,3) * t161 + t206) * t176 + t197 + t201;
t38 = t178 * t184 + (-t178 * rSges(4,3) + t176 * t199) * t176 + t237;
t35 = t180 * t178;
t34 = t180 * t176;
t33 = (-t176 * t62 + t178 * t60) * t161;
t28 = t179 + t239;
t27 = (-rSges(7,2) * t161 + t206) * t176 + t201 + t266;
t26 = -t163 * t239 - t224 * t233;
t25 = t163 * t240 + t225 * t233;
t20 = t176 * (rSges(5,3) * t225 - t198) + t178 * t211 + t205;
t11 = (-t176 * t239 + t178 * t240) * t161;
t10 = t176 * t60 + t178 * t62 + t185;
t9 = t176 * t240 + t178 * t239 + t185;
t8 = t19 * t176 - t178 * t18;
t7 = -t16 * t178 + t17 * t176;
t6 = -t14 * t178 + t15 * t176;
t5 = -t12 * t178 + t13 * t176;
t1 = [t177 * (Icges(3,2) * t177 + t232) + t175 * (Icges(3,1) * t175 + t231) + Icges(2,3) + (t230 - (Icges(5,5) * t172 - Icges(5,6) * t171) * t161 + (Icges(4,2) + Icges(5,3)) * t163 + t263) * t163 + (Icges(4,1) * t161 - t160 * t85 - t171 * t95 + t172 * t96 + t229) * t161 + m(7) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t79 ^ 2 + t80 ^ 2) + m(3) * (t92 ^ 2 + t93 ^ 2) + m(2) * (t143 ^ 2 + t144 ^ 2) + t255; m(7) * (t27 * t35 + t28 * t34) + m(6) * (t39 * t44 + t40 * t43) + m(5) * (t45 * t66 + t46 * t65) + m(4) * (t79 * t99 + t80 * t98) + (-t127 * t95 / 0.2e1 - t128 * t96 / 0.2e1 - t92 * t245 - t177 * (-Icges(3,6) * t178 + t176 * t193) / 0.2e1 - t175 * (-Icges(3,5) * t178 + t176 * t195) / 0.2e1 + t210 * t178 + (Icges(4,6) * t265 - t176 * t192 / 0.2e1 + t67 / 0.2e1) * t163 - t181) * t178 + (t129 * t95 / 0.2e1 + t130 * t96 / 0.2e1 - t93 * t245 + (Icges(3,6) * t176 + t178 * t193) * t258 + (Icges(3,5) * t176 + t178 * t195) * t259 + t210 * t176 + (Icges(4,6) * t247 + t192 * t265 - t68 / 0.2e1) * t163 + t182) * t176 + ((Icges(4,5) * t176 - t171 * t70 + t172 * t72 + t178 * t194) * t247 + (-Icges(4,5) * t178 - t171 * t69 + t172 * t71 + t176 * t194) * t246) * t161; m(7) * (t34 ^ 2 + t35 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t20 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t38 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(3) * (t142 ^ 2 * t212 + t75 ^ 2) + (-t5 - t6 + (t127 * t69 + t128 * t71 + t67 * t225) * t178 + t254 * t170) * t178 + (t7 + t8 + (t129 * t70 + t130 * t72 + t68 * t224) * t176 + t253 * t169 + (-t127 * t70 - t128 * t72 - t129 * t69 - t130 * t71 + t254 * t176 + t253 * t178 - t67 * t224 - t68 * t225) * t178) * t176; m(7) * (t176 * t27 - t178 * t28) + m(6) * (t176 * t39 - t178 * t40) + m(5) * (t176 * t45 - t178 * t46) + m(4) * (t176 * t79 - t178 * t80); m(7) * (t176 * t35 - t178 * t34) + m(6) * (t176 * t44 - t178 * t43) + m(5) * (t176 * t66 - t178 * t65) + m(4) * (t176 * t99 - t178 * t98); 0.2e1 * (m(4) / 0.2e1 + t204) * t212; 0.2e1 * ((t176 * t28 + t178 * t27) * t249 + (t176 * t40 + t178 * t39) * t250 + (t176 * t46 + t178 * t45) * t251) * t161; m(7) * (-t163 * t9 + (t176 * t34 + t178 * t35) * t161) + m(6) * (-t163 * t10 + (t176 * t43 + t178 * t44) * t161) + m(5) * (-t163 * t20 + (t176 * t65 + t178 * t66) * t161); 0; 0.2e1 * t204 * (t159 * t212 + t252); -t241 * t163 + m(7) * (t25 * t27 + t26 * t28) + m(6) * (t39 * t41 + t40 * t42) + (t176 * t181 + t178 * t182) * t161; m(7) * (t11 * t9 + t25 * t35 + t26 * t34) + m(6) * (t10 * t33 + t41 * t44 + t42 * t43) + ((t8 / 0.2e1 + t7 / 0.2e1) * t178 + (t6 / 0.2e1 + t5 / 0.2e1) * t176) * t161 - (t256 * t176 + t257 * t178) * t163 / 0.2e1 + t260 * t247 + t261 * t246; m(6) * (t41 * t176 - t178 * t42) + m(7) * (t25 * t176 - t178 * t26); m(6) * (-t33 * t163 + (t176 * t42 + t178 * t41) * t161) + m(7) * (-t11 * t163 + (t176 * t26 + t178 * t25) * t161); m(7) * (t11 ^ 2 + t25 ^ 2 + t26 ^ 2) + m(6) * (t33 ^ 2 + t41 ^ 2 + t42 ^ 2) + t241 * t252 + (t260 * t178 + t261 * t176 + (t257 * t176 - t256 * t178) * t163) * t161; m(7) * (t115 * t28 + t117 * t27); m(7) * (t115 * t34 + t117 * t35 + t227 * t9); m(7) * (-t115 * t178 + t117 * t176); m(7) * (t115 * t176 + t117 * t178 - t160 * t163) * t161; m(7) * (t11 * t227 + t115 * t26 + t117 * t25); m(7) * (t159 * t160 ^ 2 + t115 ^ 2 + t117 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
