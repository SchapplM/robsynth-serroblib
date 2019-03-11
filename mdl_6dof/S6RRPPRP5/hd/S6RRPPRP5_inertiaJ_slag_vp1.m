% Calculate joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:31
% EndTime: 2019-03-09 08:41:38
% DurationCPUTime: 3.36s
% Computational Cost: add. (4507->417), mult. (6762->612), div. (0->0), fcn. (7123->8), ass. (0->193)
t276 = Icges(4,1) + Icges(3,3);
t177 = sin(qJ(2));
t179 = cos(qJ(2));
t275 = (-Icges(4,4) + Icges(3,5)) * t179 + (Icges(4,5) - Icges(3,6)) * t177;
t169 = pkin(9) + qJ(5);
t160 = sin(t169);
t180 = cos(qJ(1));
t161 = cos(t169);
t178 = sin(qJ(1));
t232 = t178 * t161;
t105 = t160 * t180 + t177 * t232;
t235 = t160 * t178;
t107 = -t161 * t180 + t177 * t235;
t268 = rSges(7,3) + qJ(6);
t269 = rSges(7,1) + pkin(5);
t274 = t268 * t105 - t269 * t107;
t273 = -t178 / 0.2e1;
t272 = t178 / 0.2e1;
t271 = -t180 / 0.2e1;
t270 = t180 / 0.2e1;
t234 = t161 * t179;
t85 = Icges(7,6) * t177 + (-Icges(7,5) * t160 + Icges(7,3) * t161) * t179;
t86 = Icges(6,3) * t177 + (-Icges(6,5) * t160 - Icges(6,6) * t161) * t179;
t87 = Icges(7,2) * t177 + (-Icges(7,4) * t160 + Icges(7,6) * t161) * t179;
t267 = t85 * t234 + (t86 + t87) * t177;
t88 = Icges(6,6) * t177 + (-Icges(6,4) * t160 - Icges(6,2) * t161) * t179;
t89 = Icges(7,4) * t177 + (-Icges(7,1) * t160 + Icges(7,5) * t161) * t179;
t90 = Icges(6,5) * t177 + (-Icges(6,1) * t160 - Icges(6,4) * t161) * t179;
t266 = -t161 * t88 + (-t89 - t90) * t160;
t233 = t177 * t180;
t103 = -t161 * t233 + t235;
t104 = t160 * t233 + t232;
t228 = t179 * t180;
t47 = Icges(7,5) * t104 + Icges(7,6) * t228 + Icges(7,3) * t103;
t51 = Icges(7,4) * t104 + Icges(7,2) * t228 + Icges(7,6) * t103;
t55 = Icges(7,1) * t104 + Icges(7,4) * t228 + Icges(7,5) * t103;
t12 = t103 * t47 + t104 * t55 + t51 * t228;
t229 = t178 * t179;
t48 = Icges(7,5) * t107 + Icges(7,6) * t229 - Icges(7,3) * t105;
t52 = Icges(7,4) * t107 + Icges(7,2) * t229 - Icges(7,6) * t105;
t56 = Icges(7,1) * t107 + Icges(7,4) * t229 - Icges(7,5) * t105;
t13 = t103 * t48 + t104 * t56 + t52 * t228;
t49 = Icges(6,5) * t104 - Icges(6,6) * t103 + Icges(6,3) * t228;
t53 = Icges(6,4) * t104 - Icges(6,2) * t103 + Icges(6,6) * t228;
t57 = Icges(6,1) * t104 - Icges(6,4) * t103 + Icges(6,5) * t228;
t14 = -t103 * t53 + t104 * t57 + t49 * t228;
t50 = Icges(6,5) * t107 + Icges(6,6) * t105 + Icges(6,3) * t229;
t54 = Icges(6,4) * t107 + Icges(6,2) * t105 + Icges(6,6) * t229;
t58 = Icges(6,1) * t107 + Icges(6,4) * t105 + Icges(6,5) * t229;
t15 = -t103 * t54 + t104 * t58 + t50 * t228;
t29 = t103 * t85 + t104 * t89 + t87 * t228;
t30 = -t103 * t88 + t104 * t90 + t86 * t228;
t265 = ((t12 + t14) * t180 + (t13 + t15) * t178) * t179 + (t29 + t30) * t177;
t16 = -t105 * t47 + t107 * t55 + t51 * t229;
t17 = -t105 * t48 + t107 * t56 + t52 * t229;
t18 = t105 * t53 + t107 * t57 + t49 * t229;
t19 = t105 * t54 + t107 * t58 + t50 * t229;
t31 = -t105 * t85 + t107 * t89 + t87 * t229;
t32 = t105 * t88 + t107 * t90 + t86 * t229;
t264 = ((t16 + t18) * t180 + (t17 + t19) * t178) * t179 + (t31 + t32) * t177;
t263 = t177 / 0.2e1;
t20 = t177 * t51 + (-t160 * t55 + t161 * t47) * t179;
t22 = t177 * t49 + (-t160 * t57 - t161 * t53) * t179;
t262 = t20 + t22;
t21 = t177 * t52 + (-t160 * t56 + t161 * t48) * t179;
t23 = t177 * t50 + (-t160 * t58 - t161 * t54) * t179;
t261 = t21 + t23;
t260 = t178 * t276 + t275 * t180;
t259 = -t275 * t178 + t180 * t276;
t171 = t178 ^ 2;
t173 = t180 ^ 2;
t258 = 0.2e1 * t179;
t257 = m(4) / 0.2e1;
t256 = m(5) / 0.2e1;
t255 = m(6) / 0.2e1;
t254 = m(7) / 0.2e1;
t145 = rSges(3,1) * t177 + rSges(3,2) * t179;
t250 = m(3) * t145;
t174 = sin(pkin(9));
t249 = pkin(4) * t174;
t248 = (t179 * t266 + t267) * t177;
t247 = rSges(7,2) * t228 + t268 * t103 + t269 * t104;
t246 = rSges(7,2) * t229 - t274;
t243 = t180 * rSges(4,1);
t242 = t180 * rSges(3,3);
t241 = t177 * rSges(7,2) + (-t269 * t160 + t268 * t161) * t179;
t240 = Icges(3,4) * t177;
t239 = Icges(3,4) * t179;
t238 = Icges(4,6) * t177;
t237 = Icges(4,6) * t179;
t236 = qJ(3) * t177;
t231 = t178 * t174;
t175 = cos(pkin(9));
t230 = t178 * t175;
t224 = pkin(2) * t228 + qJ(3) * t233;
t227 = t171 * (pkin(2) * t179 + t236) + t180 * t224;
t143 = pkin(2) * t177 - qJ(3) * t179;
t226 = rSges(4,2) * t177 + rSges(4,3) * t179 - t143;
t159 = pkin(4) * t175 + pkin(3);
t176 = -pkin(8) - qJ(4);
t225 = -t180 * t159 - t176 * t229;
t223 = t180 * pkin(1) + t178 * pkin(7);
t222 = t178 * pkin(3) + qJ(4) * t228;
t221 = t171 + t173;
t60 = t104 * rSges(6,1) - t103 * rSges(6,2) + rSges(6,3) * t228;
t220 = t174 * t233;
t128 = t175 * t233 - t231;
t129 = t220 + t230;
t219 = t129 * rSges(5,1) + t128 * rSges(5,2) + rSges(5,3) * t228;
t166 = t180 * pkin(7);
t218 = t166 - t225;
t217 = -Icges(4,4) * t177 / 0.2e1 + Icges(3,5) * t263 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t179;
t216 = -qJ(4) * t177 - t143;
t215 = t256 + t255 + t254;
t167 = t180 * pkin(3);
t214 = t178 * (qJ(4) * t229 - t167) + t180 * t222 + t227;
t213 = t223 + t224;
t212 = -t177 * rSges(5,3) - (-rSges(5,1) * t174 - rSges(5,2) * t175) * t179 + t216;
t211 = t179 * t249 - (-qJ(4) - t176) * t177 + t216;
t210 = rSges(3,1) * t179 - rSges(3,2) * t177;
t130 = t174 * t180 + t177 * t230;
t131 = -t175 * t180 + t177 * t231;
t209 = -t131 * rSges(5,1) - t130 * rSges(5,2);
t208 = -t107 * rSges(6,1) - t105 * rSges(6,2);
t25 = -t246 * t177 + t241 * t229;
t26 = t247 * t177 - t241 * t228;
t207 = t178 * t26 + t180 * t25;
t183 = t211 - t241;
t38 = t183 * t178;
t39 = t183 * t180;
t206 = t178 * t38 + t180 * t39;
t62 = rSges(6,3) * t229 - t208;
t93 = t177 * rSges(6,3) + (-rSges(6,1) * t160 - rSges(6,2) * t161) * t179;
t40 = -t177 * t62 + t93 * t229;
t41 = t177 * t60 - t93 * t228;
t205 = t178 * t41 + t180 * t40;
t189 = t211 - t93;
t45 = t189 * t178;
t46 = t189 * t180;
t204 = t178 * t45 + t180 * t46;
t64 = t212 * t178;
t65 = t212 * t180;
t203 = t178 * t64 + t180 * t65;
t202 = Icges(3,1) * t179 - t240;
t201 = -Icges(3,2) * t177 + t239;
t198 = -Icges(4,2) * t179 + t238;
t197 = Icges(4,3) * t177 - t237;
t196 = t103 * t180 - t105 * t178;
t191 = rSges(3,1) * t228 - rSges(3,2) * t233 + t178 * rSges(3,3);
t190 = t178 * rSges(4,1) - rSges(4,2) * t228 + rSges(4,3) * t233;
t188 = pkin(4) * t220 + t178 * t159 - t176 * t228;
t187 = t22 / 0.2e1 + t20 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t186 = t23 / 0.2e1 + t21 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1;
t185 = t178 * (t167 + (-qJ(4) * t179 + t177 * t249) * t178 + t225) + t180 * (t188 - t222) + t214;
t184 = (-qJ(3) - t249) * t177 - pkin(1);
t182 = t188 + t213;
t27 = ((-rSges(7,2) - pkin(2)) * t179 + t184) * t178 + t218 + t274;
t28 = t182 + t247;
t34 = ((-rSges(6,3) - pkin(2)) * t179 + t184) * t178 + t208 + t218;
t35 = t182 + t60;
t43 = t166 + t167 + (-t236 - pkin(1) + (-rSges(5,3) - pkin(2) - qJ(4)) * t179) * t178 + t209;
t44 = t213 + t219 + t222;
t181 = (t178 * t28 + t180 * t27) * t254 + (t178 * t35 + t180 * t34) * t255 + (t178 * t44 + t180 * t43) * t256;
t172 = t179 ^ 2;
t170 = t177 ^ 2;
t147 = rSges(2,1) * t180 - t178 * rSges(2,2);
t146 = -t178 * rSges(2,1) - rSges(2,2) * t180;
t101 = Icges(5,5) * t177 + (-Icges(5,1) * t174 - Icges(5,4) * t175) * t179;
t100 = Icges(5,6) * t177 + (-Icges(5,4) * t174 - Icges(5,2) * t175) * t179;
t84 = t226 * t180;
t83 = t226 * t178;
t80 = t191 + t223;
t79 = t242 + t166 + (-pkin(1) - t210) * t178;
t77 = t190 + t213;
t76 = t243 + t166 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t179 + (-rSges(4,3) - qJ(3)) * t177) * t178;
t73 = Icges(5,1) * t131 + Icges(5,4) * t130 + Icges(5,5) * t229;
t72 = Icges(5,1) * t129 + Icges(5,4) * t128 + Icges(5,5) * t228;
t71 = Icges(5,4) * t131 + Icges(5,2) * t130 + Icges(5,6) * t229;
t70 = Icges(5,4) * t129 + Icges(5,2) * t128 + Icges(5,6) * t228;
t69 = Icges(5,5) * t131 + Icges(5,6) * t130 + Icges(5,3) * t229;
t68 = Icges(5,5) * t129 + Icges(5,6) * t128 + Icges(5,3) * t228;
t63 = t180 * t191 + (t210 * t178 - t242) * t178;
t42 = t180 * t190 + (-t243 + (-rSges(4,2) * t179 + rSges(4,3) * t177) * t178) * t178 + t227;
t33 = (-t178 * t60 + t180 * t62) * t179;
t24 = t178 * (rSges(5,3) * t229 - t209) + t180 * t219 + t214;
t11 = (-t247 * t178 + t246 * t180) * t179;
t10 = t178 * t62 + t180 * t60 + t185;
t9 = t246 * t178 + t247 * t180 + t185;
t8 = t18 * t178 - t180 * t19;
t7 = t16 * t178 - t17 * t180;
t6 = t14 * t178 - t15 * t180;
t5 = t12 * t178 - t13 * t180;
t1 = [Icges(2,3) + (t239 + t237 + (-Icges(5,5) * t174 - Icges(5,6) * t175) * t179 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t177) * t177 + (-t175 * t100 - t174 * t101 + t238 + t240 + (Icges(4,3) + Icges(3,2)) * t179 + t266) * t179 + m(7) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t76 ^ 2 + t77 ^ 2) + m(3) * (t79 ^ 2 + t80 ^ 2) + m(2) * (t146 ^ 2 + t147 ^ 2) + t267; m(4) * (t76 * t84 + t77 * t83) + m(5) * (t43 * t65 + t44 * t64) + m(6) * (t34 * t46 + t35 * t45) + m(7) * (t27 * t39 + t28 * t38) + (-t130 * t100 / 0.2e1 - t131 * t101 / 0.2e1 - t79 * t250 + t217 * t180 + (-t69 / 0.2e1 + Icges(3,5) * t270 + t202 * t273 + Icges(4,4) * t271 + t198 * t272) * t177 - t186) * t180 + (t128 * t100 / 0.2e1 + t129 * t101 / 0.2e1 - t80 * t250 + t217 * t178 + (t68 / 0.2e1 + Icges(3,5) * t272 + t202 * t270 + Icges(4,4) * t273 + t198 * t271) * t177 + t187) * t178 + ((t174 * t73 / 0.2e1 + t175 * t71 / 0.2e1 + Icges(3,6) * t270 + t201 * t273 + Icges(4,5) * t271 + t197 * t272) * t180 + (-t174 * t72 / 0.2e1 - t175 * t70 / 0.2e1 + Icges(3,6) * t272 + t201 * t270 + Icges(4,5) * t273 + t197 * t271) * t178) * t179; m(7) * (t38 ^ 2 + t39 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t24 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t42 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(3) * (t221 * t145 ^ 2 + t63 ^ 2) + (-t7 - t8 + (t130 * t71 + t131 * t73 + t69 * t229) * t180 + t259 * t173) * t180 + (t6 + t5 + (t128 * t70 + t129 * t72 + t68 * t228) * t178 + t260 * t171 + (-t128 * t71 - t129 * t73 - t130 * t70 - t131 * t72 + t178 * t259 + t180 * t260 - t69 * t228 - t68 * t229) * t180) * t178; 0.2e1 * ((t178 * t77 + t180 * t76) * t257 + t181) * t177; m(7) * (t206 * t177 - t179 * t9) + m(6) * (-t179 * t10 + t204 * t177) + m(5) * (t203 * t177 - t179 * t24) + m(4) * (-t179 * t42 + (t178 * t83 + t180 * t84) * t177); 0.2e1 * (t257 + t215) * (t221 * t170 + t172); t181 * t258; m(7) * (t177 * t9 + t206 * t179) + m(6) * (t177 * t10 + t204 * t179) + m(5) * (t177 * t24 + t203 * t179); t215 * (-0.1e1 + t221) * t177 * t258; 0.2e1 * t215 * (t221 * t172 + t170); m(7) * (t25 * t27 + t26 * t28) + m(6) * (t34 * t40 + t35 * t41) + (t186 * t178 + t187 * t180) * t179 + t248; m(7) * (t11 * t9 + t25 * t39 + t26 * t38) + m(6) * (t10 * t33 + t40 * t46 + t41 * t45) + ((t6 / 0.2e1 + t5 / 0.2e1) * t180 + (t8 / 0.2e1 + t7 / 0.2e1) * t178) * t179 + (t178 * t262 - t180 * t261) * t263 + t265 * t272 + t264 * t271; m(6) * (t205 * t177 - t33 * t179) + m(7) * (-t11 * t179 + t207 * t177); m(6) * (t33 * t177 + t205 * t179) + m(7) * (t11 * t177 + t207 * t179); t248 * t177 + m(7) * (t11 ^ 2 + t25 ^ 2 + t26 ^ 2) + m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t265 * t180 + t264 * t178 + (t178 * t261 + t180 * t262) * t177) * t179; m(7) * (t103 * t27 - t105 * t28); m(7) * (t103 * t39 - t105 * t38 + t9 * t234); m(7) * (-t172 * t161 + t196 * t177); m(7) * (t161 * t177 + t196) * t179; m(7) * (t103 * t25 - t105 * t26 + t11 * t234); m(7) * (t161 ^ 2 * t172 + t103 ^ 2 + t105 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
