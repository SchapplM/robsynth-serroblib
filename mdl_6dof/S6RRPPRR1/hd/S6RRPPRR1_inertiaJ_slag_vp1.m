% Calculate joint inertia matrix for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:54
% EndTime: 2019-03-09 08:46:02
% DurationCPUTime: 3.50s
% Computational Cost: add. (7333->369), mult. (9653->554), div. (0->0), fcn. (11227->10), ass. (0->186)
t276 = Icges(4,1) + Icges(5,1);
t274 = Icges(5,4) + Icges(4,5);
t153 = qJ(2) + pkin(10);
t147 = sin(t153);
t275 = (Icges(4,4) - Icges(5,5)) * t147;
t148 = cos(t153);
t273 = t276 * t148 - t275;
t272 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t159 = sin(qJ(2));
t162 = cos(qJ(2));
t271 = Icges(3,5) * t162 - Icges(3,6) * t159 + t274 * t148 + (-Icges(4,6) + Icges(5,6)) * t147;
t160 = sin(qJ(1));
t270 = -t160 / 0.2e1;
t269 = t160 / 0.2e1;
t163 = cos(qJ(1));
t268 = -t163 / 0.2e1;
t261 = t163 / 0.2e1;
t267 = t160 * pkin(7);
t158 = sin(qJ(5));
t253 = cos(qJ(5));
t210 = t147 * t253;
t219 = t148 * t160;
t104 = t158 * t219 - t160 * t210;
t116 = t147 * t158 + t148 * t253;
t105 = t116 * t160;
t155 = t163 ^ 2;
t224 = Icges(6,6) * t163;
t225 = Icges(6,2) * t104;
t226 = Icges(6,5) * t163;
t229 = Icges(6,4) * t105;
t234 = Icges(6,1) * t105;
t218 = t148 * t163;
t106 = t158 * t218 - t163 * t210;
t107 = t116 * t163;
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t77 = -t107 * t157 - t160 * t161;
t78 = t107 * t161 - t157 * t160;
t40 = Icges(7,5) * t78 + Icges(7,6) * t77 + Icges(7,3) * t106;
t41 = Icges(7,4) * t78 + Icges(7,2) * t77 + Icges(7,6) * t106;
t42 = Icges(7,1) * t78 + Icges(7,4) * t77 + Icges(7,5) * t106;
t75 = -t105 * t157 + t161 * t163;
t76 = t105 * t161 + t157 * t163;
t11 = t104 * t40 + t41 * t75 + t42 * t76;
t222 = Icges(7,3) * t104;
t239 = Icges(7,5) * t76;
t191 = t222 + 0.2e1 * t239;
t238 = Icges(7,2) * t75;
t240 = Icges(7,4) * t76;
t197 = t238 + 0.2e1 * t240;
t237 = Icges(7,6) * t75;
t241 = Icges(7,1) * t76 ^ 2;
t3 = t11 * t160 - (t241 + t197 * t75 + (t191 + 0.2e1 * t237) * t104) * t163;
t55 = Icges(6,5) * t107 - Icges(6,6) * t106 - Icges(6,3) * t160;
t56 = Icges(6,4) * t107 - Icges(6,2) * t106 - Icges(6,6) * t160;
t57 = Icges(6,1) * t107 - Icges(6,4) * t106 - Icges(6,5) * t160;
t255 = -(-t104 * t56 + t105 * t57 + t163 * t55) * t160 + (Icges(6,3) * t155 + (0.2e1 * t226 + t234) * t105 + (-0.2e1 * t224 + t225 - 0.2e1 * t229) * t104) * t163 - t3;
t167 = t224 - t225 + t229;
t168 = -Icges(6,4) * t104 + t226 + t234;
t12 = t106 * t40 + t41 * t77 + t42 * t78;
t169 = t222 + t237 + t239;
t223 = Icges(7,6) * t104;
t170 = t223 + t238 + t240;
t171 = Icges(7,1) * t76 + Icges(7,4) * t75 + Icges(7,5) * t104;
t164 = t106 * t169 + t77 * t170 + t78 * t171;
t4 = t12 * t160 - t164 * t163;
t254 = (-t106 * t56 + t107 * t57 - t160 * t55) * t160 - (t107 * t168 - t106 * t167 - t160 * (Icges(6,5) * t105 - Icges(6,6) * t104 + Icges(6,3) * t163)) * t163 + t4;
t266 = -t104 / 0.2e1;
t265 = -t116 / 0.2e1;
t263 = t159 / 0.2e1;
t262 = t162 / 0.2e1;
t260 = -t160 * t271 + t272 * t163;
t259 = t160 * t272 + t271 * t163;
t44 = t78 * rSges(7,1) + t77 * rSges(7,2) + t106 * rSges(7,3);
t248 = t107 * pkin(5) + pkin(9) * t106 + t44;
t251 = t105 * pkin(5);
t202 = -t76 * rSges(7,1) - t75 * rSges(7,2);
t43 = t104 * rSges(7,3) - t202;
t15 = -(t104 * pkin(9) + t251 + t43) * t160 - t248 * t163;
t258 = -t254 * t160 - t255 * t163;
t154 = t160 ^ 2;
t257 = m(5) / 0.2e1;
t256 = m(7) / 0.2e1;
t252 = pkin(2) * t159;
t156 = -qJ(3) - pkin(7);
t250 = -pkin(8) - t156;
t117 = -t148 * t158 + t210;
t49 = Icges(7,3) * t116 + (Icges(7,5) * t161 - Icges(7,6) * t157) * t117;
t51 = Icges(7,5) * t116 + (Icges(7,1) * t161 - Icges(7,4) * t157) * t117;
t247 = t117 * t161 * t51 + t116 * t49;
t52 = rSges(7,3) * t116 + (rSges(7,1) * t161 - rSges(7,2) * t157) * t117;
t246 = pkin(5) * t117 + pkin(9) * t116 + t52;
t145 = pkin(2) * t162 + pkin(1);
t142 = t163 * t145;
t152 = t163 * pkin(7);
t245 = t160 * (t152 + (-pkin(1) + t145) * t160) + t163 * (-pkin(1) * t163 + t142 - t267);
t244 = t107 * rSges(6,1) - t106 * rSges(6,2);
t243 = rSges(3,1) * t162;
t242 = rSges(3,2) * t159;
t50 = Icges(7,6) * t116 + (Icges(7,4) * t161 - Icges(7,2) * t157) * t117;
t236 = t157 * t50;
t235 = t163 * rSges(3,3);
t233 = Icges(3,4) * t159;
t232 = Icges(3,4) * t162;
t230 = Icges(4,4) * t148;
t227 = Icges(5,5) * t148;
t221 = qJ(4) * t147;
t220 = t147 * t163;
t217 = pkin(3) * t218 + qJ(4) * t220;
t216 = t160 * rSges(3,3) + t163 * t243;
t215 = t154 + t155;
t214 = -rSges(6,3) + t250;
t13 = t116 * t169 + (-t157 * t170 + t161 * t171) * t117;
t16 = t104 * t49 + t50 * t75 + t51 * t76;
t213 = t16 / 0.2e1 + t13 / 0.2e1;
t14 = t116 * t40 + (-t157 * t41 + t161 * t42) * t117;
t17 = t106 * t49 + t50 * t77 + t51 * t78;
t212 = t17 / 0.2e1 + t14 / 0.2e1;
t211 = rSges(5,1) * t218 + t160 * rSges(5,2) + rSges(5,3) * t220;
t209 = -pkin(3) * t147 + qJ(4) * t148 - t252;
t208 = -rSges(4,1) * t147 - rSges(4,2) * t148 - t252;
t207 = -t160 * t156 + t142;
t206 = t154 * (pkin(3) * t148 + t221) + t163 * t217 + t245;
t205 = t257 + m(6) / 0.2e1 + t256;
t140 = pkin(4) * t218;
t204 = t140 + t142 + t217;
t203 = -rSges(5,1) * t147 + rSges(5,3) * t148 + t209;
t201 = Icges(3,5) * t263 + Icges(3,6) * t262 + t274 * t147 / 0.2e1 + (-Icges(5,6) / 0.2e1 + Icges(4,6) / 0.2e1) * t148;
t200 = -t242 + t243;
t199 = rSges(4,1) * t148 - rSges(4,2) * t147;
t198 = -t105 * rSges(6,1) + t104 * rSges(6,2);
t178 = -pkin(4) * t147 + t209;
t71 = rSges(6,1) * t117 - rSges(6,2) * t116;
t173 = t178 - t71;
t47 = t173 * t160;
t48 = t173 * t163;
t192 = t160 * t47 + t163 * t48;
t34 = t160 * t198 - t163 * t244;
t190 = Icges(3,1) * t162 - t233;
t187 = -Icges(3,2) * t159 + t232;
t186 = -Icges(4,2) * t147 + t230;
t182 = Icges(5,3) * t147 + t227;
t179 = rSges(4,1) * t218 - rSges(4,2) * t220 + t160 * rSges(4,3);
t68 = Icges(6,5) * t117 - Icges(6,6) * t116;
t69 = Icges(6,4) * t117 - Icges(6,2) * t116;
t70 = Icges(6,1) * t117 - Icges(6,4) * t116;
t177 = t167 * t265 + t117 * t168 / 0.2e1 + t69 * t266 + t105 * t70 / 0.2e1 + t68 * t261 + t213;
t176 = t106 * t69 / 0.2e1 - t107 * t70 / 0.2e1 + t68 * t269 + t116 * t56 / 0.2e1 - t117 * t57 / 0.2e1 - t212;
t175 = pkin(4) * t160 * t219 + t140 * t163 + t206;
t166 = (-t221 - t145 + (-pkin(3) - pkin(4)) * t148) * t160;
t38 = t214 * t163 + t166 + t198;
t39 = t214 * t160 + t204 + t244;
t174 = m(6) * (t160 * t39 + t163 * t38);
t172 = t178 - t246;
t1 = t11 * t106 + t16 * t116 + (t241 + t191 * t104 + (t197 + 0.2e1 * t223) * t75) * t104;
t2 = t164 * t104 + t12 * t106 + t17 * t116;
t165 = t3 * t266 - t106 * t4 / 0.2e1 + (-t13 * t163 + t14 * t160) * t265 + t2 * t270 + t1 * t261;
t135 = rSges(2,1) * t163 - t160 * rSges(2,2);
t134 = -t160 * rSges(2,1) - rSges(2,2) * t163;
t133 = rSges(3,1) * t159 + rSges(3,2) * t162;
t84 = t208 * t163;
t83 = t208 * t160;
t82 = t267 + (pkin(1) - t242) * t163 + t216;
t81 = t235 + t152 + (-pkin(1) - t200) * t160;
t67 = t179 + t207;
t66 = (rSges(4,3) - t156) * t163 + (-t145 - t199) * t160;
t64 = t203 * t163;
t63 = t203 * t160;
t62 = t163 * (-t163 * t242 + t216) + (t200 * t160 - t235) * t160;
t54 = t207 + t211 + t217;
t53 = (rSges(5,2) - t156) * t163 + (-t145 + (-rSges(5,1) - pkin(3)) * t148 + (-rSges(5,3) - qJ(4)) * t147) * t160;
t37 = t246 * t163;
t36 = t246 * t160;
t35 = t163 * t179 + (-t163 * rSges(4,3) + t199 * t160) * t160 + t245;
t33 = t172 * t163;
t32 = t172 * t160;
t29 = t163 * t211 + (-t163 * rSges(5,2) + (rSges(5,1) * t148 + rSges(5,3) * t147) * t160) * t160 + t206;
t26 = t250 * t160 + t204 + t248;
t25 = -t251 + t250 * t163 + (-rSges(7,3) - pkin(9)) * t104 + t166 + t202;
t24 = -t106 * t52 + t116 * t44;
t23 = t104 * t52 - t116 * t43;
t20 = -t104 * t44 + t106 * t43;
t19 = t175 - t34;
t18 = (-t117 * t236 + t247) * t116;
t6 = -t15 + t175;
t5 = [-t116 * t69 + t162 * (Icges(3,2) * t162 + t233) + t159 * (Icges(3,1) * t159 + t232) + Icges(2,3) + (t70 - t236) * t117 + m(7) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t53 ^ 2 + t54 ^ 2) + m(3) * (t81 ^ 2 + t82 ^ 2) + m(4) * (t66 ^ 2 + t67 ^ 2) + m(2) * (t134 ^ 2 + t135 ^ 2) + t247 + ((Icges(4,2) + Icges(5,3)) * t148 + t275) * t148 + (t276 * t147 - t227 + t230) * t147; (-t162 * (-Icges(3,6) * t163 + t187 * t160) / 0.2e1 - t159 * (-Icges(3,5) * t163 + t190 * t160) / 0.2e1 + t201 * t163 + (Icges(4,6) * t261 + Icges(5,6) * t268 + t182 * t269 + t186 * t270) * t148 + (t261 * t274 + t273 * t270) * t147 - t177) * t163 + ((Icges(3,6) * t160 + t187 * t163) * t262 + (Icges(3,5) * t160 + t190 * t163) * t263 + t201 * t160 + (Icges(4,6) * t269 + Icges(5,6) * t270 + t182 * t268 + t186 * t261) * t148 + (t273 * t261 + t269 * t274) * t147 - t176) * t160 + m(7) * (t33 * t25 + t32 * t26) + m(6) * (t48 * t38 + t47 * t39) + m(5) * (t64 * t53 + t63 * t54) + m(4) * (t84 * t66 + t83 * t67) + m(3) * (-t160 * t82 - t163 * t81) * t133; m(7) * (t32 ^ 2 + t33 ^ 2 + t6 ^ 2) + m(6) * (t19 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(5) * (t29 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(4) * (t35 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(3) * (t215 * t133 ^ 2 + t62 ^ 2) + (t155 * t260 + t255) * t163 + (t259 * t154 + (t160 * t260 + t163 * t259) * t163 + t254) * t160; m(7) * (t160 * t25 - t163 * t26) + m(6) * (t160 * t38 - t163 * t39) + m(5) * (t160 * t53 - t163 * t54) + m(4) * (t160 * t66 - t163 * t67); m(7) * (t160 * t33 - t163 * t32) + m(6) * (t160 * t48 - t163 * t47) + m(5) * (t160 * t64 - t163 * t63) + m(4) * (t160 * t84 - t163 * t83); 0.2e1 * (m(4) / 0.2e1 + t205) * t215; 0.2e1 * ((t160 * t26 + t163 * t25) * t256 + t174 / 0.2e1 + (t160 * t54 + t163 * t53) * t257) * t147; m(7) * (-t148 * t6 + (t160 * t32 + t163 * t33) * t147) + m(6) * (t147 * t192 - t148 * t19) + m(5) * (-t148 * t29 + (t160 * t63 + t163 * t64) * t147); 0; 0.2e1 * t205 * (t147 ^ 2 * t215 + t148 ^ 2); t177 * t163 + t176 * t160 + m(7) * (t25 * t37 + t26 * t36) + t71 * t174; m(7) * (t15 * t6 + t32 * t36 + t33 * t37) + m(6) * (t34 * t19 + t192 * t71) + t258; m(7) * (t37 * t160 - t163 * t36); m(6) * (t147 * t215 * t71 - t148 * t34) + m(7) * (-t15 * t148 + (t160 * t36 + t163 * t37) * t147); m(6) * (t215 * t71 ^ 2 + t34 ^ 2) + m(7) * (t15 ^ 2 + t36 ^ 2 + t37 ^ 2) - t258; m(7) * (t23 * t25 + t24 * t26) + t18 + t212 * t106 + t213 * t104; m(7) * (t20 * t6 + t23 * t33 + t24 * t32) - t165; m(7) * (t23 * t160 - t163 * t24); m(7) * (-t20 * t148 + (t160 * t24 + t163 * t23) * t147); m(7) * (t15 * t20 + t23 * t37 + t24 * t36) + t165; t106 * t2 + t104 * t1 + t116 * (t13 * t104 + t14 * t106 + t18) + m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
