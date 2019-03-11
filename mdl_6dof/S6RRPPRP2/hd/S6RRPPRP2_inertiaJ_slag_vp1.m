% Calculate joint inertia matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:51
% EndTime: 2019-03-09 08:29:58
% DurationCPUTime: 3.17s
% Computational Cost: add. (4646->392), mult. (5882->562), div. (0->0), fcn. (6081->8), ass. (0->179)
t161 = -qJ(6) - pkin(8);
t262 = rSges(7,3) - t161;
t261 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t158 = qJ(2) + pkin(9);
t150 = sin(t158);
t151 = cos(t158);
t164 = sin(qJ(2));
t167 = cos(qJ(2));
t260 = Icges(3,5) * t167 - Icges(3,6) * t164 + (-Icges(5,4) + Icges(4,5)) * t151 + (Icges(5,5) - Icges(4,6)) * t150;
t165 = sin(qJ(1));
t259 = -t165 / 0.2e1;
t258 = t165 / 0.2e1;
t168 = cos(qJ(1));
t257 = -t168 / 0.2e1;
t256 = t168 / 0.2e1;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t76 = Icges(7,3) * t150 + (-Icges(7,5) * t163 - Icges(7,6) * t166) * t151;
t77 = Icges(6,3) * t150 + (-Icges(6,5) * t163 - Icges(6,6) * t166) * t151;
t255 = (t76 + t77) * t150;
t78 = Icges(7,6) * t150 + (-Icges(7,4) * t163 - Icges(7,2) * t166) * t151;
t79 = Icges(6,6) * t150 + (-Icges(6,4) * t163 - Icges(6,2) * t166) * t151;
t80 = Icges(7,5) * t150 + (-Icges(7,1) * t163 - Icges(7,4) * t166) * t151;
t81 = Icges(6,5) * t150 + (-Icges(6,1) * t163 - Icges(6,4) * t166) * t151;
t254 = (-t78 - t79) * t166 + (-t80 - t81) * t163;
t213 = t166 * t168;
t215 = t165 * t163;
t112 = t150 * t213 - t215;
t216 = t163 * t168;
t208 = t150 * t216;
t214 = t165 * t166;
t113 = t208 + t214;
t218 = t151 * t168;
t47 = Icges(7,5) * t113 + Icges(7,6) * t112 + Icges(7,3) * t218;
t51 = Icges(7,4) * t113 + Icges(7,2) * t112 + Icges(7,6) * t218;
t55 = Icges(7,1) * t113 + Icges(7,4) * t112 + Icges(7,5) * t218;
t12 = t112 * t51 + t113 * t55 + t47 * t218;
t114 = t150 * t214 + t216;
t115 = t150 * t215 - t213;
t219 = t151 * t165;
t48 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t219;
t52 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t219;
t56 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t219;
t13 = t112 * t52 + t113 * t56 + t48 * t218;
t49 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t218;
t53 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t218;
t57 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t218;
t14 = t112 * t53 + t113 * t57 + t49 * t218;
t50 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t219;
t54 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t219;
t58 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t219;
t15 = t112 * t54 + t113 * t58 + t50 * t218;
t27 = t112 * t78 + t113 * t80 + t76 * t218;
t28 = t112 * t79 + t113 * t81 + t77 * t218;
t253 = ((t12 + t14) * t168 + (t13 + t15) * t165) * t151 + (t27 + t28) * t150;
t16 = t114 * t51 + t115 * t55 + t47 * t219;
t17 = t114 * t52 + t115 * t56 + t48 * t219;
t18 = t114 * t53 + t115 * t57 + t49 * t219;
t19 = t114 * t54 + t115 * t58 + t50 * t219;
t29 = t114 * t78 + t115 * t80 + t76 * t219;
t30 = t114 * t79 + t115 * t81 + t77 * t219;
t252 = ((t16 + t18) * t168 + (t17 + t19) * t165) * t151 + (t29 + t30) * t150;
t251 = t150 / 0.2e1;
t250 = t164 / 0.2e1;
t249 = t167 / 0.2e1;
t22 = t150 * t47 + (-t163 * t55 - t166 * t51) * t151;
t24 = t150 * t49 + (-t163 * t57 - t166 * t53) * t151;
t248 = t22 + t24;
t23 = t150 * t48 + (-t163 * t56 - t166 * t52) * t151;
t25 = t150 * t50 + (-t163 * t58 - t166 * t54) * t151;
t247 = t23 + t25;
t146 = pkin(5) * t166 + pkin(4);
t246 = t113 * rSges(7,1) + t112 * rSges(7,2) + pkin(5) * t208 + t165 * t146 + t262 * t218;
t245 = -t260 * t165 + t261 * t168;
t244 = t261 * t165 + t260 * t168;
t159 = t165 ^ 2;
t160 = t168 ^ 2;
t243 = m(5) / 0.2e1;
t242 = m(6) / 0.2e1;
t238 = pkin(2) * t164;
t237 = pkin(5) * t163;
t236 = -pkin(8) - t161;
t235 = (t151 * t254 + t255) * t150;
t210 = t165 * pkin(4) + pkin(8) * t218;
t234 = -t210 + t246;
t157 = t168 * pkin(4);
t197 = -t115 * rSges(7,1) - t114 * rSges(7,2);
t233 = rSges(7,3) * t219 - t197 - t146 * t168 + t157 + (t150 * t237 + t236 * t151) * t165;
t147 = pkin(2) * t167 + pkin(1);
t143 = t168 * t147;
t156 = t168 * pkin(7);
t162 = -qJ(3) - pkin(7);
t217 = t162 * t168;
t232 = t165 * (t217 + t156 + (-pkin(1) + t147) * t165) + t168 * (-pkin(1) * t168 + t143 + (-pkin(7) - t162) * t165);
t231 = (-rSges(7,1) * t163 - rSges(7,2) * t166 - t237) * t151 + (rSges(7,3) + t236) * t150;
t230 = rSges(3,1) * t167;
t229 = rSges(3,2) * t164;
t228 = t168 * rSges(3,3);
t227 = Icges(3,4) * t164;
t226 = Icges(3,4) * t167;
t225 = Icges(4,4) * t150;
t224 = Icges(4,4) * t151;
t223 = Icges(5,6) * t150;
t222 = Icges(5,6) * t151;
t221 = qJ(4) * t150;
t220 = t150 * t168;
t212 = pkin(3) * t218 + qJ(4) * t220;
t211 = t165 * rSges(3,3) + t168 * t230;
t209 = t159 + t160;
t60 = t113 * rSges(6,1) + t112 * rSges(6,2) + rSges(6,3) * t218;
t207 = -pkin(3) * t150 + qJ(4) * t151 - t238;
t206 = -rSges(4,1) * t150 - rSges(4,2) * t151 - t238;
t205 = -t165 * t162 + t143;
t204 = t159 * (pkin(3) * t151 + t221) + t168 * t212 + t232;
t203 = t243 + t242 + m(7) / 0.2e1;
t202 = rSges(5,2) * t150 + rSges(5,3) * t151 + t207;
t201 = -Icges(5,4) * t150 / 0.2e1 + Icges(4,5) * t251 + Icges(3,5) * t250 + Icges(3,6) * t249 + (-Icges(5,5) / 0.2e1 + Icges(4,6) / 0.2e1) * t151;
t200 = -t229 + t230;
t199 = rSges(4,1) * t151 - rSges(4,2) * t150;
t198 = -t115 * rSges(6,1) - t114 * rSges(6,2);
t20 = -t233 * t150 + t231 * t219;
t21 = t234 * t150 - t231 * t218;
t192 = t165 * t21 + t168 * t20;
t177 = -pkin(8) * t150 + t207;
t169 = t177 - t231;
t39 = t169 * t165;
t40 = t169 * t168;
t191 = t165 * t39 + t168 * t40;
t190 = Icges(3,1) * t167 - t227;
t189 = Icges(4,1) * t151 - t225;
t188 = -Icges(3,2) * t164 + t226;
t187 = -Icges(4,2) * t150 + t224;
t183 = -Icges(5,2) * t151 + t223;
t182 = Icges(5,3) * t150 - t222;
t179 = rSges(4,1) * t218 - rSges(4,2) * t220 + t165 * rSges(4,3);
t178 = t165 * rSges(5,1) - rSges(5,2) * t218 + rSges(5,3) * t220;
t175 = t205 + t212;
t174 = t27 / 0.2e1 + t24 / 0.2e1 + t22 / 0.2e1 + t28 / 0.2e1;
t173 = t30 / 0.2e1 + t29 / 0.2e1 + t25 / 0.2e1 + t23 / 0.2e1;
t172 = t165 * (pkin(8) * t219 - t157) + t168 * t210 + t204;
t31 = (t146 - t162) * t168 + (-t147 + (-qJ(4) - t237) * t150 + (-pkin(3) - t262) * t151) * t165 + t197;
t32 = t175 + t246;
t171 = m(7) * (t165 * t32 + t168 * t31);
t83 = rSges(6,3) * t150 + (-rSges(6,1) * t163 - rSges(6,2) * t166) * t151;
t170 = t177 - t83;
t149 = t151 ^ 2;
t148 = t150 ^ 2;
t134 = rSges(2,1) * t168 - t165 * rSges(2,2);
t133 = -t165 * rSges(2,1) - rSges(2,2) * t168;
t132 = rSges(3,1) * t164 + rSges(3,2) * t167;
t85 = t206 * t168;
t84 = t206 * t165;
t75 = t165 * pkin(7) + (pkin(1) - t229) * t168 + t211;
t74 = t228 + t156 + (-pkin(1) - t200) * t165;
t69 = t179 + t205;
t68 = (rSges(4,3) - t162) * t168 + (-t147 - t199) * t165;
t67 = t202 * t168;
t66 = t202 * t165;
t65 = t168 * (-t168 * t229 + t211) + (t200 * t165 - t228) * t165;
t62 = rSges(6,3) * t219 - t198;
t46 = t175 + t178;
t45 = (rSges(5,1) - t162) * t168 + (-t147 + (rSges(5,2) - pkin(3)) * t151 + (-rSges(5,3) - qJ(4)) * t150) * t165;
t44 = t170 * t168;
t43 = t170 * t165;
t42 = t150 * t60 - t83 * t218;
t41 = -t150 * t62 + t83 * t219;
t38 = t175 + t60 + t210;
t37 = -t217 + t157 + (-t221 - t147 + (-rSges(6,3) - pkin(3) - pkin(8)) * t151) * t165 + t198;
t34 = t168 * t179 + (-rSges(4,3) * t168 + t199 * t165) * t165 + t232;
t33 = (-t165 * t60 + t168 * t62) * t151;
t26 = t168 * t178 + (-rSges(5,1) * t168 + (-rSges(5,2) * t151 + rSges(5,3) * t150) * t165) * t165 + t204;
t11 = (-t234 * t165 + t233 * t168) * t151;
t10 = t165 * t62 + t168 * t60 + t172;
t9 = t233 * t165 + t234 * t168 + t172;
t8 = t18 * t165 - t168 * t19;
t7 = t16 * t165 - t168 * t17;
t6 = t14 * t165 - t15 * t168;
t5 = t12 * t165 - t13 * t168;
t1 = [t167 * (Icges(3,2) * t167 + t227) + t164 * (Icges(3,1) * t164 + t226) + Icges(2,3) + m(7) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t68 ^ 2 + t69 ^ 2) + m(5) * (t45 ^ 2 + t46 ^ 2) + m(3) * (t74 ^ 2 + t75 ^ 2) + m(2) * (t133 ^ 2 + t134 ^ 2) + (t223 + t225 + (Icges(5,3) + Icges(4,2)) * t151 + t254) * t151 + (t222 + t224 + (Icges(4,1) + Icges(5,2)) * t150) * t150 + t255; (-(-Icges(3,6) * t168 + t188 * t165) * t167 / 0.2e1 - (-Icges(3,5) * t168 + t190 * t165) * t164 / 0.2e1 + t201 * t168 + (Icges(5,5) * t257 + Icges(4,6) * t256 + t182 * t258 + t187 * t259) * t151 + (Icges(5,4) * t257 + Icges(4,5) * t256 + t183 * t258 + t189 * t259) * t150 - t173) * t168 + ((Icges(3,6) * t165 + t188 * t168) * t249 + (Icges(3,5) * t165 + t190 * t168) * t250 + t201 * t165 + (Icges(5,5) * t259 + Icges(4,6) * t258 + t182 * t257 + t187 * t256) * t151 + (Icges(5,4) * t259 + Icges(4,5) * t258 + t183 * t257 + t189 * t256) * t150 + t174) * t165 + m(4) * (t68 * t85 + t69 * t84) + m(5) * (t45 * t67 + t46 * t66) + m(6) * (t37 * t44 + t38 * t43) + m(7) * (t31 * t40 + t32 * t39) + m(3) * (-t165 * t75 - t168 * t74) * t132; m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t26 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t34 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(3) * (t209 * t132 ^ 2 + t65 ^ 2) + (t245 * t160 - t7 - t8) * t168 + (t5 + t6 + t244 * t159 + (t245 * t165 + t244 * t168) * t168) * t165; m(7) * (t165 * t31 - t168 * t32) + m(6) * (t165 * t37 - t168 * t38) + m(4) * (t165 * t68 - t168 * t69) + m(5) * (t165 * t45 - t168 * t46); m(7) * (t165 * t40 - t168 * t39) + m(6) * (t165 * t44 - t168 * t43) + m(5) * (t165 * t67 - t168 * t66) + m(4) * (t165 * t85 - t168 * t84); 0.2e1 * (m(4) / 0.2e1 + t203) * t209; 0.2e1 * (t171 / 0.2e1 + (t165 * t38 + t168 * t37) * t242 + (t165 * t46 + t168 * t45) * t243) * t150; m(7) * (t191 * t150 - t151 * t9) + m(6) * (-t151 * t10 + (t165 * t43 + t168 * t44) * t150) + m(5) * (-t151 * t26 + (t165 * t66 + t168 * t67) * t150); 0; 0.2e1 * t203 * (t209 * t148 + t149); m(7) * (t20 * t31 + t21 * t32) + m(6) * (t37 * t41 + t38 * t42) + (t173 * t165 + t174 * t168) * t151 + t235; m(7) * (t11 * t9 + t20 * t40 + t21 * t39) + m(6) * (t10 * t33 + t41 * t44 + t42 * t43) + ((t5 / 0.2e1 + t6 / 0.2e1) * t168 + (t8 / 0.2e1 + t7 / 0.2e1) * t165) * t151 + (t248 * t165 - t247 * t168) * t251 + t253 * t258 + t252 * t257; m(6) * (t41 * t165 - t168 * t42) + m(7) * (t20 * t165 - t168 * t21); m(6) * (-t33 * t151 + (t165 * t42 + t168 * t41) * t150) + m(7) * (-t11 * t151 + t192 * t150); t235 * t150 + m(7) * (t11 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t33 ^ 2 + t41 ^ 2 + t42 ^ 2) + (t253 * t168 + t252 * t165 + (t247 * t165 + t248 * t168) * t150) * t151; t151 * t171; m(7) * (t150 * t9 + t191 * t151); 0; m(7) * (-0.1e1 + t209) * t151 * t150; m(7) * (t150 * t11 + t192 * t151); m(7) * (t209 * t149 + t148);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
