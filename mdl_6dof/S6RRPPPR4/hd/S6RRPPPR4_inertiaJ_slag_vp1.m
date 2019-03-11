% Calculate joint inertia matrix for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:19
% EndTime: 2019-03-09 08:17:26
% DurationCPUTime: 3.06s
% Computational Cost: add. (2914->390), mult. (7104->562), div. (0->0), fcn. (8026->8), ass. (0->181)
t247 = Icges(4,1) + Icges(3,3);
t154 = sin(qJ(2));
t157 = cos(qJ(2));
t246 = (-Icges(4,4) + Icges(3,5)) * t157 + (Icges(4,5) - Icges(3,6)) * t154;
t155 = sin(qJ(1));
t245 = -t155 / 0.2e1;
t244 = t155 / 0.2e1;
t158 = cos(qJ(1));
t227 = -t158 / 0.2e1;
t243 = t158 / 0.2e1;
t152 = cos(pkin(9));
t212 = t154 * t158;
t151 = sin(pkin(9));
t214 = t151 * t155;
t111 = -t152 * t212 + t214;
t211 = t155 * t152;
t112 = t151 * t212 + t211;
t209 = t157 * t158;
t50 = Icges(6,5) * t112 + Icges(6,6) * t209 + Icges(6,3) * t111;
t56 = Icges(5,4) * t112 - Icges(5,2) * t111 + Icges(5,6) * t209;
t242 = t50 - t56;
t113 = t151 * t158 + t154 * t211;
t115 = -t152 * t158 + t154 * t214;
t210 = t155 * t157;
t51 = Icges(6,5) * t115 + Icges(6,6) * t210 - Icges(6,3) * t113;
t57 = Icges(5,4) * t115 + Icges(5,2) * t113 + Icges(5,6) * t210;
t241 = -t51 + t57;
t52 = Icges(5,5) * t112 - Icges(5,6) * t111 + Icges(5,3) * t209;
t54 = Icges(6,4) * t112 + Icges(6,2) * t209 + Icges(6,6) * t111;
t240 = t54 + t52;
t53 = Icges(5,5) * t115 + Icges(5,6) * t113 + Icges(5,3) * t210;
t55 = Icges(6,4) * t115 + Icges(6,2) * t210 - Icges(6,6) * t113;
t239 = -t55 - t53;
t59 = Icges(6,1) * t115 + Icges(6,4) * t210 - Icges(6,5) * t113;
t61 = Icges(5,1) * t115 + Icges(5,4) * t113 + Icges(5,5) * t210;
t238 = t59 + t61;
t58 = Icges(6,1) * t112 + Icges(6,4) * t209 + Icges(6,5) * t111;
t60 = Icges(5,1) * t112 - Icges(5,4) * t111 + Icges(5,5) * t209;
t237 = -t60 - t58;
t228 = m(7) / 0.2e1;
t229 = m(6) / 0.2e1;
t201 = t229 + t228;
t236 = 0.2e1 * t201;
t235 = -t154 / 0.2e1;
t234 = -t246 * t155 + t247 * t158;
t233 = t247 * t155 + t246 * t158;
t147 = t154 ^ 2;
t148 = t155 ^ 2;
t150 = t158 ^ 2;
t232 = 0.2e1 * t157;
t231 = m(4) / 0.2e1;
t230 = m(5) / 0.2e1;
t226 = -rSges(7,3) - pkin(8);
t225 = -pkin(2) - qJ(4);
t153 = sin(qJ(6));
t156 = cos(qJ(6));
t98 = (t151 * t153 + t152 * t156) * t157;
t99 = (-t151 * t156 + t152 * t153) * t157;
t44 = Icges(7,4) * t99 + Icges(7,2) * t98 - Icges(7,6) * t154;
t45 = Icges(7,1) * t99 + Icges(7,4) * t98 - Icges(7,5) * t154;
t224 = t98 * t44 + t99 * t45;
t69 = t111 * t156 - t112 * t153;
t70 = t111 * t153 + t112 * t156;
t223 = t70 * rSges(7,1) + t69 * rSges(7,2);
t222 = rSges(6,3) * t113;
t221 = t158 * rSges(4,1);
t220 = t158 * rSges(3,3);
t219 = Icges(3,4) * t154;
t218 = Icges(3,4) * t157;
t217 = Icges(4,6) * t154;
t216 = Icges(4,6) * t157;
t215 = qJ(3) * t154;
t213 = t152 * t157;
t206 = pkin(2) * t209 + qJ(3) * t212;
t208 = t148 * (pkin(2) * t157 + t215) + t158 * t206;
t127 = pkin(2) * t154 - qJ(3) * t157;
t207 = rSges(4,2) * t154 + rSges(4,3) * t157 - t127;
t205 = t158 * pkin(1) + t155 * pkin(7);
t204 = t155 * pkin(3) + qJ(4) * t209;
t144 = t158 * pkin(7);
t145 = t158 * pkin(3);
t203 = t144 + t145;
t202 = t148 + t150;
t26 = Icges(7,5) * t70 + Icges(7,6) * t69 - Icges(7,3) * t209;
t28 = Icges(7,4) * t70 + Icges(7,2) * t69 - Icges(7,6) * t209;
t30 = Icges(7,1) * t70 + Icges(7,4) * t69 - Icges(7,5) * t209;
t10 = -t154 * t26 + t28 * t98 + t30 * t99;
t43 = Icges(7,5) * t99 + Icges(7,6) * t98 - Icges(7,3) * t154;
t13 = -t43 * t209 + t69 * t44 + t70 * t45;
t200 = -t10 / 0.2e1 - t13 / 0.2e1;
t71 = -t113 * t156 - t115 * t153;
t72 = -t113 * t153 + t115 * t156;
t27 = Icges(7,5) * t72 + Icges(7,6) * t71 - Icges(7,3) * t210;
t29 = Icges(7,4) * t72 + Icges(7,2) * t71 - Icges(7,6) * t210;
t31 = Icges(7,1) * t72 + Icges(7,4) * t71 - Icges(7,5) * t210;
t11 = -t154 * t27 + t29 * t98 + t31 * t99;
t14 = -t43 * t210 + t44 * t71 + t45 * t72;
t199 = -t11 / 0.2e1 - t14 / 0.2e1;
t78 = Icges(6,6) * t154 + (-Icges(6,5) * t151 + Icges(6,3) * t152) * t157;
t81 = Icges(5,6) * t154 + (-Icges(5,4) * t151 - Icges(5,2) * t152) * t157;
t198 = -t81 / 0.2e1 + t78 / 0.2e1;
t82 = Icges(6,4) * t154 + (-Icges(6,1) * t151 + Icges(6,5) * t152) * t157;
t83 = Icges(5,5) * t154 + (-Icges(5,1) * t151 - Icges(5,4) * t152) * t157;
t197 = -t83 / 0.2e1 - t82 / 0.2e1;
t102 = t113 * qJ(5);
t196 = t102 + t203;
t195 = t112 * rSges(6,1) + rSges(6,2) * t209 + t111 * rSges(6,3);
t194 = t112 * rSges(5,1) - t111 * rSges(5,2) + rSges(5,3) * t209;
t193 = Icges(4,4) * t235 + Icges(3,5) * t154 / 0.2e1 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t157;
t192 = -pkin(1) - t215;
t191 = -qJ(4) * t154 - t127;
t190 = t112 * pkin(4) + t111 * qJ(5);
t189 = t230 + t201;
t188 = t155 * (qJ(4) * t210 - t145) + t158 * t204 + t208;
t187 = t205 + t206;
t186 = t191 - t154 * rSges(5,3) - (-rSges(5,1) * t151 - rSges(5,2) * t152) * t157;
t185 = -rSges(7,1) * t72 - rSges(7,2) * t71;
t184 = -(-pkin(4) * t151 + qJ(5) * t152) * t157 + t191;
t183 = rSges(3,1) * t157 - rSges(3,2) * t154;
t182 = -t115 * rSges(5,1) - t113 * rSges(5,2);
t33 = -rSges(7,3) * t210 - t185;
t46 = rSges(7,1) * t99 + rSges(7,2) * t98 - rSges(7,3) * t154;
t20 = t154 * t33 - t46 * t210;
t32 = -rSges(7,3) * t209 + t223;
t21 = -t154 * t32 + t46 * t209;
t177 = t155 * t21 + t158 * t20;
t161 = pkin(5) * t151 * t157 + pkin(8) * t154 + t184 - t46;
t24 = t161 * t155;
t25 = t161 * t158;
t176 = t155 * t24 + t158 * t25;
t164 = t184 - t154 * rSges(6,2) - (-rSges(6,1) * t151 + rSges(6,3) * t152) * t157;
t39 = t164 * t155;
t40 = t164 * t158;
t175 = t155 * t39 + t158 * t40;
t48 = t186 * t155;
t49 = t186 * t158;
t174 = t155 * t48 + t158 * t49;
t173 = Icges(3,1) * t157 - t219;
t172 = -Icges(3,2) * t154 + t218;
t169 = -Icges(4,2) * t157 + t217;
t168 = Icges(4,3) * t154 - t216;
t167 = t111 * t158 - t113 * t155;
t166 = rSges(3,1) * t209 - rSges(3,2) * t212 + t155 * rSges(3,3);
t165 = t155 * rSges(4,1) - rSges(4,2) * t209 + rSges(4,3) * t212;
t163 = t155 * (pkin(4) * t115 - t102) + t158 * t190 + t188;
t162 = t187 + t204;
t160 = t162 + t190;
t16 = (-pkin(4) - pkin(5)) * t115 + ((t225 - t226) * t157 + t192) * t155 + t185 + t196;
t107 = t112 * pkin(5);
t17 = t226 * t209 + t107 + t160 + t223;
t22 = t222 + (-rSges(6,1) - pkin(4)) * t115 + ((-rSges(6,2) + t225) * t157 + t192) * t155 + t196;
t23 = t160 + t195;
t35 = ((-rSges(5,3) + t225) * t157 + t192) * t155 + t182 + t203;
t36 = t162 + t194;
t159 = (t155 * t17 + t158 * t16) * t228 + (t155 * t23 + t158 * t22) * t229 + (t155 * t36 + t158 * t35) * t230;
t149 = t157 ^ 2;
t131 = rSges(2,1) * t158 - t155 * rSges(2,2);
t130 = -t155 * rSges(2,1) - rSges(2,2) * t158;
t129 = rSges(3,1) * t154 + rSges(3,2) * t157;
t76 = t207 * t158;
t75 = t207 * t155;
t74 = t166 + t205;
t73 = t220 + t144 + (-pkin(1) - t183) * t155;
t63 = t165 + t187;
t62 = t221 + t144 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t157 + (-rSges(4,3) - qJ(3)) * t154) * t155;
t47 = t158 * t166 + (t183 * t155 - t220) * t155;
t34 = t158 * t165 + (-t221 + (-rSges(4,2) * t157 + rSges(4,3) * t154) * t155) * t155 + t208;
t19 = t155 * (rSges(5,3) * t210 - t182) + t158 * t194 + t188;
t18 = -t154 * t43 + t224;
t15 = (t155 * t32 - t158 * t33) * t157;
t12 = t155 * (rSges(6,1) * t115 + rSges(6,2) * t210 - t222) + t158 * t195 + t163;
t9 = -t27 * t210 + t29 * t71 + t31 * t72;
t8 = -t26 * t210 + t28 * t71 + t30 * t72;
t7 = -t27 * t209 + t69 * t29 + t70 * t31;
t6 = -t26 * t209 + t69 * t28 + t70 * t30;
t5 = (-pkin(8) * t209 + t107 + t32) * t158 + (pkin(5) * t115 - pkin(8) * t210 + t33) * t155 + t163;
t4 = t8 * t155 - t158 * t9;
t3 = t6 * t155 - t158 * t7;
t2 = -t14 * t154 + (-t155 * t9 - t158 * t8) * t157;
t1 = -t13 * t154 + (-t155 * t7 - t158 * t6) * t157;
t37 = [Icges(2,3) + (-t43 + t218 + t216 + ((-Icges(5,6) + Icges(6,6)) * t152 + (-Icges(6,4) - Icges(5,5)) * t151) * t157 + (Icges(3,1) + Icges(4,2) + Icges(5,3) + Icges(6,2)) * t154) * t154 + (t217 + t219 + (t78 - t81) * t152 + (-t82 - t83) * t151 + (Icges(4,3) + Icges(3,2)) * t157) * t157 + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2) + m(3) * (t73 ^ 2 + t74 ^ 2) + m(2) * (t130 ^ 2 + t131 ^ 2) + t224; (t198 * t113 + t197 * t115 + t193 * t158 + t199) * t158 + (t198 * t111 - t197 * t112 + t193 * t155 - t200) * t155 + m(4) * (t62 * t76 + t63 * t75) + m(5) * (t35 * t49 + t36 * t48) + m(6) * (t22 * t40 + t23 * t39) + m(7) * (t16 * t25 + t17 * t24) + m(3) * (-t155 * t74 - t158 * t73) * t129 + ((Icges(3,5) * t243 + t173 * t245 + Icges(4,4) * t227 + t169 * t244 - t55 / 0.2e1 - t53 / 0.2e1) * t158 + (Icges(3,5) * t244 + t173 * t243 + Icges(4,4) * t245 + t169 * t227 + t54 / 0.2e1 + t52 / 0.2e1) * t155) * t154 + ((Icges(3,6) * t243 + t172 * t245 + Icges(4,5) * t227 + t168 * t244 + (-t51 / 0.2e1 + t57 / 0.2e1) * t152 + (t59 / 0.2e1 + t61 / 0.2e1) * t151) * t158 + (Icges(3,6) * t244 + t172 * t243 + Icges(4,5) * t245 + t168 * t227 + (t50 / 0.2e1 - t56 / 0.2e1) * t152 + (-t58 / 0.2e1 - t60 / 0.2e1) * t151) * t155) * t157; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t19 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t34 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(3) * (t202 * t129 ^ 2 + t47 ^ 2) + (-t4 + (t241 * t113 + t238 * t115 - t239 * t210) * t158 + t234 * t150) * t158 + (t3 + (t242 * t111 - t237 * t112 + t240 * t209) * t155 + t233 * t148 + (t241 * t111 - t238 * t112 + t242 * t113 + t237 * t115 + t234 * t155 + t233 * t158 + t239 * t209 - t240 * t210) * t158) * t155; 0.2e1 * ((t155 * t63 + t158 * t62) * t231 + t159) * t154; m(7) * (t176 * t154 - t157 * t5) + m(6) * (-t157 * t12 + t175 * t154) + m(5) * (t174 * t154 - t157 * t19) + m(4) * (-t157 * t34 + (t155 * t75 + t158 * t76) * t154); 0.2e1 * (t231 + t189) * (t202 * t147 + t149); t159 * t232; m(7) * (t154 * t5 + t176 * t157) + m(6) * (t154 * t12 + t175 * t157) + m(5) * (t154 * t19 + t174 * t157); t189 * (-0.1e1 + t202) * t154 * t232; 0.2e1 * t189 * (t202 * t149 + t147); m(7) * (t111 * t16 - t113 * t17) + m(6) * (t111 * t22 - t113 * t23); m(7) * (t111 * t25 - t113 * t24 + t5 * t213) + m(6) * (t111 * t40 - t113 * t39 + t12 * t213); (-t149 * t152 + t167 * t154) * t236; t201 * (t152 * t154 + t167) * t232; (t149 * t152 ^ 2 + t111 ^ 2 + t113 ^ 2) * t236; m(7) * (t16 * t20 + t17 * t21) - t18 * t154 + (t199 * t155 + t200 * t158) * t157; m(7) * (t15 * t5 + t20 * t25 + t21 * t24) + t2 * t227 + t1 * t244 + (t10 * t155 - t11 * t158) * t235 + (t3 * t227 + t4 * t245) * t157; m(7) * (-t15 * t157 + t177 * t154); m(7) * (t15 * t154 + t177 * t157); m(7) * (t111 * t20 - t113 * t21 + t15 * t213); t147 * t18 + m(7) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + (-t158 * t1 - t155 * t2 - t154 * (-t10 * t158 - t11 * t155)) * t157;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t37(1) t37(2) t37(4) t37(7) t37(11) t37(16); t37(2) t37(3) t37(5) t37(8) t37(12) t37(17); t37(4) t37(5) t37(6) t37(9) t37(13) t37(18); t37(7) t37(8) t37(9) t37(10) t37(14) t37(19); t37(11) t37(12) t37(13) t37(14) t37(15) t37(20); t37(16) t37(17) t37(18) t37(19) t37(20) t37(21);];
Mq  = res;
