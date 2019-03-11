% Calculate joint inertia matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:14:03
% EndTime: 2019-03-09 08:14:08
% DurationCPUTime: 3.18s
% Computational Cost: add. (2685->356), mult. (4539->516), div. (0->0), fcn. (4578->8), ass. (0->166)
t249 = Icges(3,1) + Icges(4,1);
t248 = Icges(5,1) + Icges(4,3);
t243 = Icges(5,5) + Icges(3,6);
t152 = cos(qJ(2));
t247 = (Icges(5,4) - Icges(4,5)) * t152;
t150 = sin(qJ(2));
t246 = (Icges(3,4) - Icges(4,5)) * t150;
t241 = Icges(5,6) + Icges(4,4) + Icges(3,5);
t245 = t150 * t248 - t247;
t244 = t152 * t249 - t246;
t242 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t240 = -t241 * t152 + (-Icges(4,6) + t243) * t150;
t239 = t241 * t150;
t151 = sin(qJ(1));
t238 = -t151 / 0.2e1;
t224 = t151 / 0.2e1;
t153 = cos(qJ(1));
t237 = -t153 / 0.2e1;
t236 = t153 / 0.2e1;
t235 = t150 / 0.2e1;
t148 = cos(pkin(9));
t134 = pkin(5) * t148 + pkin(4);
t149 = -pkin(8) - qJ(5);
t202 = t152 * t153;
t209 = t150 * t153;
t142 = pkin(9) + qJ(6);
t135 = sin(t142);
t136 = cos(t142);
t206 = t151 * t136;
t70 = -t135 * t209 - t206;
t207 = t151 * t135;
t71 = t136 * t209 - t207;
t38 = rSges(7,1) * t71 + rSges(7,2) * t70 + rSges(7,3) * t202;
t233 = t134 * t209 - t149 * t202 + t38;
t232 = t151 * t240 + t153 * t242;
t231 = t151 * t242 - t153 * t240;
t144 = t151 ^ 2;
t146 = t153 ^ 2;
t230 = 0.2e1 * t152;
t229 = m(4) / 0.2e1;
t228 = m(5) / 0.2e1;
t227 = m(6) / 0.2e1;
t226 = m(7) / 0.2e1;
t225 = -pkin(2) - pkin(3);
t117 = rSges(3,1) * t150 + rSges(3,2) * t152;
t223 = m(3) * t117;
t57 = Icges(7,3) * t150 + (-Icges(7,5) * t136 + Icges(7,6) * t135) * t152;
t58 = Icges(7,6) * t150 + (-Icges(7,4) * t136 + Icges(7,2) * t135) * t152;
t222 = t135 * t152 * t58 + t150 * t57;
t199 = pkin(2) * t202 + qJ(3) * t209;
t221 = t144 * (pkin(2) * t152 + qJ(3) * t150) + t153 * t199;
t59 = Icges(7,5) * t150 + (-Icges(7,1) * t136 + Icges(7,4) * t135) * t152;
t220 = t136 * t59;
t219 = t153 * rSges(4,2);
t218 = t153 * rSges(3,3);
t217 = -rSges(5,3) - qJ(4);
t215 = Icges(3,4) * t152;
t214 = Icges(5,4) * t150;
t210 = qJ(4) * t153;
t208 = t151 * qJ(4);
t147 = sin(pkin(9));
t205 = t151 * t147;
t204 = t151 * t148;
t203 = t151 * t152;
t115 = pkin(2) * t150 - qJ(3) * t152;
t201 = -rSges(4,1) * t150 + rSges(4,3) * t152 - t115;
t200 = pkin(4) * t209 + qJ(5) * t202;
t198 = pkin(1) * t153 + pkin(7) * t151;
t197 = t144 + t146;
t196 = t227 + t226;
t100 = -t147 * t209 - t204;
t101 = t148 * t209 - t205;
t195 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t202;
t68 = t136 * t153 - t150 * t207;
t69 = t135 * t153 + t150 * t206;
t31 = Icges(7,5) * t69 + Icges(7,6) * t68 + Icges(7,3) * t203;
t33 = Icges(7,4) * t69 + Icges(7,2) * t68 + Icges(7,6) * t203;
t35 = Icges(7,1) * t69 + Icges(7,4) * t68 + Icges(7,5) * t203;
t11 = t150 * t31 + (t135 * t33 - t136 * t35) * t152;
t13 = t203 * t57 + t58 * t68 + t59 * t69;
t194 = t11 / 0.2e1 + t13 / 0.2e1;
t32 = Icges(7,5) * t71 + Icges(7,6) * t70 + Icges(7,3) * t202;
t34 = Icges(7,4) * t71 + Icges(7,2) * t70 + Icges(7,6) * t202;
t36 = Icges(7,1) * t71 + Icges(7,4) * t70 + Icges(7,5) * t202;
t12 = t150 * t32 + (t135 * t34 - t136 * t36) * t152;
t14 = t202 * t57 + t58 * t70 + t59 * t71;
t193 = t12 / 0.2e1 + t14 / 0.2e1;
t192 = rSges(4,1) * t202 + rSges(4,2) * t151 + rSges(4,3) * t209;
t191 = -pkin(3) * t150 - t115;
t190 = -pkin(5) * t147 - qJ(4);
t132 = pkin(3) * t202;
t189 = t151 * (pkin(3) * t203 + t210) + t153 * (t132 - t208) + t221;
t188 = t228 + t196;
t187 = t198 + t199;
t186 = pkin(4) * t152 - qJ(5) * t150 + t191;
t185 = rSges(5,1) * t152 + rSges(5,2) * t150 + t191;
t184 = rSges(5,1) * t209 - rSges(5,2) * t202;
t98 = t148 * t153 - t150 * t205;
t99 = t147 * t153 + t150 * t204;
t183 = -t99 * rSges(6,1) - t98 * rSges(6,2);
t182 = -t69 * rSges(7,1) - t68 * rSges(7,2);
t180 = t241 * t235 + (-Icges(4,6) / 0.2e1 + t243 / 0.2e1) * t152;
t179 = rSges(3,1) * t152 - rSges(3,2) * t150;
t178 = pkin(4) * t150 + qJ(5) * t152;
t61 = rSges(7,3) * t150 + (-rSges(7,1) * t136 + rSges(7,2) * t135) * t152;
t155 = t186 - (pkin(4) - t134) * t152 - (-qJ(5) - t149) * t150 - t61;
t22 = t155 * t151;
t23 = t155 * t153;
t171 = t151 * t22 + t153 * t23;
t37 = rSges(7,3) * t203 - t182;
t24 = -t150 * t37 + t203 * t61;
t25 = t150 * t38 - t202 * t61;
t170 = t151 * t25 + t153 * t24;
t158 = t186 - rSges(6,3) * t150 - (-rSges(6,1) * t148 + rSges(6,2) * t147) * t152;
t27 = t158 * t151;
t28 = t158 * t153;
t169 = t151 * t27 + t153 * t28;
t168 = t132 + t187;
t164 = -Icges(3,2) * t150 + t215;
t162 = -Icges(5,2) * t152 + t214;
t157 = t144 * t178 + t153 * t200 + t189;
t156 = rSges(3,1) * t202 - rSges(3,2) * t209 + rSges(3,3) * t151;
t140 = t153 * pkin(7);
t16 = t140 + t190 * t153 + (-pkin(1) + (-qJ(3) - t134) * t150 + (-rSges(7,3) + t149 + t225) * t152) * t151 + t182;
t17 = t151 * t190 + t168 + t233;
t20 = -t210 + t140 + (-pkin(1) + (-pkin(4) - qJ(3)) * t150 + (-rSges(6,3) - qJ(5) + t225) * t152) * t151 + t183;
t21 = t168 + t195 + t200 - t208;
t154 = (t151 * t17 + t153 * t16) * t226 + (t151 * t21 + t153 * t20) * t227;
t145 = t152 ^ 2;
t143 = t150 ^ 2;
t121 = rSges(2,1) * t153 - rSges(2,2) * t151;
t118 = -rSges(2,1) * t151 - rSges(2,2) * t153;
t66 = Icges(6,5) * t150 + (-Icges(6,1) * t148 + Icges(6,4) * t147) * t152;
t65 = Icges(6,6) * t150 + (-Icges(6,4) * t148 + Icges(6,2) * t147) * t152;
t55 = t201 * t153;
t54 = t201 * t151;
t52 = t156 + t198;
t51 = t218 + t140 + (-pkin(1) - t179) * t151;
t49 = t185 * t153;
t48 = t185 * t151;
t47 = t187 + t192;
t46 = t219 + t140 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t152 + (-rSges(4,3) - qJ(3)) * t150) * t151;
t45 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t202;
t44 = Icges(6,1) * t99 + Icges(6,4) * t98 + Icges(6,5) * t203;
t43 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t202;
t42 = Icges(6,4) * t99 + Icges(6,2) * t98 + Icges(6,6) * t203;
t41 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t202;
t40 = Icges(6,5) * t99 + Icges(6,6) * t98 + Icges(6,3) * t203;
t39 = t153 * t156 + (t151 * t179 - t218) * t151;
t30 = t151 * t217 + t168 + t184;
t29 = t140 + t217 * t153 + (-pkin(1) + (-rSges(5,1) - qJ(3)) * t150 + (rSges(5,2) + t225) * t152) * t151;
t26 = t153 * t192 + (-t219 + (rSges(4,1) * t152 + rSges(4,3) * t150) * t151) * t151 + t221;
t19 = (-t152 * t220 + t222) * t150;
t18 = t153 * t184 + (rSges(5,1) * t150 - rSges(5,2) * t152) * t144 + t189;
t15 = (-t151 * t38 + t153 * t37) * t152;
t10 = t151 * (rSges(6,3) * t203 - t183) + t153 * t195 + t157;
t9 = t202 * t32 + t34 * t70 + t36 * t71;
t8 = t202 * t31 + t33 * t70 + t35 * t71;
t7 = t203 * t32 + t34 * t68 + t36 * t69;
t6 = t203 * t31 + t33 * t68 + t35 * t69;
t5 = (-t200 + t233) * t153 + (t37 + (t134 * t150 - t149 * t152 - t178) * t151) * t151 + t157;
t4 = t151 * t9 - t153 * t8;
t3 = t151 * t7 - t153 * t6;
t2 = t14 * t150 + (t151 * t8 + t153 * t9) * t152;
t1 = t13 * t150 + (t151 * t6 + t153 * t7) * t152;
t50 = [Icges(2,3) + m(7) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(3) * (t51 ^ 2 + t52 ^ 2) + m(2) * (t118 ^ 2 + t121 ^ 2) + t222 + (t147 * t65 - t148 * t66 + t214 - t220 + (Icges(3,2) + t248) * t152 + t246) * t152 + ((-Icges(6,5) * t148 + Icges(6,6) * t147) * t152 + t215 + (Icges(5,2) + Icges(6,3) + t249) * t150 + t247) * t150; m(4) * (t46 * t55 + t47 * t54) + m(5) * (t29 * t49 + t30 * t48) + m(6) * (t20 * t28 + t21 * t27) + m(7) * (t16 * t23 + t17 * t22) + (-t51 * t223 - t65 * t98 / 0.2e1 - t66 * t99 / 0.2e1 + t180 * t153 + (-t40 / 0.2e1 + t162 * t224 + t244 * t238) * t150 - t194 + t236 * t239) * t153 + (-t52 * t223 + t100 * t65 / 0.2e1 + t101 * t66 / 0.2e1 + t180 * t151 + (t41 / 0.2e1 + t162 * t237 + t244 * t236) * t150 + t193 + t224 * t239) * t151 + ((-t147 * t42 / 0.2e1 + t148 * t44 / 0.2e1 + Icges(4,6) * t237 + t164 * t238 + t243 * t236) * t153 + (t147 * t43 / 0.2e1 - t148 * t45 / 0.2e1 + Icges(4,6) * t238 + t164 * t236 + t245 * t237) * t151 + (t243 * t151 + t153 * t245) * t224) * t152; m(7) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(6) * (t10 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(5) * (t18 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t26 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(3) * (t117 ^ 2 * t197 + t39 ^ 2) + (-t3 + (t40 * t203 + t98 * t42 + t99 * t44) * t153 + t232 * t146) * t153 + (t4 + (t100 * t43 + t101 * t45 + t41 * t202) * t151 + t231 * t144 + (-t100 * t42 - t101 * t44 + t151 * t232 + t153 * t231 - t40 * t202 - t41 * t203 - t43 * t98 - t45 * t99) * t153) * t151; 0.2e1 * ((t151 * t47 + t153 * t46) * t229 + (t151 * t30 + t153 * t29) * t228 + t154) * t150; m(7) * (t150 * t171 - t152 * t5) + m(6) * (-t10 * t152 + t150 * t169) + m(5) * (-t152 * t18 + (t151 * t48 + t153 * t49) * t150) + m(4) * (-t152 * t26 + (t151 * t54 + t153 * t55) * t150); 0.2e1 * (t229 + t188) * (t143 * t197 + t145); m(7) * (-t151 * t16 + t153 * t17) + m(6) * (-t151 * t20 + t153 * t21) + m(5) * (-t151 * t29 + t153 * t30); m(7) * (-t151 * t23 + t153 * t22) + m(6) * (-t151 * t28 + t153 * t27) + m(5) * (-t151 * t49 + t153 * t48); 0; 0.2e1 * t188 * t197; t154 * t230; m(7) * (t150 * t5 + t152 * t171) + m(6) * (t150 * t10 + t152 * t169); t196 * (-0.1e1 + t197) * t150 * t230; 0; 0.2e1 * t196 * (t145 * t197 + t143); m(7) * (t16 * t24 + t17 * t25) + t19 + (t151 * t194 + t153 * t193) * t152; t2 * t224 + t1 * t237 + m(7) * (t15 * t5 + t22 * t25 + t23 * t24) + (-t11 * t153 + t12 * t151) * t235 + (t3 * t224 + t236 * t4) * t152; m(7) * (-t15 * t152 + t150 * t170); m(7) * (-t151 * t24 + t153 * t25); m(7) * (t15 * t150 + t152 * t170); t150 * t19 + m(7) * (t15 ^ 2 + t24 ^ 2 + t25 ^ 2) + (t153 * t2 + t151 * t1 + t150 * (t11 * t151 + t12 * t153)) * t152;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t50(1) t50(2) t50(4) t50(7) t50(11) t50(16); t50(2) t50(3) t50(5) t50(8) t50(12) t50(17); t50(4) t50(5) t50(6) t50(9) t50(13) t50(18); t50(7) t50(8) t50(9) t50(10) t50(14) t50(19); t50(11) t50(12) t50(13) t50(14) t50(15) t50(20); t50(16) t50(17) t50(18) t50(19) t50(20) t50(21);];
Mq  = res;
