% Calculate joint inertia matrix for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:58
% EndTime: 2019-03-09 08:10:04
% DurationCPUTime: 2.86s
% Computational Cost: add. (4102->372), mult. (4393->548), div. (0->0), fcn. (4436->10), ass. (0->174)
t245 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t146 = qJ(2) + pkin(9);
t136 = sin(t146);
t138 = cos(t146);
t153 = sin(qJ(2));
t155 = cos(qJ(2));
t244 = Icges(3,5) * t155 - Icges(3,6) * t153 + (-Icges(5,4) + Icges(4,5)) * t138 + (Icges(5,5) - Icges(4,6)) * t136;
t154 = sin(qJ(1));
t243 = -t154 / 0.2e1;
t229 = t154 / 0.2e1;
t156 = cos(qJ(1));
t242 = -t156 / 0.2e1;
t241 = t156 / 0.2e1;
t240 = t136 / 0.2e1;
t239 = t153 / 0.2e1;
t238 = t155 / 0.2e1;
t150 = cos(pkin(10));
t131 = pkin(5) * t150 + pkin(4);
t152 = -pkin(8) - qJ(5);
t149 = sin(pkin(10));
t211 = t149 * t156;
t196 = t136 * t211;
t212 = t138 * t156;
t145 = pkin(10) + qJ(6);
t135 = sin(t145);
t208 = t154 * t135;
t137 = cos(t145);
t214 = t137 * t156;
t85 = t136 * t214 - t208;
t207 = t154 * t137;
t215 = t136 * t156;
t86 = t135 * t215 + t207;
t35 = t86 * rSges(7,1) + t85 * rSges(7,2) + rSges(7,3) * t212;
t237 = pkin(5) * t196 + t154 * t131 - t152 * t212 + t35;
t236 = -qJ(5) - t152;
t235 = -t244 * t154 + t245 * t156;
t234 = t245 * t154 + t244 * t156;
t147 = t154 ^ 2;
t148 = t156 ^ 2;
t233 = 0.2e1 * t138;
t232 = m(5) / 0.2e1;
t231 = m(6) / 0.2e1;
t230 = m(7) / 0.2e1;
t228 = pkin(2) * t153;
t227 = pkin(5) * t149;
t132 = pkin(2) * t155 + pkin(1);
t128 = t156 * t132;
t143 = t156 * pkin(7);
t151 = -qJ(3) - pkin(7);
t209 = t151 * t156;
t226 = t154 * (t209 + t143 + (-pkin(1) + t132) * t154) + t156 * (-pkin(1) * t156 + t128 + (-pkin(7) - t151) * t154);
t225 = rSges(3,1) * t155;
t224 = rSges(3,2) * t153;
t223 = t156 * rSges(3,3);
t222 = Icges(3,4) * t153;
t221 = Icges(3,4) * t155;
t220 = Icges(4,4) * t136;
t219 = Icges(4,4) * t138;
t218 = Icges(5,6) * t136;
t217 = Icges(5,6) * t138;
t216 = qJ(4) * t136;
t213 = t138 * t154;
t210 = t150 * t156;
t206 = t154 * t149;
t205 = t154 * t150;
t204 = pkin(3) * t212 + qJ(4) * t215;
t203 = t154 * rSges(3,3) + t156 * t225;
t202 = t154 * pkin(4) + qJ(5) * t212;
t201 = t147 + t148;
t200 = t231 + t230;
t97 = t136 * t210 - t206;
t98 = t196 + t205;
t199 = t98 * rSges(6,1) + t97 * rSges(6,2) + rSges(6,3) * t212;
t87 = t135 * t156 + t136 * t207;
t88 = t136 * t208 - t214;
t30 = Icges(7,5) * t88 + Icges(7,6) * t87 + Icges(7,3) * t213;
t32 = Icges(7,4) * t88 + Icges(7,2) * t87 + Icges(7,6) * t213;
t34 = Icges(7,1) * t88 + Icges(7,4) * t87 + Icges(7,5) * t213;
t12 = t136 * t30 + (-t135 * t34 - t137 * t32) * t138;
t52 = Icges(7,3) * t136 + (-Icges(7,5) * t135 - Icges(7,6) * t137) * t138;
t53 = Icges(7,6) * t136 + (-Icges(7,4) * t135 - Icges(7,2) * t137) * t138;
t54 = Icges(7,5) * t136 + (-Icges(7,1) * t135 - Icges(7,4) * t137) * t138;
t15 = t213 * t52 + t53 * t87 + t54 * t88;
t198 = t12 / 0.2e1 + t15 / 0.2e1;
t29 = Icges(7,5) * t86 + Icges(7,6) * t85 + Icges(7,3) * t212;
t31 = Icges(7,4) * t86 + Icges(7,2) * t85 + Icges(7,6) * t212;
t33 = Icges(7,1) * t86 + Icges(7,4) * t85 + Icges(7,5) * t212;
t11 = t136 * t29 + (-t135 * t33 - t137 * t31) * t138;
t14 = t212 * t52 + t85 * t53 + t86 * t54;
t197 = t14 / 0.2e1 + t11 / 0.2e1;
t195 = -pkin(3) * t136 + qJ(4) * t138 - t228;
t194 = -rSges(4,1) * t136 - rSges(4,2) * t138 - t228;
t193 = t147 * (pkin(3) * t138 + t216) + t156 * t204 + t226;
t192 = -t154 * t151 + t128;
t191 = t232 + t200;
t190 = rSges(5,2) * t136 + rSges(5,3) * t138 + t195;
t189 = -t88 * rSges(7,1) - t87 * rSges(7,2);
t188 = -Icges(5,4) * t136 / 0.2e1 + Icges(4,5) * t240 + Icges(3,5) * t239 + Icges(3,6) * t238 + (-Icges(5,5) / 0.2e1 + Icges(4,6) / 0.2e1) * t138;
t100 = t136 * t206 - t210;
t99 = t136 * t205 + t211;
t187 = -t100 * rSges(6,1) - t99 * rSges(6,2);
t186 = -t224 + t225;
t185 = rSges(4,1) * t138 - rSges(4,2) * t136;
t184 = -t135 * t54 - t137 * t53;
t36 = rSges(7,3) * t213 - t189;
t55 = rSges(7,3) * t136 + (-rSges(7,1) * t135 - rSges(7,2) * t137) * t138;
t21 = -t136 * t36 + t213 * t55;
t22 = t136 * t35 - t212 * t55;
t177 = t154 * t22 + t156 * t21;
t160 = -qJ(5) * t136 + t195;
t158 = -t236 * t136 + t138 * t227 + t160 - t55;
t25 = t158 * t154;
t26 = t158 * t156;
t176 = t154 * t25 + t156 * t26;
t159 = t160 - rSges(6,3) * t136 - (-rSges(6,1) * t149 - rSges(6,2) * t150) * t138;
t27 = t159 * t154;
t28 = t159 * t156;
t175 = t154 * t27 + t156 * t28;
t174 = Icges(3,1) * t155 - t222;
t173 = Icges(4,1) * t138 - t220;
t172 = -Icges(3,2) * t153 + t221;
t171 = -Icges(4,2) * t136 + t219;
t167 = -Icges(5,2) * t138 + t218;
t166 = Icges(5,3) * t136 - t217;
t144 = t156 * pkin(4);
t165 = t154 * (qJ(5) * t213 - t144) + t156 * t202 + t193;
t164 = rSges(4,1) * t212 - rSges(4,2) * t215 + t154 * rSges(4,3);
t163 = t154 * rSges(5,1) - rSges(5,2) * t212 + rSges(5,3) * t215;
t161 = t192 + t204;
t18 = (t131 - t151) * t156 + (-t132 + (-qJ(4) - t227) * t136 + (-rSges(7,3) - pkin(3) + t152) * t138) * t154 + t189;
t19 = t161 + t237;
t23 = -t209 + t144 + (-t216 - t132 + (-rSges(6,3) - pkin(3) - qJ(5)) * t138) * t154 + t187;
t24 = t161 + t199 + t202;
t157 = (t154 * t19 + t156 * t18) * t230 + (t154 * t24 + t156 * t23) * t231;
t134 = t138 ^ 2;
t133 = t136 ^ 2;
t119 = rSges(2,1) * t156 - t154 * rSges(2,2);
t118 = -t154 * rSges(2,1) - rSges(2,2) * t156;
t117 = rSges(3,1) * t153 + rSges(3,2) * t155;
t66 = t194 * t156;
t65 = t194 * t154;
t62 = Icges(6,5) * t136 + (-Icges(6,1) * t149 - Icges(6,4) * t150) * t138;
t61 = Icges(6,6) * t136 + (-Icges(6,4) * t149 - Icges(6,2) * t150) * t138;
t59 = t154 * pkin(7) + (pkin(1) - t224) * t156 + t203;
t58 = t223 + t143 + (-pkin(1) - t186) * t154;
t51 = t136 * t52;
t49 = t164 + t192;
t48 = (rSges(4,3) - t151) * t156 + (-t132 - t185) * t154;
t47 = t190 * t156;
t46 = t190 * t154;
t45 = t156 * (-t156 * t224 + t203) + (t154 * t186 - t223) * t154;
t44 = Icges(6,1) * t100 + Icges(6,4) * t99 + Icges(6,5) * t213;
t43 = Icges(6,1) * t98 + Icges(6,4) * t97 + Icges(6,5) * t212;
t42 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t213;
t41 = Icges(6,4) * t98 + Icges(6,2) * t97 + Icges(6,6) * t212;
t40 = Icges(6,5) * t100 + Icges(6,6) * t99 + Icges(6,3) * t213;
t39 = Icges(6,5) * t98 + Icges(6,6) * t97 + Icges(6,3) * t212;
t38 = t161 + t163;
t37 = (rSges(5,1) - t151) * t156 + (-t132 + (rSges(5,2) - pkin(3)) * t138 + (-rSges(5,3) - qJ(4)) * t136) * t154;
t20 = t156 * t164 + (-t156 * rSges(4,3) + t154 * t185) * t154 + t226;
t17 = (t138 * t184 + t51) * t136;
t16 = (-t154 * t35 + t156 * t36) * t138;
t13 = t156 * t163 + (-t156 * rSges(5,1) + (-rSges(5,2) * t138 + rSges(5,3) * t136) * t154) * t154 + t193;
t10 = t213 * t30 + t32 * t87 + t34 * t88;
t9 = t213 * t29 + t31 * t87 + t33 * t88;
t8 = t212 * t30 + t85 * t32 + t86 * t34;
t7 = t212 * t29 + t85 * t31 + t86 * t33;
t6 = t154 * (rSges(6,3) * t213 - t187) + t156 * t199 + t165;
t5 = (-t202 + t237) * t156 + (-t156 * t131 + t144 + t36 + (t136 * t227 + t236 * t138) * t154) * t154 + t165;
t4 = -t10 * t156 + t9 * t154;
t3 = t7 * t154 - t156 * t8;
t2 = t15 * t136 + (t10 * t154 + t156 * t9) * t138;
t1 = t14 * t136 + (t154 * t8 + t156 * t7) * t138;
t50 = [t155 * (Icges(3,2) * t155 + t222) + t153 * (Icges(3,1) * t153 + t221) + Icges(2,3) + t51 + m(7) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t23 ^ 2 + t24 ^ 2) + m(4) * (t48 ^ 2 + t49 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t58 ^ 2 + t59 ^ 2) + m(2) * (t118 ^ 2 + t119 ^ 2) + (-t149 * t62 - t150 * t61 + t184 + t218 + t220 + (Icges(4,2) + Icges(5,3)) * t138) * t138 + (t219 + t217 + (-Icges(6,5) * t149 - Icges(6,6) * t150) * t138 + (Icges(4,1) + Icges(5,2) + Icges(6,3)) * t136) * t136; (-t153 * (-Icges(3,5) * t156 + t154 * t174) / 0.2e1 - t155 * (-Icges(3,6) * t156 + t154 * t172) / 0.2e1 - t100 * t62 / 0.2e1 - t61 * t99 / 0.2e1 + t188 * t156 - t198) * t156 + ((Icges(3,5) * t154 + t156 * t174) * t239 + (Icges(3,6) * t154 + t156 * t172) * t238 + t97 * t61 / 0.2e1 + t98 * t62 / 0.2e1 + t188 * t154 + t197) * t154 + m(4) * (t48 * t66 + t49 * t65) + m(5) * (t37 * t47 + t38 * t46) + m(6) * (t23 * t28 + t24 * t27) + m(7) * (t18 * t26 + t19 * t25) + m(3) * (-t154 * t59 - t156 * t58) * t117 + ((-t40 / 0.2e1 + Icges(4,5) * t241 + t173 * t243 + Icges(5,4) * t242 + t167 * t229) * t156 + (t39 / 0.2e1 + Icges(4,5) * t229 + t173 * t241 + Icges(5,4) * t243 + t167 * t242) * t154) * t136 + ((t149 * t44 / 0.2e1 + t150 * t42 / 0.2e1 + Icges(4,6) * t241 + t171 * t243 + Icges(5,5) * t242 + t166 * t229) * t156 + (-t149 * t43 / 0.2e1 - t150 * t41 / 0.2e1 + Icges(4,6) * t229 + t171 * t241 + Icges(5,5) * t243 + t166 * t242) * t154) * t138; m(7) * (t25 ^ 2 + t26 ^ 2 + t5 ^ 2) + m(6) * (t27 ^ 2 + t28 ^ 2 + t6 ^ 2) + m(5) * (t13 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(4) * (t20 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(3) * (t117 ^ 2 * t201 + t45 ^ 2) + (-t4 + (t100 * t44 + t40 * t213 + t42 * t99) * t156 + t235 * t148) * t156 + (t3 + (t39 * t212 + t97 * t41 + t98 * t43) * t154 + t234 * t147 + (-t100 * t43 + t235 * t154 + t234 * t156 - t212 * t40 - t213 * t39 - t41 * t99 - t97 * t42 - t98 * t44) * t156) * t154; m(7) * (t154 * t18 - t156 * t19) + m(6) * (t154 * t23 - t156 * t24) + m(4) * (t154 * t48 - t156 * t49) + m(5) * (t154 * t37 - t156 * t38); m(7) * (t154 * t26 - t156 * t25) + m(6) * (t154 * t28 - t156 * t27) + m(5) * (t154 * t47 - t156 * t46) + m(4) * (t154 * t66 - t156 * t65); 0.2e1 * (m(4) / 0.2e1 + t191) * t201; 0.2e1 * ((t154 * t38 + t156 * t37) * t232 + t157) * t136; m(7) * (t136 * t176 - t138 * t5) + m(6) * (t136 * t175 - t138 * t6) + m(5) * (-t138 * t13 + (t154 * t46 + t156 * t47) * t136); 0; 0.2e1 * t191 * (t133 * t201 + t134); t157 * t233; m(7) * (t136 * t5 + t138 * t176) + m(6) * (t136 * t6 + t138 * t175); 0; t200 * (-0.1e1 + t201) * t136 * t233; 0.2e1 * t200 * (t134 * t201 + t133); m(7) * (t18 * t21 + t19 * t22) + t17 + (t154 * t198 + t156 * t197) * t138; t2 * t242 + (t11 * t154 - t12 * t156) * t240 + m(7) * (t16 * t5 + t21 * t26 + t22 * t25) + t1 * t229 + (t4 * t229 + t3 * t241) * t138; m(7) * (t21 * t154 - t156 * t22); m(7) * (t136 * t177 - t16 * t138); m(7) * (t16 * t136 + t138 * t177); t136 * t17 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + (t156 * t1 + t154 * t2 + t136 * (t11 * t156 + t12 * t154)) * t138;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t50(1) t50(2) t50(4) t50(7) t50(11) t50(16); t50(2) t50(3) t50(5) t50(8) t50(12) t50(17); t50(4) t50(5) t50(6) t50(9) t50(13) t50(18); t50(7) t50(8) t50(9) t50(10) t50(14) t50(19); t50(11) t50(12) t50(13) t50(14) t50(15) t50(20); t50(16) t50(17) t50(18) t50(19) t50(20) t50(21);];
Mq  = res;
