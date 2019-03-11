% Calculate joint inertia matrix for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:59
% EndTime: 2019-03-09 03:31:04
% DurationCPUTime: 2.66s
% Computational Cost: add. (2316->333), mult. (5557->489), div. (0->0), fcn. (5844->6), ass. (0->159)
t242 = -Icges(4,4) - Icges(5,6);
t241 = Icges(4,1) + Icges(5,2);
t240 = Icges(4,2) + Icges(5,3);
t151 = cos(qJ(3));
t239 = t242 * t151;
t148 = sin(qJ(3));
t238 = t242 * t148;
t237 = Icges(5,1) + Icges(4,3);
t236 = t240 * t151 - t238;
t235 = t241 * t148 - t239;
t234 = (-Icges(5,5) + Icges(4,6)) * t151 + (-Icges(5,4) + Icges(4,5)) * t148;
t147 = sin(qJ(5));
t149 = sin(qJ(1));
t150 = cos(qJ(5));
t152 = cos(qJ(1));
t192 = t151 * t152;
t103 = t147 * t149 - t150 * t192;
t104 = t147 * t192 + t149 * t150;
t227 = rSges(7,3) + qJ(6);
t228 = rSges(7,1) + pkin(5);
t233 = -t227 * t103 - t228 * t104;
t232 = -t149 / 0.2e1;
t231 = t149 / 0.2e1;
t230 = -t152 / 0.2e1;
t229 = t152 / 0.2e1;
t195 = t148 * t150;
t75 = Icges(6,3) * t151 + (Icges(6,5) * t147 + Icges(6,6) * t150) * t148;
t78 = Icges(7,2) * t151 + (Icges(7,4) * t147 - Icges(7,6) * t150) * t148;
t79 = Icges(6,6) * t151 + (Icges(6,4) * t147 + Icges(6,2) * t150) * t148;
t82 = Icges(7,4) * t151 + (Icges(7,1) * t147 - Icges(7,5) * t150) * t148;
t83 = Icges(6,5) * t151 + (Icges(6,1) * t147 + Icges(6,4) * t150) * t148;
t226 = t79 * t195 + (t82 + t83) * t147 * t148 + (t75 + t78) * t151;
t194 = t148 * t152;
t47 = Icges(7,5) * t104 - Icges(7,6) * t194 + Icges(7,3) * t103;
t51 = Icges(7,4) * t104 - Icges(7,2) * t194 + Icges(7,6) * t103;
t55 = Icges(7,1) * t104 - Icges(7,4) * t194 + Icges(7,5) * t103;
t11 = t103 * t47 + t104 * t55 - t194 * t51;
t193 = t149 * t151;
t105 = t147 * t152 + t150 * t193;
t106 = -t147 * t193 + t150 * t152;
t196 = t148 * t149;
t48 = Icges(7,5) * t106 + Icges(7,6) * t196 + Icges(7,3) * t105;
t52 = Icges(7,4) * t106 + Icges(7,2) * t196 + Icges(7,6) * t105;
t56 = Icges(7,1) * t106 + Icges(7,4) * t196 + Icges(7,5) * t105;
t12 = t103 * t48 + t104 * t56 - t194 * t52;
t49 = Icges(6,5) * t104 - Icges(6,6) * t103 - Icges(6,3) * t194;
t53 = Icges(6,4) * t104 - Icges(6,2) * t103 - Icges(6,6) * t194;
t57 = Icges(6,1) * t104 - Icges(6,4) * t103 - Icges(6,5) * t194;
t13 = -t103 * t53 + t104 * t57 - t194 * t49;
t50 = Icges(6,5) * t106 - Icges(6,6) * t105 + Icges(6,3) * t196;
t54 = Icges(6,4) * t106 - Icges(6,2) * t105 + Icges(6,6) * t196;
t58 = Icges(6,1) * t106 - Icges(6,4) * t105 + Icges(6,5) * t196;
t14 = -t103 * t54 + t104 * t58 - t194 * t50;
t74 = Icges(7,6) * t151 + (Icges(7,5) * t147 - Icges(7,3) * t150) * t148;
t28 = t103 * t74 + t104 * t82 - t194 * t78;
t29 = -t103 * t79 + t104 * t83 - t194 * t75;
t225 = (t29 + t28) * t151 + ((-t11 - t13) * t152 + (t12 + t14) * t149) * t148;
t15 = t105 * t47 + t106 * t55 + t196 * t51;
t16 = t105 * t48 + t106 * t56 + t196 * t52;
t17 = -t105 * t53 + t106 * t57 + t196 * t49;
t18 = -t105 * t54 + t106 * t58 + t196 * t50;
t30 = t105 * t74 + t106 * t82 + t196 * t78;
t31 = -t105 * t79 + t106 * t83 + t196 * t75;
t224 = (t30 + t31) * t151 + ((-t15 - t17) * t152 + (t16 + t18) * t149) * t148;
t223 = t151 / 0.2e1;
t205 = rSges(5,2) * t148;
t155 = -t205 + (-rSges(5,3) - qJ(4)) * t151;
t135 = pkin(3) * t194;
t137 = t152 * qJ(2);
t189 = t135 + t137;
t214 = -pkin(1) - pkin(7);
t42 = t155 * t152 + (-rSges(5,1) + t214) * t149 + t189;
t132 = pkin(3) * t196;
t188 = t152 * pkin(1) + t149 * qJ(2);
t184 = t152 * pkin(7) + t188;
t43 = t152 * rSges(5,1) + t149 * t155 + t132 + t184;
t222 = m(5) * (t149 * t42 - t152 * t43);
t134 = pkin(8) * t194;
t154 = t134 + (-pkin(4) + t214) * t149 + t189;
t178 = -t104 * rSges(6,1) + t103 * rSges(6,2);
t198 = qJ(4) * t151;
t34 = (rSges(6,3) * t148 - t198) * t152 + t154 + t178;
t109 = -qJ(4) * t193 + t132;
t191 = -t152 * pkin(4) - pkin(8) * t196 - t109;
t153 = t184 - t191;
t62 = t106 * rSges(6,1) - t105 * rSges(6,2) + rSges(6,3) * t196;
t35 = t153 + t62;
t221 = m(6) * (t149 * t34 - t152 * t35);
t26 = (rSges(7,2) * t148 - t198) * t152 + t154 + t233;
t207 = rSges(7,2) * t196 + t105 * t227 + t106 * t228;
t27 = t153 + t207;
t220 = m(7) * (t149 * t26 - t152 * t27);
t20 = t151 * t51 + (t147 * t55 - t150 * t47) * t148;
t22 = t151 * t49 + (t147 * t57 + t150 * t53) * t148;
t219 = t20 + t22;
t21 = t151 * t52 + (t147 * t56 - t150 * t48) * t148;
t23 = t151 * t50 + (t147 * t58 + t150 * t54) * t148;
t218 = t21 + t23;
t217 = t234 * t149 + t237 * t152;
t216 = t237 * t149 - t234 * t152;
t215 = (rSges(4,1) * t148 + rSges(4,2) * t151) * t152;
t144 = t149 ^ 2;
t146 = t152 ^ 2;
t123 = rSges(4,1) * t151 - rSges(4,2) * t148;
t210 = m(4) * t123;
t209 = (-t195 * t74 + t226) * t151;
t208 = -rSges(7,2) * t194 - t233;
t96 = t152 * (qJ(4) * t192 - t135);
t204 = t152 * (t149 * pkin(4) - t134) + t96;
t203 = rSges(7,2) * t151 + (t147 * t228 - t150 * t227) * t148;
t121 = pkin(3) * t151 + qJ(4) * t148;
t110 = t149 * t121;
t190 = pkin(8) * t193 + t110;
t187 = t144 + t146;
t185 = rSges(4,1) * t196 + rSges(4,2) * t193 + t152 * rSges(4,3);
t183 = Icges(4,5) * t223 - Icges(5,4) * t151 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(5,5) / 0.2e1) * t148;
t182 = -pkin(8) * t151 - t121;
t181 = t148 * t203;
t180 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t177 = rSges(5,3) * t151 + t205;
t24 = -t149 * t181 + t151 * t207;
t25 = -t151 * t208 - t152 * t181;
t172 = t25 * t149 - t152 * t24;
t93 = rSges(6,3) * t151 + (rSges(6,1) * t147 + rSges(6,2) * t150) * t148;
t38 = t151 * t62 - t196 * t93;
t60 = -rSges(6,3) * t194 - t178;
t39 = -t151 * t60 - t194 * t93;
t169 = t39 * t149 - t152 * t38;
t40 = t149 * t203 + t190;
t41 = (t182 - t203) * t152;
t168 = t40 * t149 - t152 * t41;
t45 = t149 * t93 + t190;
t46 = (t182 - t93) * t152;
t166 = t45 * t149 - t152 * t46;
t122 = -rSges(5,2) * t151 + rSges(5,3) * t148;
t67 = t122 * t149 + t110;
t68 = (-t121 - t122) * t152;
t165 = t67 * t149 - t152 * t68;
t158 = t103 * t152 - t105 * t149;
t157 = -t22 / 0.2e1 - t20 / 0.2e1 - t29 / 0.2e1 - t28 / 0.2e1;
t156 = t23 / 0.2e1 + t21 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1;
t143 = t148 ^ 2;
t124 = rSges(2,1) * t152 - t149 * rSges(2,2);
t120 = -t149 * rSges(2,1) - rSges(2,2) * t152;
t95 = -rSges(3,2) * t152 + t149 * rSges(3,3) + t188;
t94 = rSges(3,3) * t152 + t137 + (rSges(3,2) - pkin(1)) * t149;
t66 = t184 + t185;
t65 = t137 + t215 + (-rSges(4,3) + t214) * t149;
t44 = -t149 * t185 + (t149 * rSges(4,3) - t215) * t152;
t33 = t96 + t177 * t146 + (t149 * t177 - t109) * t149;
t32 = (t149 * t60 + t152 * t62) * t148;
t19 = t152 * t60 + (-t62 + t191) * t149 + t204;
t10 = (t149 * t208 + t152 * t207) * t148;
t9 = t208 * t152 + (t191 - t207) * t149 + t204;
t8 = t17 * t149 + t152 * t18;
t7 = t15 * t149 + t152 * t16;
t6 = t13 * t149 + t14 * t152;
t5 = t11 * t149 + t12 * t152;
t1 = [Icges(3,1) + Icges(2,3) + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(3) * (t94 ^ 2 + t95 ^ 2) + m(2) * (t120 ^ 2 + t124 ^ 2) + (t241 * t151 + t238) * t151 + (t240 * t148 - t150 * t74 + t239) * t148 + t226; t220 + t221 + t222 + m(4) * (t149 * t65 - t152 * t66) + m(3) * (t149 * t94 - t152 * t95); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t180) * t187; m(7) * (t26 * t40 + t27 * t41) + m(6) * (t34 * t45 + t35 * t46) + m(5) * (t42 * t67 + t43 * t68) + (-t66 * t210 + t183 * t152 + (Icges(5,4) * t230 + Icges(4,5) * t229 + t235 * t231) * t151 + (Icges(5,5) * t229 + Icges(4,6) * t230 + t236 * t232) * t148 + t156) * t152 + (t65 * t210 + t183 * t149 + (Icges(5,4) * t232 + Icges(4,5) * t231 + t235 * t230) * t151 + (Icges(5,5) * t231 + Icges(4,6) * t232 + t236 * t229) * t148 - t157) * t149; m(5) * t165 + m(6) * t166 + m(7) * t168 + t187 * t210; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t19 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t33 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t123 ^ 2 * t187 + t44 ^ 2) + (t146 * t217 + t7 + t8) * t152 + (t5 + t6 + t216 * t144 + (t149 * t217 + t152 * t216) * t152) * t149; 0.2e1 * (-t220 / 0.2e1 - t221 / 0.2e1 - t222 / 0.2e1) * t151; -0.2e1 * t180 * t187 * t151; m(7) * (t148 * t9 - t151 * t168) + m(6) * (t148 * t19 - t151 * t166) + m(5) * (t148 * t33 - t151 * t165); 0.2e1 * t180 * (t151 ^ 2 * t187 + t143); m(7) * (t24 * t27 + t25 * t26) + m(6) * (t34 * t39 + t35 * t38) + (t149 * t156 + t152 * t157) * t148 + t209; m(6) * t169 + m(7) * t172; m(7) * (t10 * t9 + t24 * t41 + t25 * t40) + m(6) * (t19 * t32 + t38 * t46 + t39 * t45) + ((-t6 / 0.2e1 - t5 / 0.2e1) * t152 + (t8 / 0.2e1 + t7 / 0.2e1) * t149) * t148 + t225 * t231 + (t149 * t219 + t152 * t218) * t223 + t224 * t229; m(6) * (t32 * t148 - t151 * t169) + m(7) * (t10 * t148 - t151 * t172); t209 * t151 + m(7) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t38 ^ 2 + t39 ^ 2) + ((-t151 * t219 - t225) * t152 + (t151 * t218 + t224) * t149) * t148; m(7) * (t103 * t27 + t105 * t26); -m(7) * t158; m(7) * (t103 * t41 + t105 * t40 - t195 * t9); m(7) * (-t143 * t150 + t151 * t158); m(7) * (-t10 * t195 + t103 * t24 + t105 * t25); m(7) * (t143 * t150 ^ 2 + t103 ^ 2 + t105 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
