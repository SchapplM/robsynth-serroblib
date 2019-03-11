% Calculate joint inertia matrix for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:21
% EndTime: 2019-03-09 03:01:26
% DurationCPUTime: 2.08s
% Computational Cost: add. (6112->344), mult. (5097->485), div. (0->0), fcn. (5334->10), ass. (0->162)
t143 = -qJ(6) - pkin(8);
t227 = rSges(7,3) - t143;
t226 = Icges(4,3) + Icges(5,3);
t141 = qJ(3) + pkin(10);
t136 = sin(t141);
t138 = cos(t141);
t146 = sin(qJ(3));
t149 = cos(qJ(3));
t225 = Icges(4,5) * t149 + Icges(5,5) * t138 - Icges(4,6) * t146 - Icges(5,6) * t136;
t145 = sin(qJ(5));
t148 = cos(qJ(5));
t82 = -Icges(7,3) * t138 + (Icges(7,5) * t148 - Icges(7,6) * t145) * t136;
t83 = -Icges(6,3) * t138 + (Icges(6,5) * t148 - Icges(6,6) * t145) * t136;
t224 = -t82 - t83;
t86 = -Icges(7,6) * t138 + (Icges(7,4) * t148 - Icges(7,2) * t145) * t136;
t87 = -Icges(6,6) * t138 + (Icges(6,4) * t148 - Icges(6,2) * t145) * t136;
t223 = (-t86 - t87) * t145;
t142 = qJ(1) + pkin(9);
t139 = cos(t142);
t179 = t139 * t148;
t137 = sin(t142);
t184 = t137 * t145;
t100 = -t138 * t184 - t179;
t180 = t139 * t145;
t183 = t137 * t148;
t101 = t138 * t183 - t180;
t188 = t136 * t137;
t47 = Icges(7,5) * t101 + Icges(7,6) * t100 + Icges(7,3) * t188;
t51 = Icges(7,4) * t101 + Icges(7,2) * t100 + Icges(7,6) * t188;
t55 = Icges(7,1) * t101 + Icges(7,4) * t100 + Icges(7,5) * t188;
t12 = t100 * t51 + t101 * t55 + t47 * t188;
t102 = -t138 * t180 + t183;
t103 = t138 * t179 + t184;
t187 = t136 * t139;
t48 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t187;
t52 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t187;
t56 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t187;
t13 = t100 * t52 + t101 * t56 + t48 * t188;
t49 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t188;
t53 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t188;
t57 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t188;
t14 = t100 * t53 + t101 * t57 + t49 * t188;
t50 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t187;
t54 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t187;
t58 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t187;
t15 = t100 * t54 + t101 * t58 + t50 * t188;
t90 = -Icges(7,5) * t138 + (Icges(7,1) * t148 - Icges(7,4) * t145) * t136;
t28 = t100 * t86 + t101 * t90 + t82 * t188;
t91 = -Icges(6,5) * t138 + (Icges(6,1) * t148 - Icges(6,4) * t145) * t136;
t29 = t100 * t87 + t101 * t91 + t83 * t188;
t222 = (-t28 - t29) * t138 + ((t13 + t15) * t139 + (t12 + t14) * t137) * t136;
t16 = t102 * t51 + t103 * t55 + t47 * t187;
t17 = t102 * t52 + t103 * t56 + t48 * t187;
t18 = t102 * t53 + t103 * t57 + t49 * t187;
t19 = t102 * t54 + t103 * t58 + t50 * t187;
t30 = t102 * t86 + t103 * t90 + t82 * t187;
t31 = t102 * t87 + t103 * t91 + t83 * t187;
t221 = (-t30 - t31) * t138 + ((t17 + t19) * t139 + (t16 + t18) * t137) * t136;
t220 = t136 / 0.2e1;
t219 = t138 / 0.2e1;
t218 = t146 / 0.2e1;
t217 = t149 / 0.2e1;
t22 = -t138 * t47 + (-t145 * t51 + t148 * t55) * t136;
t24 = -t138 * t49 + (-t145 * t53 + t148 * t57) * t136;
t216 = -t22 - t24;
t23 = -t138 * t48 + (-t145 * t52 + t148 * t56) * t136;
t25 = -t138 * t50 + (-t145 * t54 + t148 * t58) * t136;
t215 = t23 + t25;
t214 = (t90 + t91) * t136 * t148;
t213 = -t225 * t137 + t226 * t139;
t212 = t226 * t137 + t225 * t139;
t131 = pkin(5) * t148 + pkin(4);
t182 = t138 * t139;
t211 = t103 * rSges(7,1) + t102 * rSges(7,2) + pkin(5) * t184 + t131 * t182 + t227 * t187;
t134 = t137 ^ 2;
t210 = t138 ^ 2;
t135 = t139 ^ 2;
t208 = -t138 / 0.2e1;
t123 = rSges(4,1) * t146 + rSges(4,2) * t149;
t206 = m(4) * t123;
t147 = sin(qJ(1));
t205 = pkin(1) * t147;
t204 = pkin(3) * t146;
t203 = pkin(4) * t138;
t202 = -pkin(4) + t131;
t201 = pkin(8) + t143;
t200 = t136 * t223 + t224 * t138 + t214;
t166 = -rSges(7,1) * t101 - rSges(7,2) * t100;
t199 = -pkin(5) * t180 + (-t201 * t136 + t202 * t138) * t137 + rSges(7,3) * t188 - t166;
t178 = pkin(4) * t182 + pkin(8) * t187;
t198 = -t178 + t211;
t132 = pkin(3) * t149 + pkin(2);
t118 = t139 * t132;
t130 = t139 * pkin(7);
t144 = -qJ(4) - pkin(7);
t181 = t139 * t144;
t197 = t137 * (t181 + t130 + (-pkin(2) + t132) * t137) + t139 * (-pkin(2) * t139 + t118 + (-pkin(7) - t144) * t137);
t196 = (t201 - rSges(7,3)) * t138 + (rSges(7,1) * t148 - rSges(7,2) * t145 + t202) * t136;
t195 = rSges(4,1) * t149;
t194 = rSges(4,2) * t146;
t193 = t139 * rSges(4,3);
t192 = Icges(4,4) * t146;
t191 = Icges(4,4) * t149;
t190 = Icges(5,4) * t136;
t189 = Icges(5,4) * t138;
t177 = t137 * rSges(4,3) + t139 * t195;
t176 = t134 + t135;
t62 = t103 * rSges(6,1) + t102 * rSges(6,2) + rSges(6,3) * t187;
t175 = Icges(4,5) * t218 + Icges(5,5) * t220 + Icges(4,6) * t217 + Icges(5,6) * t219;
t174 = -rSges(5,1) * t136 - rSges(5,2) * t138 - t204;
t173 = -pkin(4) * t136 + pkin(8) * t138 - t204;
t172 = t134 * (pkin(8) * t136 + t203) + t139 * t178 + t197;
t95 = -rSges(6,3) * t138 + (rSges(6,1) * t148 - rSges(6,2) * t145) * t136;
t171 = t173 - t95;
t150 = cos(qJ(1));
t140 = t150 * pkin(1);
t170 = -t137 * t144 + t118 + t140;
t169 = -t194 + t195;
t168 = rSges(5,1) * t138 - rSges(5,2) * t136;
t167 = -rSges(6,1) * t101 - rSges(6,2) * t100;
t161 = t173 - t196;
t160 = Icges(4,1) * t149 - t192;
t159 = Icges(5,1) * t138 - t190;
t158 = -Icges(4,2) * t146 + t191;
t157 = -Icges(5,2) * t136 + t189;
t154 = rSges(5,1) * t182 - rSges(5,2) * t187 + t137 * rSges(5,3);
t152 = t22 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1 + t24 / 0.2e1;
t151 = t25 / 0.2e1 + t23 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1;
t125 = rSges(2,1) * t150 - t147 * rSges(2,2);
t124 = -t147 * rSges(2,1) - rSges(2,2) * t150;
t105 = rSges(3,1) * t139 - rSges(3,2) * t137 + t140;
t104 = -rSges(3,1) * t137 - rSges(3,2) * t139 - t205;
t79 = t174 * t139;
t78 = t174 * t137;
t66 = pkin(7) * t137 + t140 + (pkin(2) - t194) * t139 + t177;
t65 = t193 - t205 + t130 + (-pkin(2) - t169) * t137;
t64 = t154 + t170;
t63 = -t205 + (rSges(5,3) - t144) * t139 + (-t132 - t168) * t137;
t60 = rSges(6,3) * t188 - t167;
t46 = t171 * t139;
t45 = t171 * t137;
t42 = t139 * (-t139 * t194 + t177) + (t169 * t137 - t193) * t137;
t41 = t161 * t139;
t40 = t161 * t137;
t37 = -t138 * t62 - t95 * t187;
t36 = t138 * t60 + t95 * t188;
t35 = t170 + t62 + t178;
t34 = -t205 - t181 + (-t203 - t132 + (-rSges(6,3) - pkin(8)) * t136) * t137 + t167;
t33 = t170 + t211;
t32 = -t205 + (pkin(5) * t145 - t144) * t139 + (-t131 * t138 - t227 * t136 - t132) * t137 + t166;
t27 = (-t137 * t62 + t139 * t60) * t136;
t26 = t139 * t154 + (-rSges(5,3) * t139 + t168 * t137) * t137 + t197;
t21 = -t198 * t138 - t196 * t187;
t20 = t199 * t138 + t196 * t188;
t11 = t137 * t60 + t139 * t62 + t172;
t10 = (-t198 * t137 + t199 * t139) * t136;
t9 = t199 * t137 + t198 * t139 + t172;
t8 = t137 * t19 - t139 * t18;
t7 = t137 * t17 - t139 * t16;
t6 = t137 * t15 - t139 * t14;
t5 = -t12 * t139 + t13 * t137;
t1 = [t149 * (Icges(4,2) * t149 + t192) + t146 * (Icges(4,1) * t146 + t191) + Icges(2,3) + Icges(3,3) + (Icges(5,2) * t138 + t190 + t224) * t138 + (Icges(5,1) * t136 + t189 + t223) * t136 + m(7) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t63 ^ 2 + t64 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(3) * (t104 ^ 2 + t105 ^ 2) + m(2) * (t124 ^ 2 + t125 ^ 2) + t214; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t32 * t41 + t33 * t40) + m(6) * (t34 * t46 + t35 * t45) + m(5) * (t63 * t79 + t64 * t78) + (-t65 * t206 - t146 * (-Icges(4,5) * t139 + t160 * t137) / 0.2e1 - t149 * (-Icges(4,6) * t139 + t158 * t137) / 0.2e1 - t136 * (-Icges(5,5) * t139 + t159 * t137) / 0.2e1 + (-Icges(5,6) * t139 + t157 * t137) * t208 + t175 * t139 - t152) * t139 + (t175 * t137 - t66 * t206 + (Icges(4,6) * t137 + t158 * t139) * t217 + (Icges(4,5) * t137 + t160 * t139) * t218 + (Icges(5,6) * t137 + t157 * t139) * t219 + (Icges(5,5) * t137 + t159 * t139) * t220 + t151) * t137; m(4) * t42 + m(5) * t26 + m(6) * t11 + m(7) * t9; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t11 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t26 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t176 * t123 ^ 2 + t42 ^ 2) + (t213 * t135 - t5 - t6) * t139 + (t7 + t8 + t212 * t134 + (t213 * t137 + t212 * t139) * t139) * t137; m(7) * (t137 * t32 - t139 * t33) + m(6) * (t137 * t34 - t139 * t35) + m(5) * (t137 * t63 - t139 * t64); 0; m(7) * (t137 * t41 - t139 * t40) + m(6) * (t137 * t46 - t139 * t45) + m(5) * (t137 * t79 - t139 * t78); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t176; -t200 * t138 + m(7) * (t20 * t32 + t21 * t33) + m(6) * (t34 * t36 + t35 * t37) + (t152 * t137 + t151 * t139) * t136; m(6) * t27 + m(7) * t10; m(7) * (t10 * t9 + t20 * t41 + t21 * t40) + m(6) * (t11 * t27 + t36 * t46 + t37 * t45) + ((t8 / 0.2e1 + t7 / 0.2e1) * t139 + (t6 / 0.2e1 + t5 / 0.2e1) * t137) * t136 + t221 * t137 / 0.2e1 + (t215 * t137 + t216 * t139) * t208 - t222 * t139 / 0.2e1; m(6) * (t137 * t36 - t139 * t37) + m(7) * (t137 * t20 - t139 * t21); m(7) * (t10 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t27 ^ 2 + t36 ^ 2 + t37 ^ 2) + t200 * t210 + ((-t215 * t138 + t221) * t139 + (t216 * t138 + t222) * t137) * t136; m(7) * (t137 * t33 + t139 * t32) * t136; -m(7) * t138; m(7) * (-t138 * t9 + (t137 * t40 + t139 * t41) * t136); 0; m(7) * (-t10 * t138 + (t137 * t21 + t139 * t20) * t136); m(7) * (t176 * t136 ^ 2 + t210);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
