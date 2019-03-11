% Calculate joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:19
% EndTime: 2019-03-09 02:20:23
% DurationCPUTime: 1.94s
% Computational Cost: add. (8693->332), mult. (6277->478), div. (0->0), fcn. (6740->12), ass. (0->171)
t146 = pkin(11) + qJ(4);
t139 = sin(t146);
t212 = Icges(5,5) * t139;
t211 = t212 / 0.2e1;
t154 = cos(qJ(5));
t136 = t154 * pkin(5) + pkin(4);
t156 = -pkin(9) - pkin(8);
t141 = cos(t146);
t147 = qJ(1) + pkin(10);
t142 = cos(t147);
t185 = t141 * t142;
t140 = sin(t147);
t152 = sin(qJ(5));
t187 = t140 * t152;
t190 = t139 * t142;
t148 = qJ(5) + qJ(6);
t143 = sin(t148);
t184 = t142 * t143;
t144 = cos(t148);
t188 = t140 * t144;
t107 = -t141 * t184 + t188;
t183 = t142 * t144;
t189 = t140 * t143;
t108 = t141 * t183 + t189;
t66 = t108 * rSges(7,1) + t107 * rSges(7,2) + rSges(7,3) * t190;
t210 = pkin(5) * t187 + t136 * t185 - t156 * t190 + t66;
t137 = t140 ^ 2;
t138 = t142 ^ 2;
t191 = t139 * t140;
t105 = -t141 * t189 - t183;
t106 = t141 * t188 - t184;
t59 = Icges(7,5) * t106 + Icges(7,6) * t105 + Icges(7,3) * t191;
t61 = Icges(7,4) * t106 + Icges(7,2) * t105 + Icges(7,6) * t191;
t63 = Icges(7,1) * t106 + Icges(7,4) * t105 + Icges(7,5) * t191;
t19 = t105 * t61 + t106 * t63 + t59 * t191;
t60 = Icges(7,5) * t108 + Icges(7,6) * t107 + Icges(7,3) * t190;
t62 = Icges(7,4) * t108 + Icges(7,2) * t107 + Icges(7,6) * t190;
t64 = Icges(7,1) * t108 + Icges(7,4) * t107 + Icges(7,5) * t190;
t20 = t105 * t62 + t106 * t64 + t60 * t191;
t93 = -Icges(7,3) * t141 + (Icges(7,5) * t144 - Icges(7,6) * t143) * t139;
t94 = -Icges(7,6) * t141 + (Icges(7,4) * t144 - Icges(7,2) * t143) * t139;
t95 = -Icges(7,5) * t141 + (Icges(7,1) * t144 - Icges(7,4) * t143) * t139;
t39 = t105 * t94 + t106 * t95 + t93 * t191;
t5 = -t39 * t141 + (t140 * t19 + t142 * t20) * t139;
t21 = t107 * t61 + t108 * t63 + t59 * t190;
t22 = t107 * t62 + t108 * t64 + t60 * t190;
t40 = t107 * t94 + t108 * t95 + t93 * t190;
t6 = -t40 * t141 + (t140 * t21 + t142 * t22) * t139;
t209 = t6 * t190 + t5 * t191;
t208 = t140 / 0.2e1;
t207 = -t141 / 0.2e1;
t206 = -t142 / 0.2e1;
t205 = t142 / 0.2e1;
t122 = t139 * rSges(5,1) + t141 * rSges(5,2);
t204 = m(5) * t122;
t203 = pkin(4) * t141;
t153 = sin(qJ(1));
t202 = t153 * pkin(1);
t201 = -pkin(4) + t136;
t200 = pkin(8) + t156;
t178 = pkin(4) * t185 + pkin(8) * t190;
t199 = -t178 + t210;
t165 = -t106 * rSges(7,1) - t105 * rSges(7,2);
t65 = rSges(7,3) * t191 - t165;
t96 = -t141 * rSges(7,3) + (rSges(7,1) * t144 - rSges(7,2) * t143) * t139;
t47 = t141 * t65 + t96 * t191;
t86 = t201 * t139 + t200 * t141;
t198 = -t86 - t96;
t197 = t137 * (pkin(8) * t139 + t203) + t142 * t178;
t196 = t143 * t94;
t83 = t139 * t144 * t95;
t46 = -t139 * t196 - t141 * t93 + t83;
t195 = t46 * t141;
t194 = rSges(4,3) + qJ(3);
t192 = Icges(5,4) * t141;
t186 = t140 * t154;
t182 = t142 * t152;
t181 = t142 * t154;
t102 = -Icges(6,6) * t141 + (Icges(6,4) * t154 - Icges(6,2) * t152) * t139;
t180 = t152 * t102;
t104 = -t141 * rSges(6,3) + (rSges(6,1) * t154 - rSges(6,2) * t152) * t139;
t123 = t139 * pkin(4) - t141 * pkin(8);
t179 = -t104 - t123;
t177 = t137 + t138;
t176 = -t123 + t198;
t111 = -t141 * t187 - t181;
t112 = t141 * t186 - t182;
t69 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t191;
t71 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t191;
t73 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t191;
t33 = -t141 * t69 + (-t152 * t71 + t154 * t73) * t139;
t101 = -Icges(6,3) * t141 + (Icges(6,5) * t154 - Icges(6,6) * t152) * t139;
t103 = -Icges(6,5) * t141 + (Icges(6,1) * t154 - Icges(6,4) * t152) * t139;
t42 = t101 * t191 + t111 * t102 + t112 * t103;
t175 = t33 / 0.2e1 + t42 / 0.2e1;
t113 = -t141 * t182 + t186;
t114 = t141 * t181 + t187;
t70 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t190;
t72 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t190;
t74 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t190;
t34 = -t141 * t70 + (-t152 * t72 + t154 * t74) * t139;
t43 = t101 * t190 + t113 * t102 + t114 * t103;
t174 = t43 / 0.2e1 + t34 / 0.2e1;
t76 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t190;
t173 = t191 / 0.2e1;
t172 = t190 / 0.2e1;
t29 = -t141 * t59 + (-t143 * t61 + t144 * t63) * t139;
t30 = -t141 * t60 + (-t143 * t62 + t144 * t64) * t139;
t171 = (t29 + t39) * t173 + (t30 + t40) * t172;
t9 = -t195 + (t140 * t29 + t142 * t30) * t139;
t170 = -t141 * t9 + t209;
t12 = t20 * t140 - t19 * t142;
t13 = t22 * t140 - t21 * t142;
t169 = t12 * t173 + t13 * t172 + t5 * t206 + t6 * t208 + (t30 * t140 - t29 * t142) * t207;
t150 = cos(pkin(11));
t135 = t150 * pkin(3) + pkin(2);
t155 = cos(qJ(1));
t145 = t155 * pkin(1);
t151 = -pkin(7) - qJ(3);
t168 = t142 * t135 - t140 * t151 + t145;
t167 = rSges(5,1) * t141 - rSges(5,2) * t139;
t166 = -t112 * rSges(6,1) - t111 * rSges(6,2);
t161 = -Icges(5,2) * t139 + t192;
t160 = Icges(5,5) * t141 - Icges(5,6) * t139;
t159 = rSges(5,1) * t185 - rSges(5,2) * t190 + t140 * rSges(5,3);
t149 = sin(pkin(11));
t158 = rSges(4,1) * t150 - rSges(4,2) * t149 + pkin(2);
t132 = t155 * rSges(2,1) - t153 * rSges(2,2);
t131 = -t153 * rSges(2,1) - t155 * rSges(2,2);
t119 = Icges(5,6) * t141 + t212;
t116 = t142 * rSges(3,1) - t140 * rSges(3,2) + t145;
t115 = -t140 * rSges(3,1) - t142 * rSges(3,2) - t202;
t88 = Icges(5,3) * t140 + t160 * t142;
t87 = -Icges(5,3) * t142 + t160 * t140;
t85 = t139 * t154 * t103;
t82 = t194 * t140 + t158 * t142 + t145;
t81 = -t158 * t140 + t194 * t142 - t202;
t80 = t179 * t142;
t79 = t179 * t140;
t78 = t159 + t168;
t77 = -t202 + (rSges(5,3) - t151) * t142 + (-t135 - t167) * t140;
t75 = rSges(6,3) * t191 - t166;
t67 = -pkin(5) * t182 + (-t200 * t139 + t201 * t141) * t140;
t58 = t142 * t159 + (-t142 * rSges(5,3) + t167 * t140) * t140;
t56 = t65 * t190;
t55 = t176 * t142;
t54 = t176 * t140;
t53 = -t141 * t101 - t139 * t180 + t85;
t52 = -t104 * t190 - t141 * t76;
t51 = t104 * t191 + t141 * t75;
t50 = t168 + t76 + t178;
t49 = -t202 - t142 * t151 + (-t203 - t135 + (-rSges(6,3) - pkin(8)) * t139) * t140 + t166;
t48 = -t141 * t66 - t96 * t190;
t45 = t168 + t210;
t44 = -t202 + (pkin(5) * t152 - t151) * t142 + (-t136 * t141 - t135 + (-rSges(7,3) + t156) * t139) * t140 + t165;
t41 = (-t140 * t76 + t142 * t75) * t139;
t38 = -t66 * t191 + t56;
t35 = t140 * t75 + t142 * t76 + t197;
t32 = -t199 * t141 + t198 * t190;
t31 = t141 * t67 + t86 * t191 + t47;
t28 = t113 * t72 + t114 * t74 + t70 * t190;
t27 = t113 * t71 + t114 * t73 + t69 * t190;
t26 = t111 * t72 + t112 * t74 + t70 * t191;
t25 = t111 * t71 + t112 * t73 + t69 * t191;
t18 = t56 + (-t199 * t140 + t142 * t67) * t139;
t17 = t199 * t142 + (t65 + t67) * t140 + t197;
t15 = t28 * t140 - t27 * t142;
t14 = t26 * t140 - t25 * t142;
t8 = -t43 * t141 + (t140 * t27 + t142 * t28) * t139;
t7 = -t42 * t141 + (t140 * t25 + t142 * t26) * t139;
t1 = [Icges(4,2) * t150 ^ 2 + Icges(2,3) + Icges(3,3) + t83 + t85 + (Icges(4,1) * t149 + 0.2e1 * Icges(4,4) * t150) * t149 + (Icges(5,4) * t139 + Icges(5,2) * t141 - t101 - t93) * t141 + (Icges(5,1) * t139 - t180 + t192 - t196) * t139 + m(7) * (t44 ^ 2 + t45 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2) + m(5) * (t77 ^ 2 + t78 ^ 2) + m(4) * (t81 ^ 2 + t82 ^ 2) + m(3) * (t115 ^ 2 + t116 ^ 2) + m(2) * (t131 ^ 2 + t132 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t140 * t44 - t142 * t45) + m(6) * (t140 * t49 - t142 * t50) + m(5) * (t140 * t77 - t142 * t78) + m(4) * (t140 * t81 - t142 * t82); 0; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t177; m(7) * (t55 * t44 + t54 * t45) + m(6) * (t80 * t49 + t79 * t50) + (-t39 / 0.2e1 - t29 / 0.2e1 + t142 * t211 + (-Icges(5,6) * t142 + t161 * t140) * t207 - t77 * t204 + t119 * t205 - t175) * t142 + (t40 / 0.2e1 + t30 / 0.2e1 + t140 * t211 + t141 * (Icges(5,6) * t140 + t161 * t142) / 0.2e1 - t78 * t204 + t119 * t208 + t174) * t140; m(5) * t58 + m(6) * t35 + m(7) * t17; m(6) * (t80 * t140 - t79 * t142) + m(7) * (t55 * t140 - t54 * t142); m(7) * (t17 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t35 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t177 * t122 ^ 2 + t58 ^ 2) + (t137 * t88 + t13 + t15) * t140 + (-t138 * t87 - t12 - t14 + (-t140 * t87 + t142 * t88) * t140) * t142; (-t46 - t53) * t141 + m(7) * (t31 * t44 + t32 * t45) + m(6) * (t51 * t49 + t52 * t50) + (t175 * t140 + t174 * t142) * t139 + t171; m(6) * t41 + m(7) * t18; m(6) * (t51 * t140 - t52 * t142) + m(7) * (t31 * t140 - t32 * t142); t8 * t208 + t7 * t206 + (t34 * t140 - t33 * t142) * t207 + (t14 * t208 + t15 * t205) * t139 + m(7) * (t18 * t17 + t31 * t55 + t32 * t54) + m(6) * (t41 * t35 + t51 * t80 + t52 * t79) + t169; (t53 * t141 - t9) * t141 + m(7) * (t18 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t41 ^ 2 + t51 ^ 2 + t52 ^ 2) + (t140 * t7 + t142 * t8 - t141 * (t140 * t33 + t142 * t34)) * t139 + t209; -t195 + m(7) * (t47 * t44 + t48 * t45) + t171; m(7) * t38; m(7) * (t47 * t140 - t48 * t142); m(7) * (t38 * t17 + t47 * t55 + t48 * t54) + t169; m(7) * (t38 * t18 + t47 * t31 + t48 * t32) + t170; m(7) * (t38 ^ 2 + t47 ^ 2 + t48 ^ 2) + t170;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
