% Calculate joint inertia matrix for
% S6RPRPRP2
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
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:41
% EndTime: 2019-03-09 03:04:45
% DurationCPUTime: 1.99s
% Computational Cost: add. (5922->333), mult. (5080->480), div. (0->0), fcn. (5398->10), ass. (0->160)
t224 = Icges(4,3) + Icges(5,3);
t140 = qJ(3) + pkin(10);
t135 = sin(t140);
t137 = cos(t140);
t144 = sin(qJ(3));
t147 = cos(qJ(3));
t223 = Icges(4,5) * t147 + Icges(5,5) * t137 - Icges(4,6) * t144 - Icges(5,6) * t135;
t141 = qJ(1) + pkin(9);
t138 = cos(t141);
t146 = cos(qJ(5));
t178 = t138 * t146;
t136 = sin(t141);
t143 = sin(qJ(5));
t183 = t136 * t143;
t101 = t137 * t183 + t178;
t179 = t138 * t143;
t182 = t136 * t146;
t102 = t137 * t182 - t179;
t219 = rSges(7,3) + qJ(6);
t221 = rSges(7,1) + pkin(5);
t222 = -t219 * t101 - t221 * t102;
t83 = -Icges(6,3) * t137 + (Icges(6,5) * t146 - Icges(6,6) * t143) * t135;
t86 = -Icges(7,2) * t137 + (Icges(7,4) * t146 + Icges(7,6) * t143) * t135;
t220 = -t83 - t86;
t187 = t135 * t136;
t45 = Icges(7,5) * t102 + Icges(7,6) * t187 + Icges(7,3) * t101;
t49 = Icges(7,4) * t102 + Icges(7,2) * t187 + Icges(7,6) * t101;
t53 = Icges(7,1) * t102 + Icges(7,4) * t187 + Icges(7,5) * t101;
t12 = t101 * t45 + t102 * t53 + t187 * t49;
t103 = t137 * t179 - t182;
t104 = t137 * t178 + t183;
t186 = t135 * t138;
t46 = Icges(7,5) * t104 + Icges(7,6) * t186 + Icges(7,3) * t103;
t50 = Icges(7,4) * t104 + Icges(7,2) * t186 + Icges(7,6) * t103;
t54 = Icges(7,1) * t104 + Icges(7,4) * t186 + Icges(7,5) * t103;
t13 = t101 * t46 + t102 * t54 + t187 * t50;
t47 = Icges(6,5) * t102 - Icges(6,6) * t101 + Icges(6,3) * t187;
t51 = Icges(6,4) * t102 - Icges(6,2) * t101 + Icges(6,6) * t187;
t55 = Icges(6,1) * t102 - Icges(6,4) * t101 + Icges(6,5) * t187;
t14 = -t101 * t51 + t102 * t55 + t187 * t47;
t48 = Icges(6,5) * t104 - Icges(6,6) * t103 + Icges(6,3) * t186;
t52 = Icges(6,4) * t104 - Icges(6,2) * t103 + Icges(6,6) * t186;
t56 = Icges(6,1) * t104 - Icges(6,4) * t103 + Icges(6,5) * t186;
t15 = -t101 * t52 + t102 * t56 + t187 * t48;
t82 = -Icges(7,6) * t137 + (Icges(7,5) * t146 + Icges(7,3) * t143) * t135;
t90 = -Icges(7,4) * t137 + (Icges(7,1) * t146 + Icges(7,5) * t143) * t135;
t30 = t101 * t82 + t102 * t90 + t187 * t86;
t87 = -Icges(6,6) * t137 + (Icges(6,4) * t146 - Icges(6,2) * t143) * t135;
t91 = -Icges(6,5) * t137 + (Icges(6,1) * t146 - Icges(6,4) * t143) * t135;
t31 = -t101 * t87 + t102 * t91 + t187 * t83;
t218 = (-t30 - t31) * t137 + ((t13 + t15) * t138 + (t12 + t14) * t136) * t135;
t16 = t103 * t45 + t104 * t53 + t186 * t49;
t17 = t103 * t46 + t104 * t54 + t186 * t50;
t18 = -t103 * t51 + t104 * t55 + t186 * t47;
t19 = -t103 * t52 + t104 * t56 + t186 * t48;
t32 = t103 * t82 + t104 * t90 + t186 * t86;
t33 = -t103 * t87 + t104 * t91 + t186 * t83;
t217 = (-t32 - t33) * t137 + ((t17 + t19) * t138 + (t16 + t18) * t136) * t135;
t216 = t135 / 0.2e1;
t215 = t137 / 0.2e1;
t214 = t144 / 0.2e1;
t213 = t147 / 0.2e1;
t20 = -t137 * t49 + (t143 * t45 + t146 * t53) * t135;
t22 = -t137 * t47 + (-t143 * t51 + t146 * t55) * t135;
t212 = -t20 - t22;
t21 = -t137 * t50 + (t143 * t46 + t146 * t54) * t135;
t23 = -t137 * t48 + (-t143 * t52 + t146 * t56) * t135;
t211 = t21 + t23;
t210 = t224 * t136 + t223 * t138;
t209 = -t223 * t136 + t224 * t138;
t185 = t135 * t143;
t208 = t82 * t185 + (t90 + t91) * t135 * t146;
t133 = t136 ^ 2;
t134 = t138 ^ 2;
t206 = -t137 / 0.2e1;
t124 = rSges(4,1) * t144 + rSges(4,2) * t147;
t204 = m(4) * t124;
t145 = sin(qJ(1));
t203 = pkin(1) * t145;
t202 = pkin(3) * t144;
t201 = pkin(4) * t137;
t200 = t220 * t137 - t185 * t87 + t208;
t199 = rSges(7,2) * t187 - t222;
t198 = rSges(7,2) * t186 + t219 * t103 + t221 * t104;
t132 = pkin(3) * t147 + pkin(2);
t119 = t138 * t132;
t131 = t138 * pkin(7);
t142 = -qJ(4) - pkin(7);
t180 = t138 * t142;
t197 = t136 * (t180 + t131 + (-pkin(2) + t132) * t136) + t138 * (-pkin(2) * t138 + t119 + (-pkin(7) - t142) * t136);
t195 = rSges(4,1) * t147;
t194 = rSges(4,2) * t144;
t193 = t138 * rSges(4,3);
t192 = -rSges(7,2) * t137 + (t219 * t143 + t221 * t146) * t135;
t191 = Icges(4,4) * t144;
t190 = Icges(4,4) * t147;
t189 = Icges(5,4) * t135;
t188 = Icges(5,4) * t137;
t181 = t137 * t138;
t177 = pkin(4) * t181 + pkin(8) * t186;
t176 = t136 * rSges(4,3) + t138 * t195;
t175 = t133 + t134;
t60 = t104 * rSges(6,1) - t103 * rSges(6,2) + rSges(6,3) * t186;
t174 = Icges(4,5) * t214 + Icges(5,5) * t216 + Icges(4,6) * t213 + Icges(5,6) * t215;
t173 = -rSges(5,1) * t135 - rSges(5,2) * t137 - t202;
t172 = -pkin(4) * t135 + pkin(8) * t137 - t202;
t171 = -t132 - t201;
t170 = t133 * (pkin(8) * t135 + t201) + t138 * t177 + t197;
t95 = -rSges(6,3) * t137 + (rSges(6,1) * t146 - rSges(6,2) * t143) * t135;
t169 = t172 - t95;
t148 = cos(qJ(1));
t139 = t148 * pkin(1);
t168 = -t136 * t142 + t119 + t139;
t167 = -t194 + t195;
t166 = rSges(5,1) * t137 - rSges(5,2) * t135;
t165 = -rSges(6,1) * t102 + rSges(6,2) * t101;
t164 = -t180 - t203;
t159 = Icges(4,1) * t147 - t191;
t158 = Icges(5,1) * t137 - t189;
t157 = -Icges(4,2) * t144 + t190;
t156 = -Icges(5,2) * t135 + t188;
t153 = t172 - t192;
t152 = rSges(5,1) * t181 - rSges(5,2) * t186 + t136 * rSges(5,3);
t151 = t31 / 0.2e1 + t30 / 0.2e1 + t20 / 0.2e1 + t22 / 0.2e1;
t150 = t33 / 0.2e1 + t32 / 0.2e1 + t23 / 0.2e1 + t21 / 0.2e1;
t149 = t168 + t177;
t126 = rSges(2,1) * t148 - t145 * rSges(2,2);
t125 = -t145 * rSges(2,1) - rSges(2,2) * t148;
t107 = rSges(3,1) * t138 - rSges(3,2) * t136 + t139;
t106 = -rSges(3,1) * t136 - rSges(3,2) * t138 - t203;
t79 = t173 * t138;
t78 = t173 * t136;
t66 = pkin(7) * t136 + t139 + (pkin(2) - t194) * t138 + t176;
t65 = t193 - t203 + t131 + (-pkin(2) - t167) * t136;
t62 = t152 + t168;
t61 = -t203 + (rSges(5,3) - t142) * t138 + (-t132 - t166) * t136;
t58 = rSges(6,3) * t187 - t165;
t44 = t169 * t138;
t43 = t169 * t136;
t42 = t138 * (-t138 * t194 + t176) + (t136 * t167 - t193) * t136;
t41 = t153 * t138;
t40 = t153 * t136;
t37 = -t137 * t60 - t186 * t95;
t36 = t137 * t58 + t187 * t95;
t35 = t149 + t60;
t34 = ((-rSges(6,3) - pkin(8)) * t135 + t171) * t136 + t164 + t165;
t29 = (-t136 * t60 + t138 * t58) * t135;
t28 = t138 * t152 + (-rSges(5,3) * t138 + t136 * t166) * t136 + t197;
t27 = t149 + t198;
t26 = ((-rSges(7,2) - pkin(8)) * t135 + t171) * t136 + t164 + t222;
t25 = -t198 * t137 - t192 * t186;
t24 = t199 * t137 + t192 * t187;
t11 = (-t198 * t136 + t199 * t138) * t135;
t10 = t136 * t58 + t138 * t60 + t170;
t9 = t199 * t136 + t198 * t138 + t170;
t8 = t136 * t19 - t138 * t18;
t7 = t136 * t17 - t138 * t16;
t6 = t136 * t15 - t138 * t14;
t5 = -t12 * t138 + t13 * t136;
t1 = [t147 * (Icges(4,2) * t147 + t191) + t144 * (Icges(4,1) * t144 + t190) + Icges(2,3) + Icges(3,3) + (Icges(5,1) * t135 - t143 * t87 + t188) * t135 + (Icges(5,2) * t137 + t189 + t220) * t137 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(3) * (t106 ^ 2 + t107 ^ 2) + m(2) * (t125 ^ 2 + t126 ^ 2) + t208; 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t26 * t41 + t27 * t40) + m(6) * (t34 * t44 + t35 * t43) + m(5) * (t61 * t79 + t62 * t78) + (-t65 * t204 - t144 * (-Icges(4,5) * t138 + t136 * t159) / 0.2e1 - t147 * (-Icges(4,6) * t138 + t136 * t157) / 0.2e1 - t135 * (-Icges(5,5) * t138 + t136 * t158) / 0.2e1 + (-Icges(5,6) * t138 + t136 * t156) * t206 + t174 * t138 - t151) * t138 + (t174 * t136 - t66 * t204 + (Icges(4,6) * t136 + t138 * t157) * t213 + (Icges(4,5) * t136 + t138 * t159) * t214 + (Icges(5,6) * t136 + t138 * t156) * t215 + (Icges(5,5) * t136 + t138 * t158) * t216 + t150) * t136; m(4) * t42 + m(5) * t28 + m(6) * t10 + m(7) * t9; m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t28 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t124 ^ 2 * t175 + t42 ^ 2) + (t209 * t134 - t5 - t6) * t138 + (t7 + t8 + t210 * t133 + (t209 * t136 + t210 * t138) * t138) * t136; m(7) * (t136 * t26 - t138 * t27) + m(6) * (t136 * t34 - t138 * t35) + m(5) * (t136 * t61 - t138 * t62); 0; m(7) * (t136 * t41 - t138 * t40) + m(6) * (t136 * t44 - t138 * t43) + m(5) * (t136 * t79 - t138 * t78); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t175; -t200 * t137 + m(7) * (t24 * t26 + t25 * t27) + m(6) * (t34 * t36 + t35 * t37) + (t136 * t151 + t138 * t150) * t135; m(6) * t29 + m(7) * t11; m(7) * (t11 * t9 + t24 * t41 + t25 * t40) + m(6) * (t10 * t29 + t36 * t44 + t37 * t43) + ((t8 / 0.2e1 + t7 / 0.2e1) * t138 + (t6 / 0.2e1 + t5 / 0.2e1) * t136) * t135 + t217 * t136 / 0.2e1 + (t211 * t136 + t212 * t138) * t206 - t218 * t138 / 0.2e1; m(6) * (t136 * t36 - t138 * t37) + m(7) * (t136 * t24 - t138 * t25); m(7) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t29 ^ 2 + t36 ^ 2 + t37 ^ 2) + t200 * t137 ^ 2 + ((-t211 * t137 + t217) * t138 + (t212 * t137 + t218) * t136) * t135; m(7) * (t101 * t27 + t103 * t26); m(7) * t185; m(7) * (t101 * t40 + t103 * t41 + t185 * t9); m(7) * (-t101 * t138 + t103 * t136); m(7) * (t101 * t25 + t103 * t24 + t11 * t185); m(7) * (t135 ^ 2 * t143 ^ 2 + t101 ^ 2 + t103 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
