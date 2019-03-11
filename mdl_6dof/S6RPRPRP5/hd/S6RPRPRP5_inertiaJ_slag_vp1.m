% Calculate joint inertia matrix for
% S6RPRPRP5
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:33
% EndTime: 2019-03-09 03:14:40
% DurationCPUTime: 2.70s
% Computational Cost: add. (6270->382), mult. (5991->577), div. (0->0), fcn. (6369->10), ass. (0->176)
t150 = pkin(9) + qJ(3);
t147 = cos(t150);
t149 = pkin(10) + qJ(5);
t146 = cos(t149);
t160 = cos(qJ(1));
t190 = t160 * t146;
t144 = sin(t149);
t159 = sin(qJ(1));
t195 = t159 * t144;
t110 = t147 * t195 + t190;
t191 = t160 * t144;
t194 = t159 * t146;
t111 = t147 * t194 - t191;
t227 = rSges(7,3) + qJ(6);
t229 = rSges(7,1) + pkin(5);
t231 = -t110 * t227 - t111 * t229;
t216 = t159 / 0.2e1;
t230 = t160 / 0.2e1;
t145 = sin(t150);
t82 = -Icges(6,3) * t147 + (Icges(6,5) * t146 - Icges(6,6) * t144) * t145;
t83 = -Icges(7,2) * t147 + (Icges(7,4) * t146 + Icges(7,6) * t144) * t145;
t228 = -t82 - t83;
t198 = t145 * t159;
t46 = Icges(7,5) * t111 + Icges(7,6) * t198 + Icges(7,3) * t110;
t50 = Icges(7,4) * t111 + Icges(7,2) * t198 + Icges(7,6) * t110;
t54 = Icges(7,1) * t111 + Icges(7,4) * t198 + Icges(7,5) * t110;
t12 = t110 * t46 + t111 * t54 + t198 * t50;
t112 = t147 * t191 - t194;
t113 = t147 * t190 + t195;
t197 = t145 * t160;
t47 = Icges(7,5) * t113 + Icges(7,6) * t197 + Icges(7,3) * t112;
t51 = Icges(7,4) * t113 + Icges(7,2) * t197 + Icges(7,6) * t112;
t55 = Icges(7,1) * t113 + Icges(7,4) * t197 + Icges(7,5) * t112;
t13 = t110 * t47 + t111 * t55 + t198 * t51;
t48 = Icges(6,5) * t111 - Icges(6,6) * t110 + Icges(6,3) * t198;
t52 = Icges(6,4) * t111 - Icges(6,2) * t110 + Icges(6,6) * t198;
t56 = Icges(6,1) * t111 - Icges(6,4) * t110 + Icges(6,5) * t198;
t14 = -t110 * t52 + t111 * t56 + t198 * t48;
t49 = Icges(6,5) * t113 - Icges(6,6) * t112 + Icges(6,3) * t197;
t53 = Icges(6,4) * t113 - Icges(6,2) * t112 + Icges(6,6) * t197;
t57 = Icges(6,1) * t113 - Icges(6,4) * t112 + Icges(6,5) * t197;
t15 = -t110 * t53 + t111 * t57 + t198 * t49;
t81 = -Icges(7,6) * t147 + (Icges(7,5) * t146 + Icges(7,3) * t144) * t145;
t85 = -Icges(7,4) * t147 + (Icges(7,1) * t146 + Icges(7,5) * t144) * t145;
t29 = t110 * t81 + t111 * t85 + t198 * t83;
t84 = -Icges(6,6) * t147 + (Icges(6,4) * t146 - Icges(6,2) * t144) * t145;
t86 = -Icges(6,5) * t147 + (Icges(6,1) * t146 - Icges(6,4) * t144) * t145;
t30 = -t110 * t84 + t111 * t86 + t198 * t82;
t226 = (-t29 - t30) * t147 + ((t13 + t15) * t160 + (t12 + t14) * t159) * t145;
t16 = t112 * t46 + t113 * t54 + t197 * t50;
t17 = t112 * t47 + t113 * t55 + t197 * t51;
t18 = -t112 * t52 + t113 * t56 + t197 * t48;
t19 = -t112 * t53 + t113 * t57 + t197 * t49;
t31 = t112 * t81 + t113 * t85 + t197 * t83;
t32 = -t112 * t84 + t113 * t86 + t197 * t82;
t225 = (-t31 - t32) * t147 + ((t17 + t19) * t160 + (t16 + t18) * t159) * t145;
t20 = -t147 * t50 + (t144 * t46 + t146 * t54) * t145;
t22 = -t147 * t48 + (-t144 * t52 + t146 * t56) * t145;
t224 = -t20 - t22;
t21 = -t147 * t51 + (t144 * t47 + t146 * t55) * t145;
t23 = -t147 * t49 + (-t144 * t53 + t146 * t57) * t145;
t223 = t21 + t23;
t200 = t144 * t145;
t222 = t81 * t200 + (t85 + t86) * t145 * t146;
t221 = t147 ^ 2;
t151 = t159 ^ 2;
t152 = t160 ^ 2;
t220 = m(5) / 0.2e1;
t219 = m(6) / 0.2e1;
t218 = m(7) / 0.2e1;
t215 = -t160 / 0.2e1;
t126 = rSges(4,1) * t145 + rSges(4,2) * t147;
t214 = m(4) * t126;
t213 = pkin(3) * t147;
t155 = cos(pkin(10));
t141 = pkin(4) * t155 + pkin(3);
t212 = -pkin(3) + t141;
t211 = t147 * t228 - t200 * t84 + t222;
t210 = rSges(7,2) * t198 - t231;
t209 = rSges(7,2) * t197 + t112 * t227 + t113 * t229;
t207 = rSges(3,3) + qJ(2);
t206 = -t147 * rSges(7,2) + (t144 * t227 + t146 * t229) * t145;
t125 = t145 * pkin(3) - t147 * qJ(4);
t157 = -pkin(8) - qJ(4);
t205 = -t125 - (qJ(4) + t157) * t147 - t212 * t145;
t153 = sin(pkin(10));
t204 = -t125 + t147 * rSges(5,3) - (rSges(5,1) * t155 - rSges(5,2) * t153) * t145;
t203 = Icges(4,4) * t145;
t202 = Icges(4,4) * t147;
t201 = qJ(4) * t145;
t196 = t147 * t160;
t193 = t159 * t153;
t192 = t159 * t155;
t189 = t160 * t153;
t188 = t160 * t155;
t158 = -pkin(7) - qJ(2);
t187 = t160 * t158;
t184 = pkin(3) * t196 + qJ(4) * t197;
t186 = t151 * (t201 + t213) + t160 * t184;
t185 = -pkin(4) * t189 - t157 * t198;
t183 = t151 + t152;
t88 = -t147 * rSges(6,3) + (rSges(6,1) * t146 - rSges(6,2) * t144) * t145;
t182 = -t88 + t205;
t61 = rSges(6,1) * t113 - rSges(6,2) * t112 + rSges(6,3) * t197;
t118 = -t147 * t189 + t192;
t119 = t147 * t188 + t193;
t181 = rSges(5,1) * t119 + rSges(5,2) * t118 + rSges(5,3) * t197;
t156 = cos(pkin(9));
t142 = pkin(2) * t156 + pkin(1);
t180 = t142 * t160 - t158 * t159;
t179 = -t141 * t147 - t142;
t178 = t205 - t206;
t164 = pkin(4) * t193 + t141 * t196 - t157 * t197;
t177 = t159 * ((t147 * t212 - t201) * t159 + t185) + t160 * (t164 - t184) + t186;
t176 = t220 + t219 + t218;
t175 = -t185 - t187;
t174 = rSges(4,1) * t147 - rSges(4,2) * t145;
t116 = -t147 * t193 - t188;
t117 = t147 * t192 - t189;
t173 = -t117 * rSges(5,1) - t116 * rSges(5,2);
t172 = -t111 * rSges(6,1) + t110 * rSges(6,2);
t169 = Icges(4,1) * t147 - t203;
t168 = -Icges(4,2) * t145 + t202;
t167 = Icges(4,5) * t147 - Icges(4,6) * t145;
t166 = rSges(4,1) * t196 - rSges(4,2) * t197 + rSges(4,3) * t159;
t154 = sin(pkin(9));
t165 = rSges(3,1) * t156 - rSges(3,2) * t154 + pkin(1);
t163 = t22 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1 + t20 / 0.2e1;
t162 = t23 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1 + t21 / 0.2e1;
t161 = t164 + t180;
t143 = t145 ^ 2;
t129 = rSges(2,1) * t160 - rSges(2,2) * t159;
t128 = -rSges(2,1) * t159 - rSges(2,2) * t160;
t121 = Icges(4,5) * t145 + Icges(4,6) * t147;
t96 = Icges(4,3) * t159 + t160 * t167;
t95 = -Icges(4,3) * t160 + t159 * t167;
t93 = -Icges(5,5) * t147 + (Icges(5,1) * t155 - Icges(5,4) * t153) * t145;
t92 = -Icges(5,6) * t147 + (Icges(5,4) * t155 - Icges(5,2) * t153) * t145;
t90 = t159 * t207 + t160 * t165;
t89 = -t159 * t165 + t160 * t207;
t79 = t166 + t180;
t78 = (rSges(4,3) - t158) * t160 + (-t142 - t174) * t159;
t74 = t204 * t160;
t73 = t204 * t159;
t70 = Icges(5,1) * t119 + Icges(5,4) * t118 + Icges(5,5) * t197;
t69 = Icges(5,1) * t117 + Icges(5,4) * t116 + Icges(5,5) * t198;
t68 = Icges(5,4) * t119 + Icges(5,2) * t118 + Icges(5,6) * t197;
t67 = Icges(5,4) * t117 + Icges(5,2) * t116 + Icges(5,6) * t198;
t66 = Icges(5,5) * t119 + Icges(5,6) * t118 + Icges(5,3) * t197;
t65 = Icges(5,5) * t117 + Icges(5,6) * t116 + Icges(5,3) * t198;
t64 = t160 * t166 + (-t160 * rSges(4,3) + t159 * t174) * t159;
t59 = rSges(6,3) * t198 - t172;
t45 = t180 + t181 + t184;
t44 = -t187 + (-t213 - t142 + (-rSges(5,3) - qJ(4)) * t145) * t159 + t173;
t43 = t182 * t160;
t42 = t182 * t159;
t41 = -t147 * t61 - t197 * t88;
t40 = t147 * t59 + t198 * t88;
t39 = t161 + t61;
t38 = (-rSges(6,3) * t145 + t179) * t159 + t172 + t175;
t37 = t178 * t160;
t36 = t178 * t159;
t33 = (-t159 * t61 + t160 * t59) * t145;
t28 = t159 * (rSges(5,3) * t198 - t173) + t160 * t181 + t186;
t27 = t161 + t209;
t26 = (-rSges(7,2) * t145 + t179) * t159 + t175 + t231;
t25 = -t147 * t209 - t197 * t206;
t24 = t147 * t210 + t198 * t206;
t11 = (-t159 * t209 + t160 * t210) * t145;
t10 = t159 * t59 + t160 * t61 + t177;
t9 = t159 * t210 + t160 * t209 + t177;
t8 = t159 * t19 - t160 * t18;
t7 = t159 * t17 - t16 * t160;
t6 = -t14 * t160 + t15 * t159;
t5 = -t12 * t160 + t13 * t159;
t1 = [Icges(3,2) * t156 ^ 2 + Icges(2,3) + (Icges(3,1) * t154 + 0.2e1 * Icges(3,4) * t156) * t154 + (t203 - (Icges(5,5) * t155 - Icges(5,6) * t153) * t145 + (Icges(4,2) + Icges(5,3)) * t147 + t228) * t147 + (Icges(4,1) * t145 - t144 * t84 - t153 * t92 + t155 * t93 + t202) * t145 + m(6) * (t38 ^ 2 + t39 ^ 2) + m(7) * (t26 ^ 2 + t27 ^ 2) + m(5) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t89 ^ 2 + t90 ^ 2) + m(2) * (t128 ^ 2 + t129 ^ 2) + t222; m(6) * (t159 * t38 - t160 * t39) + m(7) * (t159 * t26 - t160 * t27) + m(5) * (t159 * t44 - t160 * t45) + m(4) * (t159 * t78 - t160 * t79) + m(3) * (t159 * t89 - t160 * t90); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t176) * t183; m(6) * (t38 * t43 + t39 * t42) + m(7) * (t26 * t37 + t27 * t36) + m(5) * (t44 * t74 + t45 * t73) + (-t116 * t92 / 0.2e1 - t117 * t93 / 0.2e1 - t78 * t214 + t121 * t230 + (t65 / 0.2e1 + Icges(4,6) * t230 - t159 * t168 / 0.2e1) * t147 - t163) * t160 + (t118 * t92 / 0.2e1 + t119 * t93 / 0.2e1 - t79 * t214 + t121 * t216 + (Icges(4,6) * t216 + t168 * t230 - t66 / 0.2e1) * t147 + t162) * t159 + ((Icges(4,5) * t159 - t153 * t68 + t155 * t70 + t160 * t169) * t216 + (-Icges(4,5) * t160 - t153 * t67 + t155 * t69 + t159 * t169) * t215) * t145; m(5) * (t159 * t74 - t160 * t73) + m(6) * (t159 * t43 - t160 * t42) + m(7) * (t159 * t37 - t160 * t36); m(7) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t28 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t126 ^ 2 * t183 + t64 ^ 2) + (-t152 * t95 - t5 - t6 + (t116 * t67 + t117 * t69 + t65 * t198) * t160) * t160 + (t8 + t7 + t151 * t96 + (t118 * t68 + t119 * t70 + t66 * t197) * t159 + (-t116 * t68 - t117 * t70 - t118 * t67 - t119 * t69 - t159 * t95 + t160 * t96 - t197 * t65 - t198 * t66) * t160) * t159; 0.2e1 * ((t159 * t39 + t160 * t38) * t219 + (t159 * t27 + t160 * t26) * t218 + (t159 * t45 + t160 * t44) * t220) * t145; 0; m(7) * (-t147 * t9 + (t159 * t36 + t160 * t37) * t145) + m(6) * (-t147 * t10 + (t159 * t42 + t160 * t43) * t145) + m(5) * (-t147 * t28 + (t159 * t73 + t160 * t74) * t145); 0.2e1 * t176 * (t143 * t183 + t221); -t211 * t147 + m(6) * (t38 * t40 + t39 * t41) + m(7) * (t24 * t26 + t25 * t27) + (t159 * t163 + t160 * t162) * t145; m(6) * (t159 * t40 - t160 * t41) + m(7) * (t159 * t24 - t160 * t25); m(7) * (t11 * t9 + t24 * t37 + t25 * t36) + m(6) * (t10 * t33 + t40 * t43 + t41 * t42) + ((t8 / 0.2e1 + t7 / 0.2e1) * t160 + (t6 / 0.2e1 + t5 / 0.2e1) * t159) * t145 - (t159 * t223 + t160 * t224) * t147 / 0.2e1 + t225 * t216 + t226 * t215; m(6) * (-t33 * t147 + (t159 * t41 + t160 * t40) * t145) + m(7) * (-t11 * t147 + (t159 * t25 + t160 * t24) * t145); m(7) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + t211 * t221 + (t225 * t160 + t226 * t159 + (t159 * t224 - t160 * t223) * t147) * t145; m(7) * (t110 * t27 + t112 * t26); m(7) * (-t110 * t160 + t112 * t159); m(7) * (t110 * t36 + t112 * t37 + t200 * t9); m(7) * (t110 * t159 + t112 * t160 - t144 * t147) * t145; m(7) * (t11 * t200 + t110 * t25 + t112 * t24); m(7) * (t143 * t144 ^ 2 + t110 ^ 2 + t112 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
