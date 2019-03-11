% Calculate joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:26
% EndTime: 2019-03-08 18:45:35
% DurationCPUTime: 4.10s
% Computational Cost: add. (30038->463), mult. (78359->675), div. (0->0), fcn. (103509->16), ass. (0->215)
t208 = m(6) / 0.2e1 + m(7) / 0.2e1;
t235 = 0.2e1 * t208;
t185 = sin(pkin(11));
t188 = cos(pkin(11));
t189 = cos(pkin(6));
t217 = sin(pkin(12));
t201 = t189 * t217;
t219 = cos(pkin(12));
t174 = t185 * t219 + t188 * t201;
t192 = sin(qJ(3));
t202 = t189 * t219;
t196 = t185 * t217 - t188 * t202;
t220 = cos(pkin(7));
t194 = t196 * t220;
t186 = sin(pkin(6));
t218 = sin(pkin(7));
t203 = t186 * t218;
t228 = cos(qJ(3));
t158 = t174 * t228 + (-t188 * t203 - t194) * t192;
t204 = t186 * t220;
t167 = -t188 * t204 + t196 * t218;
t191 = sin(qJ(4));
t227 = cos(qJ(4));
t148 = t158 * t191 - t167 * t227;
t175 = -t185 * t201 + t188 * t219;
t195 = t185 * t202 + t188 * t217;
t193 = t195 * t220;
t160 = t175 * t228 + (t185 * t203 - t193) * t192;
t168 = t185 * t204 + t195 * t218;
t150 = t160 * t191 - t168 * t227;
t198 = t220 * t219;
t166 = t189 * t218 * t192 + (t192 * t198 + t228 * t217) * t186;
t173 = t189 * t220 - t219 * t203;
t161 = t166 * t191 - t173 * t227;
t149 = t158 * t227 + t167 * t191;
t199 = t228 * t218;
t197 = t186 * t199;
t157 = t174 * t192 + t188 * t197 + t228 * t194;
t183 = pkin(13) + qJ(6);
t180 = sin(t183);
t181 = cos(t183);
t117 = -t149 * t180 + t157 * t181;
t118 = t149 * t181 + t157 * t180;
t72 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t148;
t74 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t148;
t76 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t148;
t28 = t117 * t74 + t118 * t76 + t148 * t72;
t151 = t160 * t227 + t168 * t191;
t159 = t175 * t192 - t185 * t197 + t228 * t193;
t119 = -t151 * t180 + t159 * t181;
t120 = t151 * t181 + t159 * t180;
t73 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t150;
t75 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t150;
t77 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t150;
t29 = t117 * t75 + t118 * t77 + t148 * t73;
t162 = t166 * t227 + t173 * t191;
t165 = -t189 * t199 + (t192 * t217 - t198 * t228) * t186;
t143 = -t162 * t180 + t165 * t181;
t144 = t162 * t181 + t165 * t180;
t93 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t161;
t94 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t161;
t95 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t161;
t44 = t117 * t94 + t118 * t95 + t148 * t93;
t1 = t28 * t148 + t29 * t150 + t44 * t161;
t234 = t1 / 0.2e1;
t30 = t119 * t74 + t120 * t76 + t150 * t72;
t31 = t119 * t75 + t120 * t77 + t150 * t73;
t45 = t119 * t94 + t120 * t95 + t150 * t93;
t2 = t30 * t148 + t31 * t150 + t45 * t161;
t233 = t2 / 0.2e1;
t37 = t143 * t74 + t144 * t76 + t161 * t72;
t38 = t143 * t75 + t144 * t77 + t161 * t73;
t51 = t143 * t94 + t144 * t95 + t161 * t93;
t11 = t37 * t148 + t38 * t150 + t51 * t161;
t232 = t11 / 0.2e1;
t231 = t148 / 0.2e1;
t230 = t150 / 0.2e1;
t229 = t161 / 0.2e1;
t187 = cos(pkin(13));
t226 = t187 * pkin(5);
t184 = sin(pkin(13));
t216 = t157 * t184;
t78 = t118 * rSges(7,1) + t117 * rSges(7,2) + t148 * rSges(7,3);
t225 = pkin(5) * t216 + pkin(10) * t148 + t226 * t149 + t78;
t215 = t159 * t184;
t79 = t120 * rSges(7,1) + t119 * rSges(7,2) + t150 * rSges(7,3);
t224 = pkin(5) * t215 + pkin(10) * t150 + t226 * t151 + t79;
t214 = t165 * t184;
t96 = t144 * rSges(7,1) + t143 * rSges(7,2) + t161 * rSges(7,3);
t223 = pkin(5) * t214 + pkin(10) * t161 + t226 * t162 + t96;
t114 = t149 * pkin(4) + t148 * qJ(5);
t121 = -t149 * t184 + t157 * t187;
t122 = t149 * t187 + t216;
t86 = t122 * rSges(6,1) + t121 * rSges(6,2) + t148 * rSges(6,3);
t222 = -t114 - t86;
t115 = t151 * pkin(4) + t150 * qJ(5);
t123 = -t151 * t184 + t159 * t187;
t124 = t151 * t187 + t215;
t87 = t124 * rSges(6,1) + t123 * rSges(6,2) + t150 * rSges(6,3);
t221 = -t115 - t87;
t146 = -t162 * t184 + t165 * t187;
t147 = t162 * t187 + t214;
t106 = t147 * rSges(6,1) + t146 * rSges(6,2) + t161 * rSges(6,3);
t142 = t162 * pkin(4) + t161 * qJ(5);
t212 = -t106 - t142;
t140 = t158 * pkin(3) + t157 * pkin(9);
t137 = t168 * t140;
t211 = t168 * t114 + t137;
t141 = t160 * pkin(3) + t159 * pkin(9);
t139 = t173 * t141;
t210 = t173 * t115 + t139;
t156 = t166 * pkin(3) + t165 * pkin(9);
t145 = t167 * t156;
t209 = t167 * t142 + t145;
t207 = -t114 - t225;
t206 = -t115 - t224;
t205 = -t142 - t223;
t200 = m(3) + m(4) + m(5) + m(6) + m(7);
t155 = t166 * rSges(4,1) - t165 * rSges(4,2) + t173 * rSges(4,3);
t154 = Icges(4,1) * t166 - Icges(4,4) * t165 + Icges(4,5) * t173;
t153 = Icges(4,4) * t166 - Icges(4,2) * t165 + Icges(4,6) * t173;
t152 = Icges(4,5) * t166 - Icges(4,6) * t165 + Icges(4,3) * t173;
t136 = t162 * rSges(5,1) - t161 * rSges(5,2) + t165 * rSges(5,3);
t135 = Icges(5,1) * t162 - Icges(5,4) * t161 + Icges(5,5) * t165;
t134 = Icges(5,4) * t162 - Icges(5,2) * t161 + Icges(5,6) * t165;
t133 = Icges(5,5) * t162 - Icges(5,6) * t161 + Icges(5,3) * t165;
t132 = t160 * rSges(4,1) - t159 * rSges(4,2) + t168 * rSges(4,3);
t131 = t158 * rSges(4,1) - t157 * rSges(4,2) + t167 * rSges(4,3);
t130 = Icges(4,1) * t160 - Icges(4,4) * t159 + Icges(4,5) * t168;
t129 = Icges(4,1) * t158 - Icges(4,4) * t157 + Icges(4,5) * t167;
t128 = Icges(4,4) * t160 - Icges(4,2) * t159 + Icges(4,6) * t168;
t127 = Icges(4,4) * t158 - Icges(4,2) * t157 + Icges(4,6) * t167;
t126 = Icges(4,5) * t160 - Icges(4,6) * t159 + Icges(4,3) * t168;
t125 = Icges(4,5) * t158 - Icges(4,6) * t157 + Icges(4,3) * t167;
t116 = t157 * t142;
t111 = t165 * t115;
t109 = t159 * t114;
t108 = t151 * rSges(5,1) - t150 * rSges(5,2) + t159 * rSges(5,3);
t107 = t149 * rSges(5,1) - t148 * rSges(5,2) + t157 * rSges(5,3);
t105 = Icges(5,1) * t151 - Icges(5,4) * t150 + Icges(5,5) * t159;
t104 = Icges(5,1) * t149 - Icges(5,4) * t148 + Icges(5,5) * t157;
t103 = Icges(5,4) * t151 - Icges(5,2) * t150 + Icges(5,6) * t159;
t102 = Icges(5,4) * t149 - Icges(5,2) * t148 + Icges(5,6) * t157;
t101 = Icges(5,5) * t151 - Icges(5,6) * t150 + Icges(5,3) * t159;
t100 = Icges(5,5) * t149 - Icges(5,6) * t148 + Icges(5,3) * t157;
t99 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t161;
t98 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t161;
t97 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t161;
t90 = t173 * t132 - t168 * t155;
t89 = -t173 * t131 + t167 * t155;
t88 = t168 * t131 - t167 * t132;
t85 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t150;
t84 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t148;
t83 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t150;
t82 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t148;
t81 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t150;
t80 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t148;
t69 = t165 * t108 - t159 * t136;
t68 = -t165 * t107 + t157 * t136;
t67 = t165 * t133 - t161 * t134 + t162 * t135;
t66 = t159 * t107 - t157 * t108;
t65 = t173 * t108 + t139 + (-t136 - t156) * t168;
t64 = t167 * t136 + t145 + (-t107 - t140) * t173;
t63 = t159 * t133 - t150 * t134 + t151 * t135;
t62 = t157 * t133 - t148 * t134 + t149 * t135;
t61 = -t150 * t96 + t161 * t79;
t60 = t148 * t96 - t161 * t78;
t59 = t168 * t107 + t137 + (-t108 - t141) * t167;
t58 = t165 * t101 - t161 * t103 + t162 * t105;
t57 = t165 * t100 - t161 * t102 + t162 * t104;
t56 = t159 * t101 - t150 * t103 + t151 * t105;
t55 = t159 * t100 - t150 * t102 + t151 * t104;
t54 = t157 * t101 - t148 * t103 + t149 * t105;
t53 = t157 * t100 - t148 * t102 + t149 * t104;
t52 = t146 * t98 + t147 * t99 + t161 * t97;
t50 = -t148 * t79 + t150 * t78;
t49 = t212 * t159 + t165 * t87 + t111;
t48 = t157 * t106 + t222 * t165 + t116;
t47 = t123 * t98 + t124 * t99 + t150 * t97;
t46 = t121 * t98 + t122 * t99 + t148 * t97;
t43 = t173 * t87 + (-t156 + t212) * t168 + t210;
t42 = t167 * t106 + (-t140 + t222) * t173 + t209;
t41 = t221 * t157 + t159 * t86 + t109;
t40 = t146 * t83 + t147 * t85 + t161 * t81;
t39 = t146 * t82 + t147 * t84 + t161 * t80;
t36 = t168 * t86 + (-t141 + t221) * t167 + t211;
t35 = t123 * t83 + t124 * t85 + t150 * t81;
t34 = t123 * t82 + t124 * t84 + t150 * t80;
t33 = t121 * t83 + t122 * t85 + t148 * t81;
t32 = t121 * t82 + t122 * t84 + t148 * t80;
t27 = t205 * t159 + t224 * t165 + t111;
t26 = t223 * t157 + t207 * t165 + t116;
t25 = t224 * t173 + (-t156 + t205) * t168 + t210;
t24 = t223 * t167 + (-t140 + t207) * t173 + t209;
t23 = t57 * t167 + t58 * t168 + t67 * t173;
t22 = t206 * t157 + t225 * t159 + t109;
t21 = t57 * t157 + t58 * t159 + t67 * t165;
t20 = t225 * t168 + (-t141 + t206) * t167 + t211;
t19 = t55 * t167 + t56 * t168 + t63 * t173;
t18 = t53 * t167 + t54 * t168 + t62 * t173;
t17 = t55 * t157 + t56 * t159 + t63 * t165;
t16 = t53 * t157 + t54 * t159 + t62 * t165;
t15 = t39 * t167 + t40 * t168 + t52 * t173;
t14 = t39 * t157 + t40 * t159 + t52 * t165;
t13 = t37 * t167 + t38 * t168 + t51 * t173;
t12 = t37 * t157 + t38 * t159 + t51 * t165;
t10 = t34 * t167 + t35 * t168 + t47 * t173;
t9 = t32 * t167 + t33 * t168 + t46 * t173;
t8 = t34 * t157 + t35 * t159 + t47 * t165;
t7 = t32 * t157 + t33 * t159 + t46 * t165;
t6 = t30 * t167 + t31 * t168 + t45 * t173;
t5 = t28 * t167 + t29 * t168 + t44 * t173;
t4 = t30 * t157 + t31 * t159 + t45 * t165;
t3 = t28 * t157 + t29 * t159 + t44 * t165;
t70 = [m(2) + t200; t200 * t189; 0.2e1 * (m(5) / 0.2e1 + m(4) / 0.2e1 + m(3) / 0.2e1 + t208) * (t189 ^ 2 + (t185 ^ 2 + t188 ^ 2) * t186 ^ 2); m(4) * t88 + m(5) * t59 + m(6) * t36 + m(7) * t20; m(6) * (t36 * t189 + (t185 * t42 - t188 * t43) * t186) + m(7) * (t20 * t189 + (t185 * t24 - t188 * t25) * t186) + m(5) * (t59 * t189 + (t185 * t64 - t188 * t65) * t186) + m(4) * (t88 * t189 + (t185 * t89 - t188 * t90) * t186); (t13 + t15 + t23 + (t173 * t152 - t165 * t153 + t166 * t154) * t173) * t173 + (t6 + t19 + t10 + (t168 * t126 - t159 * t128 + t160 * t130) * t168 + (t173 * t126 - t165 * t128 + t166 * t130 + t168 * t152 - t159 * t153 + t160 * t154) * t173) * t168 + (t5 + t9 + t18 + (t167 * t125 - t157 * t127 + t158 * t129) * t167 + (t173 * t125 - t165 * t127 + t166 * t129 + t167 * t152 - t157 * t153 + t158 * t154) * t173 + (t168 * t125 + t167 * t126 - t159 * t127 - t157 * t128 + t160 * t129 + t158 * t130) * t168) * t167 + m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t59 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2); m(5) * t66 + m(6) * t41 + m(7) * t22; m(6) * (t41 * t189 + (t185 * t48 - t188 * t49) * t186) + m(7) * (t22 * t189 + (t185 * t26 - t188 * t27) * t186) + m(5) * (t66 * t189 + (t185 * t68 - t188 * t69) * t186); (t12 / 0.2e1 + t14 / 0.2e1 + t21 / 0.2e1) * t173 + (t4 / 0.2e1 + t8 / 0.2e1 + t17 / 0.2e1) * t168 + (t3 / 0.2e1 + t7 / 0.2e1 + t16 / 0.2e1) * t167 + (t13 / 0.2e1 + t15 / 0.2e1 + t23 / 0.2e1) * t165 + (t6 / 0.2e1 + t19 / 0.2e1 + t10 / 0.2e1) * t159 + (t5 / 0.2e1 + t9 / 0.2e1 + t18 / 0.2e1) * t157 + m(7) * (t20 * t22 + t24 * t26 + t25 * t27) + m(6) * (t36 * t41 + t42 * t48 + t43 * t49) + m(5) * (t66 * t59 + t68 * t64 + t69 * t65); (t12 + t14 + t21) * t165 + (t4 + t17 + t8) * t159 + (t3 + t16 + t7) * t157 + m(7) * (t22 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t66 ^ 2 + t68 ^ 2 + t69 ^ 2); t161 * t235; (t161 * t189 + (-t148 * t188 + t150 * t185) * t186) * t235; m(7) * (t148 * t25 + t150 * t24 + t161 * t20) + m(6) * (t148 * t43 + t150 * t42 + t161 * t36); m(7) * (t148 * t27 + t150 * t26 + t161 * t22) + m(6) * (t148 * t49 + t150 * t48 + t161 * t41); (t148 ^ 2 + t150 ^ 2 + t161 ^ 2) * t235; m(7) * t50; m(7) * (t50 * t189 + (t185 * t60 - t188 * t61) * t186); m(7) * (t50 * t20 + t60 * t24 + t61 * t25) + t168 * t233 + t167 * t234 + t13 * t229 + t173 * t232 + t5 * t231 + t6 * t230; m(7) * (t50 * t22 + t60 * t26 + t61 * t27) + t12 * t229 + t157 * t234 + t165 * t232 + t159 * t233 + t3 * t231 + t4 * t230; m(7) * (t61 * t148 + t60 * t150 + t50 * t161); t148 * t1 + t161 * t11 + t150 * t2 + m(7) * (t50 ^ 2 + t60 ^ 2 + t61 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t70(1) t70(2) t70(4) t70(7) t70(11) t70(16); t70(2) t70(3) t70(5) t70(8) t70(12) t70(17); t70(4) t70(5) t70(6) t70(9) t70(13) t70(18); t70(7) t70(8) t70(9) t70(10) t70(14) t70(19); t70(11) t70(12) t70(13) t70(14) t70(15) t70(20); t70(16) t70(17) t70(18) t70(19) t70(20) t70(21);];
Mq  = res;
