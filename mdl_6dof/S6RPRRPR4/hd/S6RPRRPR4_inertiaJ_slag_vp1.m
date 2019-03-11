% Calculate joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:04
% EndTime: 2019-03-09 05:09:10
% DurationCPUTime: 2.44s
% Computational Cost: add. (8186->378), mult. (6253->562), div. (0->0), fcn. (6528->12), ass. (0->181)
t175 = sin(qJ(1));
t167 = t175 ^ 2;
t171 = cos(pkin(11));
t154 = pkin(5) * t171 + pkin(4);
t173 = -pkin(9) - qJ(5);
t169 = sin(pkin(11));
t218 = t175 * t169;
t166 = pkin(10) + qJ(3);
t160 = qJ(4) + t166;
t153 = cos(t160);
t176 = cos(qJ(1));
t221 = t153 * t176;
t152 = sin(t160);
t222 = t152 * t176;
t165 = pkin(11) + qJ(6);
t156 = sin(t165);
t216 = t176 * t156;
t158 = cos(t165);
t219 = t175 * t158;
t119 = -t153 * t216 + t219;
t215 = t176 * t158;
t220 = t175 * t156;
t120 = t153 * t215 + t220;
t63 = rSges(7,1) * t120 + rSges(7,2) * t119 + rSges(7,3) * t222;
t246 = pkin(5) * t218 + t154 * t221 - t173 * t222 + t63;
t212 = qJ(5) + t173;
t236 = -pkin(4) + t154;
t89 = -t153 * rSges(7,3) + (rSges(7,1) * t158 - rSges(7,2) * t156) * t152;
t245 = -t152 * t236 - t153 * t212 - t89;
t187 = Icges(5,5) * t153 - Icges(5,6) * t152;
t100 = -Icges(5,3) * t176 + t175 * t187;
t101 = Icges(5,3) * t175 + t176 * t187;
t213 = t176 * t171;
t123 = -t153 * t218 - t213;
t214 = t176 * t169;
t217 = t175 * t171;
t124 = t153 * t217 - t214;
t168 = t176 ^ 2;
t224 = Icges(5,4) * t153;
t189 = -Icges(5,2) * t152 + t224;
t103 = Icges(5,6) * t175 + t176 * t189;
t225 = Icges(5,4) * t152;
t191 = Icges(5,1) * t153 - t225;
t105 = Icges(5,5) * t175 + t176 * t191;
t185 = -t103 * t152 + t105 * t153;
t102 = -Icges(5,6) * t176 + t175 * t189;
t104 = -Icges(5,5) * t176 + t175 * t191;
t186 = t102 * t152 - t104 * t153;
t223 = t152 * t175;
t70 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t223;
t125 = -t153 * t214 + t217;
t126 = t153 * t213 + t218;
t71 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t222;
t72 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t223;
t73 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t222;
t74 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t223;
t75 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t222;
t117 = -t153 * t220 - t215;
t118 = t153 * t219 - t216;
t56 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t223;
t58 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t223;
t60 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t223;
t16 = t117 * t58 + t118 * t60 + t223 * t56;
t57 = Icges(7,5) * t120 + Icges(7,6) * t119 + Icges(7,3) * t222;
t59 = Icges(7,4) * t120 + Icges(7,2) * t119 + Icges(7,6) * t222;
t61 = Icges(7,1) * t120 + Icges(7,4) * t119 + Icges(7,5) * t222;
t17 = t117 * t59 + t118 * t61 + t223 * t57;
t8 = -t16 * t176 + t17 * t175;
t244 = -t8 - t168 * t100 + (t123 * t72 + t124 * t74 + t223 * t70) * t176 + (-(-t101 + t186) * t176 - t123 * t73 - t124 * t75 - t223 * t71 - t185 * t175) * t175;
t243 = t153 ^ 2;
t242 = m(6) / 0.2e1;
t241 = m(7) / 0.2e1;
t240 = t175 / 0.2e1;
t239 = -t176 / 0.2e1;
t157 = sin(t166);
t238 = pkin(3) * t157;
t237 = pkin(4) * t153;
t174 = -pkin(7) - qJ(2);
t172 = cos(pkin(10));
t155 = pkin(2) * t172 + pkin(1);
t159 = cos(t166);
t140 = pkin(3) * t159 + t155;
t133 = t176 * t140;
t235 = t176 * (-t155 * t176 + t133) + (t140 - t155) * t167;
t181 = rSges(5,1) * t221 - rSges(5,2) * t222 + rSges(5,3) * t175;
t196 = rSges(5,1) * t153 - rSges(5,2) * t152;
t66 = t175 * (-t176 * rSges(5,3) + t175 * t196) + t176 * t181;
t234 = rSges(4,1) * t159;
t233 = rSges(4,2) * t157;
t87 = -Icges(7,6) * t153 + (Icges(7,4) * t158 - Icges(7,2) * t156) * t152;
t232 = t156 * t87;
t22 = -t153 * t56 + (-t156 * t58 + t158 * t60) * t152;
t231 = t22 * t176;
t23 = -t153 * t57 + (-t156 * t59 + t158 * t61) * t152;
t230 = t23 * t175;
t229 = rSges(3,3) + qJ(2);
t131 = t152 * pkin(4) - t153 * qJ(5);
t93 = -t153 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t169) * t152;
t228 = -t131 - t93;
t227 = Icges(4,4) * t157;
t226 = Icges(4,4) * t159;
t210 = pkin(4) * t221 + qJ(5) * t222;
t211 = t167 * (qJ(5) * t152 + t237) + t176 * t210;
t209 = rSges(4,3) * t175 + t176 * t234;
t207 = t167 + t168;
t18 = t119 * t58 + t120 * t60 + t222 * t56;
t19 = t119 * t59 + t120 * t61 + t222 * t57;
t9 = t175 * t19 - t176 * t18;
t206 = (t9 + t167 * t101 + (t125 * t73 + t126 * t75 + t222 * t71) * t175 + (-t125 * t72 - t126 * t74 - t222 * t70 + (-t100 + t185) * t175 + t186 * t176) * t176) * t175;
t205 = t242 + t241;
t204 = -t131 + t245;
t203 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t222;
t202 = -t131 - t238;
t132 = rSges(5,1) * t152 + rSges(5,2) * t153;
t201 = -t132 - t238;
t164 = -pkin(8) + t174;
t200 = -t164 * t175 + t133;
t195 = -t124 * rSges(6,1) - t123 * rSges(6,2);
t30 = t175 * (rSges(6,3) * t223 - t195) + t176 * t203 + t211;
t86 = -Icges(7,3) * t153 + (Icges(7,5) * t158 - Icges(7,6) * t156) * t152;
t88 = -Icges(7,5) * t153 + (Icges(7,1) * t158 - Icges(7,4) * t156) * t152;
t31 = t117 * t87 + t118 * t88 + t223 * t86;
t3 = -t31 * t153 + (t16 * t175 + t17 * t176) * t152;
t32 = t119 * t87 + t120 * t88 + t222 * t86;
t4 = -t32 * t153 + (t175 * t18 + t176 * t19) * t152;
t199 = t3 * t239 + t4 * t240 - t153 * (t230 - t231) / 0.2e1 + t8 * t223 / 0.2e1 + t9 * t222 / 0.2e1;
t198 = t202 - t93;
t197 = -t233 + t234;
t194 = -t118 * rSges(7,1) - t117 * rSges(7,2);
t193 = t202 + t245;
t192 = Icges(4,1) * t159 - t227;
t190 = -Icges(4,2) * t157 + t226;
t188 = Icges(4,5) * t159 - Icges(4,6) * t157;
t129 = Icges(5,2) * t153 + t225;
t130 = Icges(5,1) * t152 + t224;
t182 = -t129 * t152 + t130 * t153;
t62 = rSges(7,3) * t223 - t194;
t14 = t211 + (-t210 + t246) * t176 + (-pkin(5) * t214 + t62 + (-t152 * t212 + t153 * t236) * t175) * t175;
t170 = sin(pkin(10));
t180 = rSges(3,1) * t172 - rSges(3,2) * t170 + pkin(1);
t178 = t176 * t244 + t206;
t128 = Icges(5,5) * t152 + Icges(5,6) * t153;
t90 = -Icges(6,3) * t153 + (Icges(6,5) * t171 - Icges(6,6) * t169) * t152;
t91 = -Icges(6,6) * t153 + (Icges(6,4) * t171 - Icges(6,2) * t169) * t152;
t92 = -Icges(6,5) * t153 + (Icges(6,1) * t171 - Icges(6,4) * t169) * t152;
t177 = -t231 / 0.2e1 + t230 / 0.2e1 + (t125 * t91 + t126 * t92 + t175 * t128 + t176 * t182 + t222 * t90 + t32 + (t103 - t71) * t153 + (-t169 * t73 + t171 * t75 + t105) * t152) * t240 + (t123 * t91 + t124 * t92 - t176 * t128 + t175 * t182 + t223 * t90 + t31 + (-t70 + t102) * t153 + (-t169 * t72 + t171 * t74 + t104) * t152) * t239;
t147 = rSges(2,1) * t176 - rSges(2,2) * t175;
t146 = -rSges(2,1) * t175 - rSges(2,2) * t176;
t139 = rSges(4,1) * t157 + rSges(4,2) * t159;
t112 = Icges(4,3) * t175 + t176 * t188;
t111 = -Icges(4,3) * t176 + t175 * t188;
t99 = t175 * t229 + t176 * t180;
t98 = -t175 * t180 + t176 * t229;
t95 = t201 * t176;
t94 = t201 * t175;
t84 = -t175 * t174 + (t155 - t233) * t176 + t209;
t83 = (rSges(4,3) - t174) * t176 + (-t155 - t197) * t175;
t80 = t152 * t158 * t88;
t79 = t181 + t200;
t78 = (rSges(5,3) - t164) * t176 + (-t140 - t196) * t175;
t77 = t228 * t176;
t76 = t228 * t175;
t69 = t176 * (-t176 * t233 + t209) + (-t176 * rSges(4,3) + t175 * t197) * t175;
t55 = t198 * t176;
t54 = t198 * t175;
t47 = t200 + t203 + t210;
t46 = -t176 * t164 + (-t237 - t140 + (-rSges(6,3) - qJ(5)) * t152) * t175 + t195;
t45 = t204 * t176;
t44 = t204 * t175;
t43 = t193 * t176;
t42 = t193 * t175;
t41 = -t153 * t63 - t222 * t89;
t40 = t153 * t62 + t223 * t89;
t39 = t200 + t246;
t38 = (pkin(5) * t169 - t164) * t176 + (-t153 * t154 - t140 + (-rSges(7,3) + t173) * t152) * t175 + t194;
t37 = -t152 * t232 - t153 * t86 + t80;
t36 = (-t175 * t63 + t176 * t62) * t152;
t35 = t66 + t235;
t15 = t30 + t235;
t13 = t14 + t235;
t1 = [Icges(3,2) * t172 ^ 2 + t159 * (Icges(4,2) * t159 + t227) + t157 * (Icges(4,1) * t157 + t226) + Icges(2,3) + t80 + (Icges(3,1) * t170 + 0.2e1 * Icges(3,4) * t172) * t170 + (-t86 + t129 - t90) * t153 + (-t169 * t91 + t171 * t92 + t130 - t232) * t152 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t83 ^ 2 + t84 ^ 2) + m(3) * (t98 ^ 2 + t99 ^ 2) + m(2) * (t146 ^ 2 + t147 ^ 2); m(7) * (t175 * t38 - t176 * t39) + m(6) * (t175 * t46 - t176 * t47) + m(5) * (t175 * t78 - t176 * t79) + m(4) * (t175 * t83 - t176 * t84) + m(3) * (t175 * t98 - t176 * t99); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + t205) * t207; m(7) * (t38 * t43 + t39 * t42) + m(6) * (t46 * t55 + t47 * t54) + m(5) * (t78 * t95 + t79 * t94) + m(4) * (-t175 * t84 - t176 * t83) * t139 + t177 + (t167 / 0.2e1 + t168 / 0.2e1) * (Icges(4,5) * t157 + Icges(4,6) * t159) + (t159 * (Icges(4,6) * t175 + t176 * t190) + t157 * (Icges(4,5) * t175 + t176 * t192)) * t240 + (t159 * (-Icges(4,6) * t176 + t175 * t190) + t157 * (-Icges(4,5) * t176 + t175 * t192)) * t239; m(5) * (t175 * t95 - t176 * t94) + m(6) * (t175 * t55 - t176 * t54) + m(7) * (t175 * t43 - t176 * t42); m(7) * (t13 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t15 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t35 ^ 2 + t94 ^ 2 + t95 ^ 2) + t175 * t167 * t112 + m(4) * (t139 ^ 2 * t207 + t69 ^ 2) + t206 + (-t168 * t111 + (-t175 * t111 + t176 * t112) * t175 + t244) * t176; m(7) * (t38 * t45 + t39 * t44) + m(6) * (t46 * t77 + t47 * t76) + m(5) * (-t175 * t79 - t176 * t78) * t132 + t177; m(6) * (t175 * t77 - t176 * t76) + m(7) * (t175 * t45 - t176 * t44); m(7) * (t13 * t14 + t42 * t44 + t43 * t45) + m(6) * (t15 * t30 + t54 * t76 + t55 * t77) + m(5) * (t66 * t35 + (-t175 * t94 - t176 * t95) * t132) + t178; m(7) * (t14 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t30 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t132 ^ 2 * t207 + t66 ^ 2) + t178; 0.2e1 * ((t175 * t39 + t176 * t38) * t241 + (t175 * t47 + t176 * t46) * t242) * t152; 0; m(7) * (-t153 * t13 + (t175 * t42 + t176 * t43) * t152) + m(6) * (-t153 * t15 + (t175 * t54 + t176 * t55) * t152); m(7) * (-t153 * t14 + (t175 * t44 + t176 * t45) * t152) + m(6) * (-t153 * t30 + (t175 * t76 + t176 * t77) * t152); 0.2e1 * t205 * (t152 ^ 2 * t207 + t243); m(7) * (t38 * t40 + t39 * t41) - t37 * t153 + ((t23 / 0.2e1 + t32 / 0.2e1) * t176 + (t22 / 0.2e1 + t31 / 0.2e1) * t175) * t152; m(7) * (t175 * t40 - t176 * t41); m(7) * (t13 * t36 + t40 * t43 + t41 * t42) + t199; m(7) * (t14 * t36 + t40 * t45 + t41 * t44) + t199; m(7) * (-t36 * t153 + (t175 * t41 + t176 * t40) * t152); t243 * t37 + m(7) * (t36 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t176 * t4 + t175 * t3 - t153 * (t175 * t22 + t176 * t23)) * t152;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
