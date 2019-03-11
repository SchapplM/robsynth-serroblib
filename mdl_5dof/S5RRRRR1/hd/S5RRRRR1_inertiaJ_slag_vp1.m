% Calculate joint inertia matrix for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:56
% EndTime: 2019-03-08 18:37:00
% DurationCPUTime: 1.83s
% Computational Cost: add. (6929->328), mult. (6617->484), div. (0->0), fcn. (6954->10), ass. (0->186)
t158 = qJ(2) + qJ(3);
t154 = qJ(4) + t158;
t149 = sin(t154);
t150 = cos(t154);
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t77 = t150 * rSges(6,3) + (-rSges(6,1) * t162 + rSges(6,2) * t159) * t149;
t236 = -t149 * pkin(4) + t150 * pkin(6) + t77;
t161 = sin(qJ(1));
t156 = t161 ^ 2;
t164 = cos(qJ(1));
t152 = sin(t158);
t153 = cos(t158);
t211 = Icges(4,4) * t153;
t179 = -Icges(4,2) * t152 + t211;
t95 = -Icges(4,6) * t161 + t164 * t179;
t212 = Icges(4,4) * t152;
t182 = Icges(4,1) * t153 - t212;
t97 = -Icges(4,5) * t161 + t164 * t182;
t184 = t152 * t95 - t153 * t97;
t94 = Icges(4,6) * t164 + t161 * t179;
t96 = Icges(4,5) * t164 + t161 * t182;
t185 = -t152 * t94 + t153 * t96;
t209 = Icges(5,4) * t150;
t178 = -Icges(5,2) * t149 + t209;
t83 = -Icges(5,6) * t161 + t164 * t178;
t210 = Icges(5,4) * t149;
t181 = Icges(5,1) * t150 - t210;
t85 = -Icges(5,5) * t161 + t164 * t181;
t186 = t149 * t83 - t150 * t85;
t82 = Icges(5,6) * t164 + t161 * t178;
t84 = Icges(5,5) * t164 + t161 * t181;
t187 = -t149 * t82 + t150 * t84;
t175 = Icges(5,5) * t150 - Icges(5,6) * t149;
t80 = Icges(5,3) * t164 + t161 * t175;
t81 = -Icges(5,3) * t161 + t164 * t175;
t202 = t164 * t159;
t203 = t161 * t162;
t118 = -t150 * t202 - t203;
t201 = t164 * t162;
t204 = t161 * t159;
t119 = t150 * t201 - t204;
t207 = t149 * t164;
t116 = -t150 * t204 + t201;
t117 = t150 * t203 + t202;
t208 = t149 * t161;
t53 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t208;
t55 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t208;
t57 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t208;
t17 = t118 * t55 + t119 * t57 + t207 * t53;
t54 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t207;
t56 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t207;
t58 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t207;
t18 = t118 * t56 + t119 * t58 + t207 * t54;
t9 = -t18 * t161 + t17 * t164;
t231 = -t156 * t81 - (t187 * t164 + (t186 - t80) * t161) * t164 - t9;
t176 = Icges(4,5) * t153 - Icges(4,6) * t152;
t92 = Icges(4,3) * t164 + t161 * t176;
t93 = -Icges(4,3) * t161 + t164 * t176;
t235 = -t156 * t93 - (t185 * t164 + (t184 - t92) * t161) * t164 + t231;
t157 = t164 ^ 2;
t234 = -t161 / 0.2e1;
t233 = t164 / 0.2e1;
t15 = t116 * t55 + t117 * t57 + t208 * t53;
t16 = t116 * t56 + t117 * t58 + t208 * t54;
t8 = t15 * t164 - t16 * t161;
t232 = (t157 * t80 + t8 + (t186 * t161 + (t187 - t81) * t164) * t161) * t164;
t160 = sin(qJ(2));
t230 = pkin(2) * t160;
t229 = pkin(3) * t152;
t228 = pkin(4) * t150;
t163 = cos(qJ(2));
t151 = t163 * pkin(2) + pkin(1);
t74 = Icges(6,3) * t150 + (-Icges(6,5) * t162 + Icges(6,6) * t159) * t149;
t75 = Icges(6,6) * t150 + (-Icges(6,4) * t162 + Icges(6,2) * t159) * t149;
t227 = t149 * t159 * t75 + t150 * t74;
t189 = rSges(5,1) * t150 - rSges(5,2) * t149;
t219 = t164 * rSges(5,3);
t86 = t161 * t189 + t219;
t134 = pkin(3) * t153 + t151;
t88 = (t134 - t151) * t161;
t226 = -t86 - t88;
t206 = t150 * t164;
t87 = rSges(5,1) * t206 - rSges(5,2) * t207 - t161 * rSges(5,3);
t126 = t164 * t134;
t146 = t164 * t151;
t89 = -t146 + t126;
t225 = -t87 - t89;
t224 = rSges(3,1) * t163;
t223 = rSges(4,1) * t153;
t222 = rSges(3,2) * t160;
t76 = Icges(6,5) * t150 + (-Icges(6,1) * t162 + Icges(6,4) * t159) * t149;
t221 = t162 * t76;
t220 = t164 * rSges(4,3);
t23 = t150 * t53 + (t159 * t55 - t162 * t57) * t149;
t218 = t23 * t164;
t24 = t150 * t54 + (t159 * t56 - t162 * t58) * t149;
t217 = t24 * t161;
t188 = -t117 * rSges(6,1) - t116 * rSges(6,2);
t61 = rSges(6,3) * t208 - t188;
t216 = -(pkin(6) * t149 + t228) * t161 - t61;
t62 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t207;
t215 = -pkin(4) * t206 - pkin(6) * t207 - t62;
t59 = t236 * t161;
t60 = t236 * t164;
t214 = Icges(3,4) * t160;
t213 = Icges(3,4) * t163;
t205 = t152 * t164;
t200 = t156 + t157;
t199 = t164 * (t157 * t92 + (t184 * t161 + (t185 - t93) * t164) * t161) + t232;
t198 = t161 * t229;
t197 = pkin(3) * t205;
t196 = -t88 + t216;
t195 = -t89 + t215;
t133 = -t152 * rSges(4,1) - t153 * rSges(4,2);
t194 = t133 - t230;
t29 = t116 * t75 + t117 * t76 + t208 * t74;
t3 = t29 * t150 + (t15 * t161 + t16 * t164) * t149;
t30 = t118 * t75 + t119 * t76 + t207 * t74;
t4 = t30 * t150 + (t161 * t17 + t164 * t18) * t149;
t193 = t3 * t233 + t4 * t234 + t150 * (-t217 + t218) / 0.2e1 + t8 * t208 / 0.2e1 + t9 * t207 / 0.2e1;
t192 = -t229 - t230;
t191 = t222 - t224;
t190 = -rSges(4,2) * t152 + t223;
t183 = Icges(3,1) * t163 - t214;
t180 = -Icges(3,2) * t160 + t213;
t177 = Icges(3,5) * t163 - Icges(3,6) * t160;
t122 = -Icges(5,2) * t150 - t210;
t123 = -Icges(5,1) * t149 - t209;
t172 = -t122 * t149 + t123 * t150;
t131 = -Icges(4,2) * t153 - t212;
t132 = -Icges(4,1) * t152 - t211;
t171 = -t131 * t152 + t132 * t153;
t170 = t161 * t231 + t232;
t169 = t192 * t161;
t168 = t192 * t164;
t121 = -Icges(5,5) * t149 - Icges(5,6) * t150;
t167 = t218 / 0.2e1 - t217 / 0.2e1 + (-t161 * t121 - t149 * t85 - t150 * t83 + t164 * t172 + t30) * t234 + (t164 * t121 - t149 * t84 - t150 * t82 + t161 * t172 + t29) * t233;
t101 = -rSges(4,2) * t205 - t161 * rSges(4,3) + t164 * t223;
t166 = t235 * t161 + t199;
t130 = -Icges(4,5) * t152 - Icges(4,6) * t153;
t165 = t167 + (-t161 * t130 - t152 * t97 - t153 * t95 + t164 * t171) * t234 + (t164 * t130 - t152 * t96 - t153 * t94 + t161 * t171) * t233;
t148 = t164 * t224;
t145 = t164 * rSges(2,1) - t161 * rSges(2,2);
t144 = -t161 * rSges(2,1) - t164 * rSges(2,2);
t143 = -t160 * rSges(3,1) - t163 * rSges(3,2);
t128 = -t164 * pkin(1) + t146;
t127 = (-pkin(1) + t151) * t161;
t124 = -t149 * rSges(5,1) - t150 * rSges(5,2);
t113 = t164 * t124;
t112 = t161 * t124;
t107 = -Icges(3,3) * t161 + t164 * t177;
t106 = Icges(3,3) * t164 + t161 * t177;
t100 = t161 * t190 + t220;
t99 = -t161 * rSges(3,3) + t148 + (pkin(1) - t222) * t164;
t98 = -t164 * rSges(3,3) + (-pkin(1) + t191) * t161;
t91 = t194 * t164;
t90 = t194 * t161;
t79 = t113 - t197;
t78 = t112 - t198;
t71 = t101 + t146;
t70 = -t220 + (-t151 - t190) * t161;
t69 = t113 + t168;
t68 = t112 + t169;
t65 = t126 + t87;
t64 = -t219 + (-t134 - t189) * t161;
t63 = -t164 * (-t164 * t222 + t148) + t191 * t156;
t52 = -t161 * t100 - t164 * t101;
t51 = -t161 * t86 - t164 * t87;
t50 = -t197 + t60;
t49 = -t198 + t59;
t44 = t168 + t60;
t43 = t169 + t59;
t38 = t126 - t215;
t37 = (-t228 - t134 + (-rSges(6,3) - pkin(6)) * t149) * t161 + t188;
t36 = (-t101 - t128) * t164 + (-t100 - t127) * t161;
t35 = t150 * t62 - t207 * t77;
t34 = -t150 * t61 + t208 * t77;
t33 = t161 * t226 + t164 * t225;
t32 = (-t149 * t221 + t227) * t150;
t31 = (-t161 * t62 + t164 * t61) * t149;
t26 = (-t128 + t225) * t164 + (-t127 + t226) * t161;
t25 = t161 * t216 + t164 * t215;
t12 = t161 * t196 + t164 * t195;
t11 = (-t128 + t195) * t164 + (-t127 + t196) * t161;
t1 = [-t150 * t122 - t153 * t131 - t152 * t132 - t163 * (-Icges(3,2) * t163 - t214) - t160 * (-Icges(3,1) * t160 - t213) + Icges(2,3) + (-t123 - t221) * t149 + m(6) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) + m(4) * (t70 ^ 2 + t71 ^ 2) + m(3) * (t98 ^ 2 + t99 ^ 2) + m(2) * (t144 ^ 2 + t145 ^ 2) + t227; (-t163 * (Icges(3,6) * t164 + t161 * t180) - t160 * (Icges(3,5) * t164 + t161 * t183)) * t233 + (-t163 * (-Icges(3,6) * t161 + t164 * t180) - t160 * (-Icges(3,5) * t161 + t164 * t183)) * t234 + t165 + (t157 / 0.2e1 + t156 / 0.2e1) * (-Icges(3,5) * t160 - Icges(3,6) * t163) + m(6) * (t44 * t37 + t43 * t38) + m(5) * (t69 * t64 + t68 * t65) + m(4) * (t91 * t70 + t90 * t71) + m(3) * (t161 * t99 + t164 * t98) * t143; m(6) * (t11 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t26 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t36 ^ 2 + t90 ^ 2 + t91 ^ 2) + t164 * t157 * t106 + m(3) * (t143 ^ 2 * t200 + t63 ^ 2) + t199 + (-t156 * t107 + (t161 * t106 - t164 * t107) * t164 + t235) * t161; t165 + m(4) * (t161 * t71 + t164 * t70) * t133 + m(6) * (t50 * t37 + t49 * t38) + m(5) * (t79 * t64 + t78 * t65); m(6) * (t12 * t11 + t49 * t43 + t50 * t44) + m(5) * (t33 * t26 + t78 * t68 + t79 * t69) + m(4) * (t52 * t36 + (t161 * t90 + t164 * t91) * t133) + t166; m(6) * (t12 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t33 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t133 ^ 2 * t200 + t52 ^ 2) + t166; m(6) * (t60 * t37 + t59 * t38) + m(5) * (t161 * t65 + t164 * t64) * t124 + t167; m(6) * (t25 * t11 + t59 * t43 + t60 * t44) + m(5) * (t51 * t26 + (t161 * t68 + t164 * t69) * t124) + t170; m(6) * (t25 * t12 + t59 * t49 + t60 * t50) + m(5) * (t51 * t33 + (t161 * t78 + t164 * t79) * t124) + t170; m(5) * (t124 ^ 2 * t200 + t51 ^ 2) + m(6) * (t25 ^ 2 + t59 ^ 2 + t60 ^ 2) + t170; m(6) * (t34 * t37 + t35 * t38) + t32 + ((t24 / 0.2e1 + t30 / 0.2e1) * t164 + (t29 / 0.2e1 + t23 / 0.2e1) * t161) * t149; m(6) * (t31 * t11 + t34 * t44 + t35 * t43) + t193; m(6) * (t31 * t12 + t34 * t50 + t35 * t49) + t193; m(6) * (t31 * t25 + t34 * t60 + t35 * t59) + t193; m(6) * (t31 ^ 2 + t34 ^ 2 + t35 ^ 2) + t150 * t32 + (t164 * t4 + t161 * t3 + t150 * (t161 * t23 + t164 * t24)) * t149;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11); t1(2) t1(3) t1(5) t1(8) t1(12); t1(4) t1(5) t1(6) t1(9) t1(13); t1(7) t1(8) t1(9) t1(10) t1(14); t1(11) t1(12) t1(13) t1(14) t1(15);];
Mq  = res;
