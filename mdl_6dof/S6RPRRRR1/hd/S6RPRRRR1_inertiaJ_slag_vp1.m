% Calculate joint inertia matrix for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:10
% EndTime: 2019-03-09 06:54:15
% DurationCPUTime: 1.96s
% Computational Cost: add. (10876->353), mult. (7129->504), div. (0->0), fcn. (7416->12), ass. (0->184)
t162 = qJ(1) + pkin(11);
t155 = sin(t162);
t153 = t155 ^ 2;
t246 = t155 * pkin(7);
t156 = cos(t162);
t154 = t156 ^ 2;
t164 = qJ(3) + qJ(4);
t157 = sin(t164);
t158 = cos(t164);
t223 = Icges(5,4) * t158;
t188 = -Icges(5,2) * t157 + t223;
t96 = Icges(5,6) * t155 + t188 * t156;
t224 = Icges(5,4) * t157;
t191 = Icges(5,1) * t158 - t224;
t98 = Icges(5,5) * t155 + t191 * t156;
t193 = -t157 * t96 + t158 * t98;
t95 = -Icges(5,6) * t156 + t188 * t155;
t97 = -Icges(5,5) * t156 + t191 * t155;
t194 = t157 * t95 - t158 * t97;
t159 = qJ(5) + t164;
t150 = sin(t159);
t151 = cos(t159);
t221 = Icges(6,4) * t151;
t187 = -Icges(6,2) * t150 + t221;
t88 = Icges(6,6) * t155 + t187 * t156;
t222 = Icges(6,4) * t150;
t190 = Icges(6,1) * t151 - t222;
t90 = Icges(6,5) * t155 + t190 * t156;
t195 = -t150 * t88 + t151 * t90;
t87 = -Icges(6,6) * t156 + t187 * t155;
t89 = -Icges(6,5) * t156 + t190 * t155;
t196 = t150 * t87 - t151 * t89;
t168 = cos(qJ(6));
t214 = t156 * t168;
t165 = sin(qJ(6));
t217 = t155 * t165;
t113 = -t151 * t217 - t214;
t215 = t156 * t165;
t216 = t155 * t168;
t114 = t151 * t216 - t215;
t220 = t150 * t155;
t56 = Icges(7,5) * t114 + Icges(7,6) * t113 + Icges(7,3) * t220;
t58 = Icges(7,4) * t114 + Icges(7,2) * t113 + Icges(7,6) * t220;
t60 = Icges(7,1) * t114 + Icges(7,4) * t113 + Icges(7,5) * t220;
t17 = t113 * t58 + t114 * t60 + t56 * t220;
t115 = -t151 * t215 + t216;
t116 = t151 * t214 + t217;
t219 = t150 * t156;
t57 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t219;
t59 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t219;
t61 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t219;
t18 = t113 * t59 + t114 * t61 + t57 * t220;
t8 = t18 * t155 - t17 * t156;
t184 = Icges(6,5) * t151 - Icges(6,6) * t150;
t85 = -Icges(6,3) * t156 + t184 * t155;
t86 = Icges(6,3) * t155 + t184 * t156;
t240 = -t154 * t85 - (t195 * t155 + (t196 - t86) * t156) * t155 - t8;
t185 = Icges(5,5) * t158 - Icges(5,6) * t157;
t93 = -Icges(5,3) * t156 + t185 * t155;
t94 = Icges(5,3) * t155 + t185 * t156;
t245 = -t154 * t93 - (t193 * t155 + (t194 - t94) * t156) * t155 + t240;
t218 = t151 * t156;
t65 = t116 * rSges(7,1) + t115 * rSges(7,2) + rSges(7,3) * t219;
t244 = pkin(5) * t218 + pkin(10) * t219 + t65;
t199 = rSges(5,1) * t158 - rSges(5,2) * t157;
t171 = -pkin(8) - pkin(7);
t243 = t155 / 0.2e1;
t242 = -t156 / 0.2e1;
t19 = t115 * t58 + t116 * t60 + t56 * t219;
t20 = t115 * t59 + t116 * t61 + t57 * t219;
t9 = t20 * t155 - t19 * t156;
t241 = (t153 * t86 + t9 + (t196 * t156 + (t195 - t85) * t155) * t156) * t155;
t166 = sin(qJ(3));
t239 = pkin(3) * t166;
t238 = pkin(4) * t157;
t237 = pkin(5) * t151;
t167 = sin(qJ(1));
t236 = t167 * pkin(1);
t169 = cos(qJ(3));
t152 = t169 * pkin(3) + pkin(2);
t135 = pkin(4) * t158 + t152;
t119 = t156 * t135;
t136 = t156 * t152;
t235 = t156 * (t119 - t136) + (t135 - t152) * t153;
t149 = t156 * pkin(7);
t234 = t155 * (t149 + (-pkin(2) + t152) * t155) + t156 * (-t156 * pkin(2) + t136 - t246);
t177 = rSges(6,1) * t218 - rSges(6,2) * t219 + t155 * rSges(6,3);
t198 = rSges(6,1) * t151 - rSges(6,2) * t150;
t45 = t155 * (-t156 * rSges(6,3) + t198 * t155) + t156 * t177;
t178 = t155 * rSges(5,3) + t199 * t156;
t50 = t155 * (-t156 * rSges(5,3) + t199 * t155) + t156 * t178;
t233 = rSges(4,1) * t169;
t231 = rSges(4,2) * t166;
t229 = t156 * rSges(4,3);
t24 = -t151 * t56 + (-t165 * t58 + t168 * t60) * t150;
t228 = t24 * t156;
t25 = -t151 * t57 + (-t165 * t59 + t168 * t61) * t150;
t227 = t25 * t155;
t226 = Icges(4,4) * t166;
t225 = Icges(4,4) * t169;
t100 = -Icges(7,6) * t151 + (Icges(7,4) * t168 - Icges(7,2) * t165) * t150;
t213 = t165 * t100;
t104 = -t151 * rSges(7,3) + (rSges(7,1) * t168 - rSges(7,2) * t165) * t150;
t212 = -t150 * pkin(5) + t151 * pkin(10) - t104;
t210 = t155 * rSges(4,3) + t156 * t233;
t209 = t153 + t154;
t207 = t155 * (t153 * t94 + (t194 * t156 + (t193 - t93) * t155) * t156) + t241;
t130 = t157 * rSges(5,1) + t158 * rSges(5,2);
t206 = -t130 - t239;
t124 = t150 * rSges(6,1) + t151 * rSges(6,2);
t205 = -t124 - t238;
t27 = t45 + t235;
t197 = -t114 * rSges(7,1) - t113 * rSges(7,2);
t64 = rSges(7,3) * t220 - t197;
t26 = t155 * t64 + t153 * (pkin(10) * t150 + t237) + t244 * t156;
t101 = -Icges(7,5) * t151 + (Icges(7,1) * t168 - Icges(7,4) * t165) * t150;
t99 = -Icges(7,3) * t151 + (Icges(7,5) * t168 - Icges(7,6) * t165) * t150;
t32 = t113 * t100 + t114 * t101 + t99 * t220;
t3 = -t32 * t151 + (t155 * t17 + t156 * t18) * t150;
t33 = t115 * t100 + t116 * t101 + t99 * t219;
t4 = -t33 * t151 + (t155 * t19 + t156 * t20) * t150;
t204 = t3 * t242 + t4 * t243 - t151 * (t227 - t228) / 0.2e1 + t8 * t220 / 0.2e1 + t9 * t219 / 0.2e1;
t203 = t212 - t238;
t202 = -t238 - t239;
t170 = cos(qJ(1));
t161 = t170 * pkin(1);
t163 = -pkin(9) + t171;
t201 = -t155 * t163 + t119 + t161;
t200 = -t231 + t233;
t12 = t26 + t235;
t192 = Icges(4,1) * t169 - t226;
t189 = -Icges(4,2) * t166 + t225;
t186 = Icges(4,5) * t169 - Icges(4,6) * t166;
t122 = Icges(6,2) * t151 + t222;
t123 = Icges(6,1) * t150 + t221;
t181 = -t122 * t150 + t123 * t151;
t128 = Icges(5,2) * t158 + t224;
t129 = Icges(5,1) * t157 + t223;
t180 = -t128 * t157 + t129 * t158;
t179 = t240 * t156 + t241;
t176 = -t124 + t202;
t175 = t202 + t212;
t121 = Icges(6,5) * t150 + Icges(6,6) * t151;
t174 = -t228 / 0.2e1 + t227 / 0.2e1 + (t155 * t121 + t150 * t90 + t151 * t88 + t181 * t156 + t33) * t243 + (-t156 * t121 + t150 * t89 + t151 * t87 + t181 * t155 + t32) * t242;
t173 = t245 * t156 + t207;
t127 = Icges(5,5) * t157 + Icges(5,6) * t158;
t172 = t174 + (t155 * t127 + t180 * t156 + t157 * t98 + t158 * t96) * t243 + (-t156 * t127 + t180 * t155 + t157 * t97 + t158 * t95) * t242;
t144 = t170 * rSges(2,1) - t167 * rSges(2,2);
t143 = -t167 * rSges(2,1) - t170 * rSges(2,2);
t142 = t166 * rSges(4,1) + t169 * rSges(4,2);
t118 = t156 * rSges(3,1) - t155 * rSges(3,2) + t161;
t117 = -t155 * rSges(3,1) - t156 * rSges(3,2) - t236;
t106 = Icges(4,3) * t155 + t186 * t156;
t105 = -Icges(4,3) * t156 + t186 * t155;
t103 = t206 * t156;
t102 = t206 * t155;
t84 = t205 * t156;
t83 = t205 * t155;
t78 = t150 * t168 * t101;
t75 = t176 * t156;
t74 = t176 * t155;
t73 = t246 + t161 + (pkin(2) - t231) * t156 + t210;
t72 = t229 - t236 + t149 + (-pkin(2) - t200) * t155;
t69 = -t155 * t171 + t136 + t161 + t178;
t68 = -t236 + (rSges(5,3) - t171) * t156 + (-t152 - t199) * t155;
t67 = t212 * t156;
t66 = t212 * t155;
t63 = t177 + t201;
t62 = -t236 + (rSges(6,3) - t163) * t156 + (-t135 - t198) * t155;
t55 = t156 * (-t156 * t231 + t210) + (t200 * t155 - t229) * t155;
t54 = t203 * t156;
t53 = t203 * t155;
t47 = t175 * t156;
t46 = t175 * t155;
t38 = -t104 * t219 - t151 * t65;
t37 = t104 * t220 + t151 * t64;
t36 = -t150 * t213 - t151 * t99 + t78;
t35 = t201 + t244;
t34 = -t236 - t156 * t163 + (-t237 - t135 + (-rSges(7,3) - pkin(10)) * t150) * t155 + t197;
t31 = (-t155 * t65 + t156 * t64) * t150;
t30 = t50 + t234;
t15 = t27 + t234;
t11 = t12 + t234;
t1 = [t158 * t128 + t157 * t129 + t169 * (Icges(4,2) * t169 + t226) + t166 * (Icges(4,1) * t166 + t225) + Icges(2,3) + Icges(3,3) + t78 + (-t99 + t122) * t151 + (t123 - t213) * t150 + m(7) * (t34 ^ 2 + t35 ^ 2) + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t68 ^ 2 + t69 ^ 2) + m(4) * (t72 ^ 2 + t73 ^ 2) + m(3) * (t117 ^ 2 + t118 ^ 2) + m(2) * (t143 ^ 2 + t144 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(7) * (t47 * t34 + t46 * t35) + m(6) * (t75 * t62 + t74 * t63) + m(5) * (t102 * t69 + t103 * t68) + m(4) * (-t155 * t73 - t156 * t72) * t142 + (t154 / 0.2e1 + t153 / 0.2e1) * (Icges(4,5) * t166 + Icges(4,6) * t169) + t172 + (t169 * (-Icges(4,6) * t156 + t189 * t155) + t166 * (-Icges(4,5) * t156 + t192 * t155)) * t242 + (t169 * (Icges(4,6) * t155 + t189 * t156) + t166 * (Icges(4,5) * t155 + t192 * t156)) * t243; m(4) * t55 + m(5) * t30 + m(6) * t15 + m(7) * t11; m(7) * (t11 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t15 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(5) * (t102 ^ 2 + t103 ^ 2 + t30 ^ 2) + m(4) * (t209 * t142 ^ 2 + t55 ^ 2) + t155 * t153 * t106 + t207 + (-t154 * t105 + (-t155 * t105 + t156 * t106) * t155 + t245) * t156; m(7) * (t54 * t34 + t53 * t35) + m(6) * (t84 * t62 + t83 * t63) + m(5) * (-t155 * t69 - t156 * t68) * t130 + t172; m(5) * t50 + m(6) * t27 + m(7) * t12; m(7) * (t12 * t11 + t53 * t46 + t54 * t47) + m(6) * (t27 * t15 + t83 * t74 + t84 * t75) + m(5) * (t50 * t30 + (-t102 * t155 - t103 * t156) * t130) + t173; m(7) * (t12 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t27 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(5) * (t209 * t130 ^ 2 + t50 ^ 2) + t173; m(7) * (t67 * t34 + t66 * t35) + m(6) * (-t155 * t63 - t156 * t62) * t124 + t174; m(6) * t45 + m(7) * t26; m(7) * (t26 * t11 + t66 * t46 + t67 * t47) + m(6) * (t45 * t15 + (-t155 * t74 - t156 * t75) * t124) + t179; m(7) * (t26 * t12 + t66 * t53 + t67 * t54) + m(6) * (t45 * t27 + (-t155 * t83 - t156 * t84) * t124) + t179; m(6) * (t209 * t124 ^ 2 + t45 ^ 2) + m(7) * (t26 ^ 2 + t66 ^ 2 + t67 ^ 2) + t179; -t36 * t151 + m(7) * (t37 * t34 + t38 * t35) + ((t25 / 0.2e1 + t33 / 0.2e1) * t156 + (t24 / 0.2e1 + t32 / 0.2e1) * t155) * t150; m(7) * t31; m(7) * (t31 * t11 + t37 * t47 + t38 * t46) + t204; m(7) * (t31 * t12 + t37 * t54 + t38 * t53) + t204; m(7) * (t31 * t26 + t37 * t67 + t38 * t66) + t204; t151 ^ 2 * t36 + m(7) * (t31 ^ 2 + t37 ^ 2 + t38 ^ 2) + (t156 * t4 + t155 * t3 - t151 * (t155 * t24 + t156 * t25)) * t150;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
