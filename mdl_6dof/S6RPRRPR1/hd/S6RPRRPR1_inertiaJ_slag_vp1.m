% Calculate joint inertia matrix for
% S6RPRRPR1
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:57:37
% DurationCPUTime: 2.42s
% Computational Cost: add. (7772->330), mult. (5175->468), div. (0->0), fcn. (5295->12), ass. (0->170)
t246 = Icges(5,3) + Icges(6,3);
t160 = qJ(3) + qJ(4);
t153 = pkin(11) + t160;
t146 = sin(t153);
t147 = cos(t153);
t154 = sin(t160);
t155 = cos(t160);
t245 = Icges(5,5) * t155 + Icges(6,5) * t147 - Icges(5,6) * t154 - Icges(6,6) * t146;
t159 = qJ(1) + pkin(10);
t151 = sin(t159);
t152 = cos(t159);
t244 = -t151 * t245 + t152 * t246;
t243 = t151 * t246 + t152 * t245;
t213 = Icges(6,4) * t147;
t181 = -Icges(6,2) * t146 + t213;
t83 = -Icges(6,6) * t152 + t151 * t181;
t214 = Icges(6,4) * t146;
t184 = Icges(6,1) * t147 - t214;
t85 = -Icges(6,5) * t152 + t151 * t184;
t215 = Icges(5,4) * t155;
t182 = -Icges(5,2) * t154 + t215;
t91 = -Icges(5,6) * t152 + t151 * t182;
t216 = Icges(5,4) * t154;
t185 = Icges(5,1) * t155 - t216;
t93 = -Icges(5,5) * t152 + t151 * t185;
t242 = t146 * t83 - t147 * t85 + t154 * t91 - t155 * t93;
t84 = Icges(6,6) * t151 + t152 * t181;
t86 = Icges(6,5) * t151 + t152 * t184;
t92 = Icges(5,6) * t151 + t152 * t182;
t94 = Icges(5,5) * t151 + t152 * t185;
t241 = t146 * t84 - t147 * t86 + t154 * t92 - t155 * t94;
t149 = t151 ^ 2;
t240 = t151 * pkin(7);
t210 = t147 * t152;
t211 = t146 * t152;
t161 = sin(qJ(6));
t207 = t152 * t161;
t164 = cos(qJ(6));
t208 = t151 * t164;
t111 = -t147 * t207 + t208;
t206 = t152 * t164;
t209 = t151 * t161;
t112 = t147 * t206 + t209;
t63 = rSges(7,1) * t112 + rSges(7,2) * t111 + rSges(7,3) * t211;
t239 = pkin(5) * t210 + pkin(9) * t211 + t63;
t238 = Icges(5,5) * t154 + Icges(6,5) * t146 + Icges(5,6) * t155 + Icges(6,6) * t147;
t117 = Icges(6,2) * t147 + t214;
t118 = Icges(6,1) * t146 + t213;
t123 = Icges(5,2) * t155 + t216;
t124 = Icges(5,1) * t154 + t215;
t237 = -t117 * t146 + t118 * t147 - t123 * t154 + t124 * t155;
t193 = rSges(5,1) * t155 - rSges(5,2) * t154;
t98 = -t147 * rSges(7,3) + (rSges(7,1) * t164 - rSges(7,2) * t161) * t146;
t236 = -pkin(5) * t146 + pkin(9) * t147 - t98;
t150 = t152 ^ 2;
t109 = -t147 * t209 - t206;
t110 = t147 * t208 - t207;
t212 = t146 * t151;
t53 = Icges(7,5) * t110 + Icges(7,6) * t109 + Icges(7,3) * t212;
t55 = Icges(7,4) * t110 + Icges(7,2) * t109 + Icges(7,6) * t212;
t57 = Icges(7,1) * t110 + Icges(7,4) * t109 + Icges(7,5) * t212;
t16 = t109 * t55 + t110 * t57 + t212 * t53;
t54 = Icges(7,5) * t112 + Icges(7,6) * t111 + Icges(7,3) * t211;
t56 = Icges(7,4) * t112 + Icges(7,2) * t111 + Icges(7,6) * t211;
t58 = Icges(7,1) * t112 + Icges(7,4) * t111 + Icges(7,5) * t211;
t17 = t109 * t56 + t110 * t58 + t212 * t54;
t8 = t151 * t17 - t152 * t16;
t235 = -t8 + t244 * t150 + (t241 * t151 + (-t242 + t243) * t152) * t151;
t167 = -pkin(8) - pkin(7);
t234 = t151 / 0.2e1;
t233 = -t152 / 0.2e1;
t162 = sin(qJ(3));
t232 = pkin(3) * t162;
t231 = pkin(4) * t154;
t230 = pkin(5) * t147;
t163 = sin(qJ(1));
t229 = t163 * pkin(1);
t165 = cos(qJ(3));
t148 = pkin(3) * t165 + pkin(2);
t130 = pkin(4) * t155 + t148;
t115 = t152 * t130;
t132 = t152 * t148;
t228 = t152 * (t115 - t132) + (t130 - t148) * t149;
t145 = t152 * pkin(7);
t227 = t151 * (t145 + (-pkin(2) + t148) * t151) + t152 * (-t152 * pkin(2) + t132 - t240);
t173 = t151 * rSges(5,3) + t152 * t193;
t48 = t151 * (-t152 * rSges(5,3) + t151 * t193) + t152 * t173;
t226 = rSges(4,1) * t165;
t224 = rSges(4,2) * t162;
t222 = t152 * rSges(4,3);
t96 = -Icges(7,6) * t147 + (Icges(7,4) * t164 - Icges(7,2) * t161) * t146;
t221 = t161 * t96;
t24 = -t147 * t53 + (-t161 * t55 + t164 * t57) * t146;
t220 = t24 * t152;
t25 = -t147 * t54 + (-t161 * t56 + t164 * t58) * t146;
t219 = t25 * t151;
t218 = Icges(4,4) * t162;
t217 = Icges(4,4) * t165;
t204 = rSges(4,3) * t151 + t152 * t226;
t203 = t149 + t150;
t18 = t111 * t55 + t112 * t57 + t211 * t53;
t19 = t111 * t56 + t112 * t58 + t211 * t54;
t9 = t151 * t19 - t152 * t18;
t201 = (t9 + t243 * t149 + ((-t241 + t244) * t151 + t242 * t152) * t152) * t151;
t125 = rSges(5,1) * t154 + rSges(5,2) * t155;
t200 = -t125 - t232;
t119 = rSges(6,1) * t146 + rSges(6,2) * t147;
t199 = -t119 - t231;
t172 = rSges(6,1) * t210 - rSges(6,2) * t211 + rSges(6,3) * t151;
t192 = rSges(6,1) * t147 - rSges(6,2) * t146;
t26 = t151 * (-t152 * rSges(6,3) + t151 * t192) + t152 * t172 + t228;
t95 = -Icges(7,3) * t147 + (Icges(7,5) * t164 - Icges(7,6) * t161) * t146;
t97 = -Icges(7,5) * t147 + (Icges(7,1) * t164 - Icges(7,4) * t161) * t146;
t30 = t109 * t96 + t110 * t97 + t212 * t95;
t3 = -t30 * t147 + (t151 * t16 + t152 * t17) * t146;
t31 = t111 * t96 + t112 * t97 + t211 * t95;
t4 = -t31 * t147 + (t151 * t18 + t152 * t19) * t146;
t198 = t3 * t233 + t4 * t234 - t147 * (t219 - t220) / 0.2e1 + t8 * t212 / 0.2e1 + t9 * t211 / 0.2e1;
t197 = -t231 + t236;
t196 = -t231 - t232;
t166 = cos(qJ(1));
t157 = t166 * pkin(1);
t158 = -qJ(5) + t167;
t195 = -t151 * t158 + t115 + t157;
t194 = -t224 + t226;
t191 = -t110 * rSges(7,1) - t109 * rSges(7,2);
t62 = rSges(7,3) * t212 - t191;
t12 = t151 * t62 + t228 + t149 * (pkin(9) * t146 + t230) + t239 * t152;
t186 = Icges(4,1) * t165 - t218;
t183 = -Icges(4,2) * t162 + t217;
t180 = Icges(4,5) * t165 - Icges(4,6) * t162;
t171 = -t119 + t196;
t170 = t196 + t236;
t169 = t152 * t235 + t201;
t168 = t219 / 0.2e1 - t220 / 0.2e1 + (t146 * t86 + t147 * t84 + t151 * t238 + t152 * t237 + t154 * t94 + t155 * t92 + t31) * t234 + (t146 * t85 + t147 * t83 + t151 * t237 - t152 * t238 + t154 * t93 + t155 * t91 + t30) * t233;
t140 = rSges(2,1) * t166 - rSges(2,2) * t163;
t139 = -rSges(2,1) * t163 - rSges(2,2) * t166;
t138 = rSges(4,1) * t162 + rSges(4,2) * t165;
t114 = rSges(3,1) * t152 - rSges(3,2) * t151 + t157;
t113 = -rSges(3,1) * t151 - rSges(3,2) * t152 - t229;
t102 = Icges(4,3) * t151 + t152 * t180;
t101 = -Icges(4,3) * t152 + t151 * t180;
t100 = t200 * t152;
t99 = t200 * t151;
t80 = t199 * t152;
t79 = t199 * t151;
t72 = t146 * t164 * t97;
t71 = t171 * t152;
t70 = t171 * t151;
t69 = t240 + t157 + (pkin(2) - t224) * t152 + t204;
t68 = t222 - t229 + t145 + (-pkin(2) - t194) * t151;
t65 = -t151 * t167 + t132 + t157 + t173;
t64 = -t229 + (rSges(5,3) - t167) * t152 + (-t148 - t193) * t151;
t61 = t152 * (-t152 * t224 + t204) + (t151 * t194 - t222) * t151;
t60 = t172 + t195;
t59 = -t229 + (rSges(6,3) - t158) * t152 + (-t130 - t192) * t151;
t52 = t197 * t152;
t51 = t197 * t151;
t45 = t170 * t152;
t44 = t170 * t151;
t37 = -t147 * t63 - t211 * t98;
t36 = t147 * t62 + t212 * t98;
t35 = -t146 * t221 - t147 * t95 + t72;
t34 = t195 + t239;
t33 = -t229 - t152 * t158 + (-t230 - t130 + (-rSges(7,3) - pkin(9)) * t146) * t151 + t191;
t32 = t48 + t227;
t29 = (-t151 * t63 + t152 * t62) * t146;
t15 = t26 + t227;
t11 = t12 + t227;
t1 = [t155 * t123 + t154 * t124 + t165 * (Icges(4,2) * t165 + t218) + t162 * (Icges(4,1) * t162 + t217) + Icges(2,3) + Icges(3,3) + t72 + (-t95 + t117) * t147 + (t118 - t221) * t146 + m(7) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t59 ^ 2 + t60 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) + m(4) * (t68 ^ 2 + t69 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t139 ^ 2 + t140 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); (t150 / 0.2e1 + t149 / 0.2e1) * (Icges(4,5) * t162 + Icges(4,6) * t165) + m(7) * (t33 * t45 + t34 * t44) + m(6) * (t59 * t71 + t60 * t70) + m(5) * (t100 * t64 + t65 * t99) + m(4) * (-t151 * t69 - t152 * t68) * t138 + (t165 * (Icges(4,6) * t151 + t152 * t183) + t162 * (Icges(4,5) * t151 + t152 * t186)) * t234 + (t165 * (-Icges(4,6) * t152 + t151 * t183) + t162 * (-Icges(4,5) * t152 + t151 * t186)) * t233 + t168; m(4) * t61 + m(5) * t32 + m(6) * t15 + m(7) * t11; m(7) * (t11 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t15 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t100 ^ 2 + t32 ^ 2 + t99 ^ 2) + t151 * t149 * t102 + m(4) * (t138 ^ 2 * t203 + t61 ^ 2) + t201 + (-t150 * t101 + (-t151 * t101 + t152 * t102) * t151 + t235) * t152; m(5) * (-t151 * t65 - t152 * t64) * t125 + m(7) * (t33 * t52 + t34 * t51) + m(6) * (t59 * t80 + t60 * t79) + t168; m(5) * t48 + m(6) * t26 + m(7) * t12; m(7) * (t11 * t12 + t44 * t51 + t45 * t52) + m(6) * (t15 * t26 + t70 * t79 + t71 * t80) + m(5) * (t48 * t32 + (-t100 * t152 - t151 * t99) * t125) + t169; m(7) * (t12 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t26 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t125 ^ 2 * t203 + t48 ^ 2) + t169; m(7) * (t151 * t33 - t152 * t34) + m(6) * (t151 * t59 - t152 * t60); 0; m(7) * (t151 * t45 - t152 * t44) + m(6) * (t151 * t71 - t152 * t70); m(7) * (t151 * t52 - t152 * t51) + m(6) * (t151 * t80 - t152 * t79); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t203; m(7) * (t33 * t36 + t34 * t37) - t35 * t147 + ((t31 / 0.2e1 + t25 / 0.2e1) * t152 + (t30 / 0.2e1 + t24 / 0.2e1) * t151) * t146; m(7) * t29; m(7) * (t11 * t29 + t36 * t45 + t37 * t44) + t198; m(7) * (t12 * t29 + t36 * t52 + t37 * t51) + t198; m(7) * (t151 * t36 - t152 * t37); t147 ^ 2 * t35 + m(7) * (t29 ^ 2 + t36 ^ 2 + t37 ^ 2) + (t152 * t4 + t151 * t3 - t147 * (t151 * t24 + t152 * t25)) * t146;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
