% Calculate joint inertia matrix for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:42
% EndTime: 2019-03-08 19:46:52
% DurationCPUTime: 4.65s
% Computational Cost: add. (12298->476), mult. (27095->697), div. (0->0), fcn. (34211->12), ass. (0->217)
t258 = Icges(3,1) + Icges(4,2);
t256 = Icges(3,4) + Icges(4,6);
t255 = Icges(3,5) - Icges(4,4);
t254 = Icges(4,5) - Icges(3,6);
t257 = Icges(3,2) + Icges(4,3);
t253 = Icges(3,3) + Icges(4,1);
t198 = cos(pkin(6));
t194 = sin(pkin(10));
t202 = cos(qJ(2));
t222 = t202 * t194;
t197 = cos(pkin(10));
t201 = sin(qJ(2));
t224 = t197 * t201;
t181 = t198 * t222 + t224;
t221 = t202 * t197;
t223 = t201 * t194;
t182 = -t198 * t223 + t221;
t195 = sin(pkin(6));
t229 = t194 * t195;
t252 = t257 * t181 - t256 * t182 + t254 * t229;
t179 = -t198 * t221 + t223;
t180 = t198 * t224 + t222;
t228 = t195 * t197;
t251 = t257 * t179 - t256 * t180 - t254 * t228;
t250 = t256 * t181 - t258 * t182 - t255 * t229;
t249 = t256 * t179 - t258 * t180 + t255 * t228;
t215 = m(6) / 0.2e1 + m(7) / 0.2e1;
t248 = 0.2e1 * t215;
t247 = t254 * t181 + t255 * t182 + t253 * t229;
t246 = -t254 * t179 - t255 * t180 + t253 * t228;
t245 = t253 * t198 + (t255 * t201 - t254 * t202) * t195;
t244 = t255 * t198 + (t258 * t201 + t256 * t202) * t195;
t243 = t254 * t198 + (-t256 * t201 - t257 * t202) * t195;
t200 = sin(qJ(4));
t227 = t195 * t200;
t236 = cos(qJ(4));
t161 = -t181 * t236 + t194 * t227;
t163 = t179 * t236 + t197 * t227;
t209 = t195 * t236;
t183 = t198 * t200 + t202 * t209;
t162 = t181 * t200 + t194 * t209;
t192 = pkin(11) + qJ(6);
t190 = sin(t192);
t191 = cos(t192);
t120 = -t162 * t190 + t182 * t191;
t121 = t162 * t191 + t182 * t190;
t71 = Icges(7,5) * t121 + Icges(7,6) * t120 + Icges(7,3) * t161;
t73 = Icges(7,4) * t121 + Icges(7,2) * t120 + Icges(7,6) * t161;
t75 = Icges(7,1) * t121 + Icges(7,4) * t120 + Icges(7,5) * t161;
t28 = t120 * t73 + t121 * t75 + t161 * t71;
t165 = t179 * t200 - t197 * t209;
t122 = -t165 * t190 + t180 * t191;
t123 = t165 * t191 + t180 * t190;
t72 = Icges(7,5) * t123 + Icges(7,6) * t122 - Icges(7,3) * t163;
t74 = Icges(7,4) * t123 + Icges(7,2) * t122 - Icges(7,6) * t163;
t76 = Icges(7,1) * t123 + Icges(7,4) * t122 - Icges(7,5) * t163;
t29 = t120 * t74 + t121 * t76 + t161 * t72;
t225 = t195 * t202;
t184 = t198 * t236 - t200 * t225;
t226 = t195 * t201;
t156 = -t184 * t190 + t191 * t226;
t157 = t184 * t191 + t190 * t226;
t93 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t183;
t94 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t183;
t95 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t183;
t44 = t120 * t94 + t121 * t95 + t161 * t93;
t1 = t161 * t28 - t163 * t29 + t183 * t44;
t241 = t1 / 0.2e1;
t37 = t156 * t73 + t157 * t75 + t183 * t71;
t38 = t156 * t74 + t157 * t76 + t183 * t72;
t51 = t156 * t94 + t157 * t95 + t183 * t93;
t11 = t161 * t37 - t163 * t38 + t183 * t51;
t240 = t11 / 0.2e1;
t239 = t161 / 0.2e1;
t238 = -t163 / 0.2e1;
t237 = t183 / 0.2e1;
t196 = cos(pkin(11));
t235 = pkin(5) * t196;
t193 = sin(pkin(11));
t230 = t182 * t193;
t77 = rSges(7,1) * t121 + rSges(7,2) * t120 + rSges(7,3) * t161;
t234 = pkin(5) * t230 + pkin(9) * t161 + t162 * t235 + t77;
t231 = t180 * t193;
t78 = rSges(7,1) * t123 + rSges(7,2) * t122 - rSges(7,3) * t163;
t233 = pkin(5) * t231 - pkin(9) * t163 + t165 * t235 + t78;
t213 = t193 * t226;
t96 = rSges(7,1) * t157 + rSges(7,2) * t156 + rSges(7,3) * t183;
t232 = pkin(5) * t213 + pkin(9) * t183 + t184 * t235 + t96;
t153 = pkin(2) * t180 + qJ(3) * t179;
t154 = pkin(2) * t182 + qJ(3) * t181;
t219 = t153 * t229 + t154 * t228;
t152 = t198 * t154;
t168 = pkin(3) * t229 + pkin(8) * t182;
t218 = t198 * t168 + t152;
t169 = -pkin(3) * t228 + pkin(8) * t180;
t217 = -t153 - t169;
t185 = (pkin(2) * t201 - qJ(3) * t202) * t195;
t216 = -pkin(3) * t198 - pkin(8) * t226 - t185;
t214 = -m(4) - m(5) - m(6) - m(7);
t118 = pkin(4) * t162 + qJ(5) * t161;
t212 = t198 * t118 + t218;
t119 = pkin(4) * t165 - qJ(5) * t163;
t211 = -t119 + t217;
t155 = pkin(4) * t184 + qJ(5) * t183;
t210 = -t155 + t216;
t208 = (-t198 * rSges(4,1) - (-rSges(4,2) * t201 - rSges(4,3) * t202) * t195 - t185) * t195;
t207 = t168 * t228 + t169 * t229 + t219;
t148 = rSges(5,1) * t184 - rSges(5,2) * t183 + rSges(5,3) * t226;
t206 = (-t148 + t216) * t195;
t166 = -t184 * t193 + t196 * t226;
t167 = t184 * t196 + t213;
t109 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t183;
t205 = (-t109 + t210) * t195;
t204 = t118 * t228 + t119 * t229 + t207;
t203 = (t210 - t232) * t195;
t176 = t198 * rSges(3,3) + (rSges(3,1) * t201 + rSges(3,2) * t202) * t195;
t147 = Icges(5,1) * t184 - Icges(5,4) * t183 + Icges(5,5) * t226;
t146 = Icges(5,4) * t184 - Icges(5,2) * t183 + Icges(5,6) * t226;
t145 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t226;
t144 = rSges(3,1) * t182 - rSges(3,2) * t181 + rSges(3,3) * t229;
t143 = rSges(3,1) * t180 - rSges(3,2) * t179 - rSges(3,3) * t228;
t142 = -rSges(4,1) * t228 - rSges(4,2) * t180 + rSges(4,3) * t179;
t141 = rSges(4,1) * t229 - rSges(4,2) * t182 + rSges(4,3) * t181;
t128 = t180 * t155;
t127 = t165 * t196 + t231;
t126 = -t165 * t193 + t180 * t196;
t125 = t162 * t196 + t230;
t124 = -t162 * t193 + t182 * t196;
t116 = t118 * t226;
t113 = -t143 * t198 - t176 * t228;
t112 = t144 * t198 - t176 * t229;
t110 = t182 * t119;
t108 = rSges(5,1) * t165 + rSges(5,2) * t163 + rSges(5,3) * t180;
t107 = rSges(5,1) * t162 - rSges(5,2) * t161 + rSges(5,3) * t182;
t106 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t183;
t105 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t183;
t104 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t183;
t103 = Icges(5,1) * t165 + Icges(5,4) * t163 + Icges(5,5) * t180;
t102 = Icges(5,1) * t162 - Icges(5,4) * t161 + Icges(5,5) * t182;
t101 = Icges(5,4) * t165 + Icges(5,2) * t163 + Icges(5,6) * t180;
t100 = Icges(5,4) * t162 - Icges(5,2) * t161 + Icges(5,6) * t182;
t99 = Icges(5,5) * t165 + Icges(5,6) * t163 + Icges(5,3) * t180;
t98 = Icges(5,5) * t162 - Icges(5,6) * t161 + Icges(5,3) * t182;
t92 = (t143 * t194 + t144 * t197) * t195;
t90 = (-t142 - t153) * t198 + t197 * t208;
t89 = t141 * t198 + t194 * t208 + t152;
t88 = rSges(6,1) * t127 + rSges(6,2) * t126 - rSges(6,3) * t163;
t87 = rSges(6,1) * t125 + rSges(6,2) * t124 + rSges(6,3) * t161;
t86 = Icges(6,1) * t127 + Icges(6,4) * t126 - Icges(6,5) * t163;
t85 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t161;
t84 = Icges(6,4) * t127 + Icges(6,2) * t126 - Icges(6,6) * t163;
t83 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t161;
t82 = Icges(6,5) * t127 + Icges(6,6) * t126 - Icges(6,3) * t163;
t81 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t161;
t80 = t107 * t226 - t148 * t182;
t79 = -t108 * t226 + t148 * t180;
t68 = t145 * t226 - t146 * t183 + t147 * t184;
t67 = (t141 * t197 + t142 * t194) * t195 + t219;
t66 = -t107 * t180 + t108 * t182;
t65 = t145 * t180 + t146 * t163 + t147 * t165;
t64 = t145 * t182 - t146 * t161 + t147 * t162;
t63 = (-t108 + t217) * t198 + t197 * t206;
t62 = t107 * t198 + t194 * t206 + t218;
t61 = -t101 * t183 + t103 * t184 + t226 * t99;
t60 = -t100 * t183 + t102 * t184 + t226 * t98;
t59 = -t163 * t96 - t183 * t78;
t58 = -t161 * t96 + t183 * t77;
t57 = t104 * t183 + t105 * t166 + t106 * t167;
t56 = t101 * t163 + t103 * t165 + t180 * t99;
t55 = t100 * t163 + t102 * t165 + t180 * t98;
t54 = -t101 * t161 + t103 * t162 + t182 * t99;
t53 = -t100 * t161 + t102 * t162 + t182 * t98;
t52 = (t107 * t197 + t108 * t194) * t195 + t207;
t50 = t161 * t78 + t163 * t77;
t49 = t87 * t226 + t116 + (-t109 - t155) * t182;
t48 = t109 * t180 + t128 + (-t119 - t88) * t226;
t47 = -t104 * t163 + t105 * t126 + t106 * t127;
t46 = t104 * t161 + t105 * t124 + t106 * t125;
t45 = t122 * t94 + t123 * t95 - t163 * t93;
t43 = (-t88 + t211) * t198 + t197 * t205;
t42 = t194 * t205 + t198 * t87 + t212;
t41 = t182 * t88 + t110 + (-t118 - t87) * t180;
t40 = t166 * t84 + t167 * t86 + t183 * t82;
t39 = t166 * t83 + t167 * t85 + t183 * t81;
t36 = t126 * t84 + t127 * t86 - t163 * t82;
t35 = t126 * t83 + t127 * t85 - t163 * t81;
t34 = t124 * t84 + t125 * t86 + t161 * t82;
t33 = t124 * t83 + t125 * t85 + t161 * t81;
t32 = (t194 * t88 + t197 * t87) * t195 + t204;
t31 = t122 * t74 + t123 * t76 - t163 * t72;
t30 = t122 * t73 + t123 * t75 - t163 * t71;
t27 = t116 + t234 * t226 + (-t155 - t232) * t182;
t26 = t128 + t232 * t180 + (-t119 - t233) * t226;
t25 = (t211 - t233) * t198 + t197 * t203;
t24 = t194 * t203 + t198 * t234 + t212;
t23 = t198 * t68 + (t194 * t60 - t197 * t61) * t195;
t22 = t180 * t61 + t182 * t60 + t226 * t68;
t21 = t110 + t233 * t182 + (-t118 - t234) * t180;
t20 = t198 * t65 + (t194 * t55 - t197 * t56) * t195;
t19 = t198 * t64 + (t194 * t53 - t197 * t54) * t195;
t18 = t180 * t56 + t182 * t55 + t226 * t65;
t17 = t180 * t54 + t182 * t53 + t226 * t64;
t16 = (t194 * t233 + t197 * t234) * t195 + t204;
t15 = t198 * t57 + (t194 * t39 - t197 * t40) * t195;
t14 = t180 * t40 + t182 * t39 + t226 * t57;
t13 = t198 * t51 + (t194 * t37 - t197 * t38) * t195;
t12 = t180 * t38 + t182 * t37 + t226 * t51;
t10 = t198 * t47 + (t194 * t35 - t197 * t36) * t195;
t9 = t198 * t46 + (t194 * t33 - t197 * t34) * t195;
t8 = t180 * t36 + t182 * t35 + t226 * t47;
t7 = t180 * t34 + t182 * t33 + t226 * t46;
t6 = t198 * t45 + (t194 * t30 - t197 * t31) * t195;
t5 = t198 * t44 + (t194 * t28 - t197 * t29) * t195;
t4 = t180 * t31 + t182 * t30 + t226 * t45;
t3 = t180 * t29 + t182 * t28 + t226 * t44;
t2 = t161 * t30 - t163 * t31 + t183 * t45;
t69 = [m(2) + m(3) - t214; m(3) * t92 + m(4) * t67 + m(5) * t52 + m(6) * t32 + m(7) * t16; m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t52 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t67 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(3) * (t112 ^ 2 + t113 ^ 2 + t92 ^ 2) + (t13 + t15 + t23 + t245 * t198 ^ 2 + ((t201 * t244 - t202 * t243) * t198 + (t246 * t198 + (t249 * t201 + t251 * t202) * t195) * t197 + (t247 * t198 + (-t250 * t201 - t252 * t202) * t195) * t194) * t195) * t198 + (t5 + t9 + t19 + (t252 * t181 - t250 * t182 + t247 * t229) * t229 + (t181 * t243 + t182 * t244 + t229 * t245) * t198) * t229 + (-t6 - t10 - t20 + (t251 * t179 - t249 * t180 + t246 * t228) * t228 + (-t179 * t243 - t180 * t244 + t228 * t245) * t198 + (-t252 * t179 + t250 * t180 - t251 * t181 + t249 * t182 + t247 * t228 + t246 * t229) * t229) * t228; t214 * t225; m(7) * (-t16 * t225 + t179 * t24 + t181 * t25) + m(6) * (t179 * t42 + t181 * t43 - t225 * t32) + m(5) * (t179 * t62 + t181 * t63 - t225 * t52) + m(4) * (t179 * t89 + t181 * t90 - t225 * t67); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t215) * (t195 ^ 2 * t202 ^ 2 + t179 ^ 2 + t181 ^ 2); m(5) * t66 + m(6) * t41 + m(7) * t21; (t12 / 0.2e1 + t14 / 0.2e1 + t22 / 0.2e1) * t198 + (t5 / 0.2e1 + t9 / 0.2e1 + t19 / 0.2e1) * t182 + (t6 / 0.2e1 + t10 / 0.2e1 + t20 / 0.2e1) * t180 + m(7) * (t16 * t21 + t24 * t27 + t25 * t26) + m(6) * (t32 * t41 + t42 * t49 + t43 * t48) + m(5) * (t52 * t66 + t62 * t80 + t63 * t79) + ((t13 / 0.2e1 + t23 / 0.2e1 + t15 / 0.2e1) * t201 + (-t4 / 0.2e1 - t18 / 0.2e1 - t8 / 0.2e1) * t197 + (t3 / 0.2e1 + t17 / 0.2e1 + t7 / 0.2e1) * t194) * t195; m(5) * (t80 * t179 + t79 * t181 - t225 * t66) + m(6) * (t49 * t179 + t48 * t181 - t225 * t41) + m(7) * (t27 * t179 + t26 * t181 - t21 * t225); (t12 + t14 + t22) * t226 + (t3 + t7 + t17) * t182 + (t4 + t18 + t8) * t180 + m(7) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t66 ^ 2 + t79 ^ 2 + t80 ^ 2); t183 * t248; m(7) * (t16 * t183 + t161 * t25 - t163 * t24) + m(6) * (t161 * t43 - t163 * t42 + t183 * t32); (t161 * t181 - t163 * t179 - t183 * t225) * t248; m(7) * (t161 * t26 - t163 * t27 + t183 * t21) + m(6) * (t161 * t48 - t163 * t49 + t183 * t41); (t161 ^ 2 + t163 ^ 2 + t183 ^ 2) * t248; m(7) * t50; t198 * t240 + t13 * t237 + t6 * t238 + t5 * t239 + m(7) * (t16 * t50 + t24 * t58 + t25 * t59) + (-t197 * t2 / 0.2e1 + t194 * t241) * t195; m(7) * (t58 * t179 + t59 * t181 - t225 * t50); t180 * t2 / 0.2e1 + m(7) * (t21 * t50 + t26 * t59 + t27 * t58) + t12 * t237 + t182 * t241 + t3 * t239 + t4 * t238 + t226 * t240; m(7) * (t161 * t59 - t163 * t58 + t183 * t50); t161 * t1 - t163 * t2 + t183 * t11 + m(7) * (t50 ^ 2 + t58 ^ 2 + t59 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t69(1) t69(2) t69(4) t69(7) t69(11) t69(16); t69(2) t69(3) t69(5) t69(8) t69(12) t69(17); t69(4) t69(5) t69(6) t69(9) t69(13) t69(18); t69(7) t69(8) t69(9) t69(10) t69(14) t69(19); t69(11) t69(12) t69(13) t69(14) t69(15) t69(20); t69(16) t69(17) t69(18) t69(19) t69(20) t69(21);];
Mq  = res;
