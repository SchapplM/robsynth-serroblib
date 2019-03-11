% Calculate joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:28
% EndTime: 2019-03-09 03:44:34
% DurationCPUTime: 2.61s
% Computational Cost: add. (7058->372), mult. (6896->534), div. (0->0), fcn. (7340->10), ass. (0->186)
t248 = Icges(5,1) + Icges(4,3);
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t247 = (-Icges(5,4) + Icges(4,5)) * t172 + (Icges(5,5) - Icges(4,6)) * t169;
t165 = qJ(1) + pkin(10);
t160 = sin(t165);
t246 = -t160 / 0.2e1;
t239 = t160 / 0.2e1;
t161 = cos(t165);
t238 = -t161 / 0.2e1;
t245 = t161 / 0.2e1;
t237 = t169 / 0.2e1;
t244 = -t247 * t160 + t248 * t161;
t243 = t248 * t160 + t247 * t161;
t158 = t160 ^ 2;
t159 = t161 ^ 2;
t242 = m(5) / 0.2e1;
t241 = m(6) / 0.2e1;
t240 = m(7) / 0.2e1;
t140 = t169 * rSges(4,1) + t172 * rSges(4,2);
t236 = m(4) * t140;
t168 = sin(qJ(5));
t235 = pkin(5) * t168;
t170 = sin(qJ(1));
t234 = t170 * pkin(1);
t174 = -pkin(9) - pkin(8);
t233 = -pkin(8) - t174;
t167 = qJ(5) + qJ(6);
t162 = sin(t167);
t163 = cos(t167);
t218 = t163 * t169;
t104 = -t160 * t162 + t161 * t218;
t219 = t162 * t169;
t105 = t160 * t163 + t161 * t219;
t220 = t161 * t172;
t68 = t105 * rSges(7,1) + t104 * rSges(7,2) + rSges(7,3) * t220;
t171 = cos(qJ(5));
t157 = t171 * pkin(5) + pkin(4);
t217 = t168 * t169;
t205 = t161 * t217;
t175 = pkin(5) * t205 + t160 * t157 - t174 * t220;
t211 = t160 * pkin(4) + pkin(8) * t220;
t78 = t175 - t211;
t232 = t68 + t78;
t106 = t160 * t218 + t161 * t162;
t107 = t160 * t219 - t161 * t163;
t192 = -t107 * rSges(7,1) - t106 * rSges(7,2);
t223 = t160 * t172;
t69 = rSges(7,3) * t223 - t192;
t155 = t161 * pkin(4);
t222 = t161 * t157;
t79 = -t222 + t155 + (pkin(5) * t217 + t172 * t233) * t160;
t231 = t69 + t79;
t230 = t161 * rSges(5,1);
t229 = t161 * rSges(4,3);
t228 = Icges(4,4) * t169;
t227 = Icges(4,4) * t172;
t226 = Icges(5,6) * t169;
t225 = Icges(5,6) * t172;
t224 = qJ(4) * t169;
t221 = t161 * t169;
t216 = t169 * t171;
t212 = pkin(3) * t220 + qJ(4) * t221;
t215 = t158 * (pkin(3) * t172 + t224) + t161 * t212;
t111 = t169 * rSges(7,3) + (-rSges(7,1) * t162 - rSges(7,2) * t163) * t172;
t128 = t169 * t233 - t172 * t235;
t214 = -t111 - t128;
t138 = t169 * pkin(3) - t172 * qJ(4);
t213 = t169 * rSges(5,2) + t172 * rSges(5,3) - t138;
t210 = t158 + t159;
t209 = -m(5) - m(6) - m(7);
t62 = Icges(7,5) * t105 + Icges(7,6) * t104 + Icges(7,3) * t220;
t64 = Icges(7,4) * t105 + Icges(7,2) * t104 + Icges(7,6) * t220;
t66 = Icges(7,1) * t105 + Icges(7,4) * t104 + Icges(7,5) * t220;
t30 = t169 * t62 + (-t162 * t66 - t163 * t64) * t172;
t63 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t223;
t65 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t223;
t67 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t223;
t31 = t169 * t63 + (-t162 * t67 - t163 * t65) * t172;
t19 = t104 * t64 + t105 * t66 + t220 * t62;
t20 = t104 * t65 + t105 * t67 + t220 * t63;
t108 = Icges(7,3) * t169 + (-Icges(7,5) * t162 - Icges(7,6) * t163) * t172;
t109 = Icges(7,6) * t169 + (-Icges(7,4) * t162 - Icges(7,2) * t163) * t172;
t110 = Icges(7,5) * t169 + (-Icges(7,1) * t162 - Icges(7,4) * t163) * t172;
t41 = t104 * t109 + t105 * t110 + t108 * t220;
t5 = t41 * t169 + (t160 * t20 + t161 * t19) * t172;
t180 = -t163 * t109 - t162 * t110;
t87 = t169 * t108;
t51 = (t172 * t180 + t87) * t169;
t21 = t106 * t64 + t107 * t66 + t223 * t62;
t22 = t106 * t65 + t107 * t67 + t223 * t63;
t42 = t106 * t109 + t107 * t110 + t108 * t223;
t6 = t42 * t169 + (t160 * t22 + t161 * t21) * t172;
t208 = t5 * t220 + t6 * t223 + t169 * (t51 + (t160 * t31 + t161 * t30) * t172);
t118 = -t160 * t168 + t161 * t216;
t119 = t160 * t171 + t205;
t70 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t220;
t72 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t220;
t74 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t220;
t32 = t169 * t70 + (-t168 * t74 - t171 * t72) * t172;
t122 = Icges(6,3) * t169 + (-Icges(6,5) * t168 - Icges(6,6) * t171) * t172;
t123 = Icges(6,6) * t169 + (-Icges(6,4) * t168 - Icges(6,2) * t171) * t172;
t124 = Icges(6,5) * t169 + (-Icges(6,1) * t168 - Icges(6,4) * t171) * t172;
t45 = t118 * t123 + t119 * t124 + t122 * t220;
t207 = t45 / 0.2e1 + t32 / 0.2e1;
t120 = t160 * t216 + t161 * t168;
t121 = t160 * t217 - t161 * t171;
t71 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t223;
t73 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t223;
t75 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t223;
t33 = t169 * t71 + (-t168 * t75 - t171 * t73) * t172;
t46 = t120 * t123 + t121 * t124 + t122 * t223;
t206 = t46 / 0.2e1 + t33 / 0.2e1;
t76 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t220;
t173 = cos(qJ(1));
t164 = t173 * pkin(1);
t204 = t161 * pkin(2) + t160 * pkin(7) + t164;
t203 = t223 / 0.2e1;
t202 = t220 / 0.2e1;
t201 = Icges(4,5) * t237 - Icges(5,4) * t169 / 0.2e1 + (Icges(4,6) / 0.2e1 - Icges(5,5) / 0.2e1) * t172;
t200 = t161 * pkin(7) - t234;
t199 = -t169 * pkin(8) - t138;
t198 = t160 * (pkin(8) * t223 - t155) + t161 * t211 + t215;
t12 = t19 * t160 - t20 * t161;
t13 = t21 * t160 - t22 * t161;
t197 = t12 * t202 + t13 * t203 + t6 * t238 + t5 * t239 + (t30 * t160 - t31 * t161) * t237;
t196 = t51 + (t31 + t42) * t203 + (t30 + t41) * t202;
t125 = t169 * rSges(6,3) + (-rSges(6,1) * t168 - rSges(6,2) * t171) * t172;
t195 = -t125 + t199;
t194 = rSges(4,1) * t172 - rSges(4,2) * t169;
t193 = -t121 * rSges(6,1) - t120 * rSges(6,2);
t187 = t204 + t212;
t186 = Icges(4,1) * t172 - t228;
t185 = -Icges(4,2) * t169 + t227;
t182 = -Icges(5,2) * t172 + t226;
t181 = Icges(5,3) * t169 - t225;
t179 = -t171 * t123 - t168 * t124;
t178 = t199 + t214;
t177 = rSges(4,1) * t220 - rSges(4,2) * t221 + t160 * rSges(4,3);
t176 = t160 * rSges(5,1) - rSges(5,2) * t220 + rSges(5,3) * t221;
t142 = t173 * rSges(2,1) - t170 * rSges(2,2);
t141 = -t170 * rSges(2,1) - t173 * rSges(2,2);
t127 = t161 * rSges(3,1) - t160 * rSges(3,2) + t164;
t126 = -t160 * rSges(3,1) - t161 * rSges(3,2) - t234;
t112 = t169 * t122;
t86 = t213 * t161;
t85 = t213 * t160;
t84 = t111 * t223;
t83 = t177 + t204;
t82 = t229 + (-pkin(2) - t194) * t160 + t200;
t81 = t195 * t161;
t80 = t195 * t160;
t77 = rSges(6,3) * t223 - t193;
t61 = t176 + t187;
t60 = t230 + (-pkin(2) + (rSges(5,2) - pkin(3)) * t172 + (-rSges(5,3) - qJ(4)) * t169) * t160 + t200;
t59 = t161 * t177 + (t160 * t194 - t229) * t160;
t58 = t169 * t68;
t57 = t69 * t220;
t56 = t178 * t161;
t55 = t178 * t160;
t54 = (t172 * t179 + t112) * t169;
t53 = -t125 * t220 + t169 * t76;
t52 = t125 * t223 - t169 * t77;
t50 = -t111 * t220 + t58;
t49 = -t169 * t69 + t84;
t48 = t187 + t76 + t211;
t47 = t155 + (-t224 - pkin(2) + (-rSges(6,3) - pkin(3) - pkin(8)) * t172) * t160 + t193 + t200;
t44 = t161 * t176 + (-t230 + (-rSges(5,2) * t172 + rSges(5,3) * t169) * t160) * t160 + t215;
t43 = (-t160 * t76 + t161 * t77) * t172;
t40 = t175 + t187 + t68;
t39 = t222 + (-pkin(2) + (-qJ(4) - t235) * t169 + (-rSges(7,3) - pkin(3) + t174) * t172) * t160 + t192 + t200;
t36 = -t223 * t68 + t57;
t35 = t169 * t78 + t214 * t220 + t58;
t34 = t128 * t223 - t169 * t231 + t84;
t29 = t120 * t73 + t121 * t75 + t223 * t71;
t28 = t120 * t72 + t121 * t74 + t223 * t70;
t27 = t118 * t73 + t119 * t75 + t220 * t71;
t26 = t118 * t72 + t119 * t74 + t220 * t70;
t25 = t160 * t77 + t161 * t76 + t198;
t18 = t57 + (-t160 * t232 + t161 * t79) * t172;
t17 = t160 * t231 + t161 * t232 + t198;
t15 = t28 * t160 - t29 * t161;
t14 = t26 * t160 - t27 * t161;
t8 = t46 * t169 + (t160 * t29 + t161 * t28) * t172;
t7 = t45 * t169 + (t160 * t27 + t161 * t26) * t172;
t1 = [Icges(2,3) + Icges(3,3) + t112 + t87 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2) + m(4) * (t82 ^ 2 + t83 ^ 2) + m(3) * (t126 ^ 2 + t127 ^ 2) + m(2) * (t141 ^ 2 + t142 ^ 2) + (t179 + t180 + t226 + t228 + (Icges(4,2) + Icges(5,3)) * t172) * t172 + (t225 + t227 + (Icges(4,1) + Icges(5,2)) * t169) * t169; 0; m(3) + m(4) - t209; m(7) * (t56 * t39 + t55 * t40) + m(6) * (t81 * t47 + t80 * t48) + m(5) * (t86 * t60 + t85 * t61) + (-t42 / 0.2e1 - t31 / 0.2e1 - t82 * t236 + t201 * t161 + (Icges(5,5) * t238 + Icges(4,6) * t245 + t181 * t239 + t185 * t246) * t172 + (Icges(5,4) * t238 + Icges(4,5) * t245 + t182 * t239 + t186 * t246) * t169 - t206) * t161 + (t41 / 0.2e1 + t30 / 0.2e1 - t83 * t236 + t201 * t160 + (Icges(5,5) * t246 + Icges(4,6) * t239 + t181 * t238 + t185 * t245) * t172 + (Icges(5,4) * t246 + Icges(4,5) * t239 + t182 * t238 + t186 * t245) * t169 + t207) * t160; m(4) * t59 + m(5) * t44 + m(6) * t25 + m(7) * t17; m(7) * (t17 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t25 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t44 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t140 ^ 2 * t210 + t59 ^ 2) + (t244 * t159 - t13 - t15) * t161 + (t12 + t14 + t243 * t158 + (t244 * t160 + t243 * t161) * t161) * t160; 0.2e1 * ((t160 * t40 + t161 * t39) * t240 + (t160 * t48 + t161 * t47) * t241 + (t160 * t61 + t161 * t60) * t242) * t169; t209 * t172; m(7) * (-t172 * t17 + (t160 * t55 + t161 * t56) * t169) + m(6) * (-t172 * t25 + (t160 * t80 + t161 * t81) * t169) + m(5) * (-t172 * t44 + (t160 * t85 + t161 * t86) * t169); 0.2e1 * (t242 + t241 + t240) * (t169 ^ 2 * t210 + t172 ^ 2); t54 + m(7) * (t34 * t39 + t35 * t40) + m(6) * (t52 * t47 + t53 * t48) + (t160 * t206 + t161 * t207) * t172 + t196; m(6) * t43 + m(7) * t18; t7 * t239 + t8 * t238 + (t32 * t160 - t33 * t161) * t237 + (t14 * t245 + t15 * t239) * t172 + m(7) * (t18 * t17 + t34 * t56 + t35 * t55) + m(6) * (t43 * t25 + t52 * t81 + t53 * t80) + t197; m(6) * (-t43 * t172 + (t160 * t53 + t161 * t52) * t169) + m(7) * (-t18 * t172 + (t160 * t35 + t161 * t34) * t169); t169 * t54 + m(7) * (t18 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t43 ^ 2 + t52 ^ 2 + t53 ^ 2) + (t161 * t7 + t160 * t8 + t169 * (t160 * t33 + t161 * t32)) * t172 + t208; m(7) * (t49 * t39 + t50 * t40) + t196; m(7) * t36; m(7) * (t36 * t17 + t49 * t56 + t50 * t55) + t197; m(7) * (-t36 * t172 + (t160 * t50 + t161 * t49) * t169); m(7) * (t36 * t18 + t49 * t34 + t50 * t35) + t208; m(7) * (t36 ^ 2 + t49 ^ 2 + t50 ^ 2) + t208;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
