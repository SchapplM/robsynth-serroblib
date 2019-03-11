% Calculate joint inertia matrix for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:01
% EndTime: 2019-03-08 19:43:05
% DurationCPUTime: 3.59s
% Computational Cost: add. (14634->511), mult. (24875->755), div. (0->0), fcn. (31111->12), ass. (0->221)
t192 = sin(pkin(6));
t245 = t192 ^ 2;
t218 = m(6) / 0.2e1 + m(7) / 0.2e1;
t244 = 0.2e1 * t218;
t190 = sin(pkin(11));
t193 = cos(pkin(11));
t195 = cos(pkin(6));
t198 = sin(qJ(2));
t230 = t192 * t198;
t178 = -t190 * t230 + t193 * t195;
t233 = t190 * t195;
t179 = t193 * t230 + t233;
t200 = cos(qJ(2));
t229 = t192 * t200;
t141 = Icges(4,5) * t179 + Icges(4,6) * t178 - Icges(4,3) * t229;
t173 = Icges(3,6) * t195 + (Icges(3,4) * t198 + Icges(3,2) * t200) * t192;
t243 = -t173 + t141;
t191 = sin(pkin(10));
t194 = cos(pkin(10));
t228 = t194 * t198;
t181 = t191 * t200 + t195 * t228;
t219 = pkin(11) + qJ(4);
t209 = sin(t219);
t204 = t192 * t209;
t210 = cos(t219);
t162 = t181 * t210 - t194 * t204;
t226 = t198 * t191;
t183 = t194 * t200 - t195 * t226;
t164 = t183 * t210 + t191 * t204;
t205 = t192 * t210;
t176 = t195 * t209 + t198 * t205;
t163 = t183 * t209 - t191 * t205;
t227 = t195 * t200;
t182 = t191 * t227 + t228;
t197 = sin(qJ(6));
t199 = cos(qJ(6));
t129 = t163 * t199 - t182 * t197;
t130 = t163 * t197 + t182 * t199;
t161 = t181 * t209 + t194 * t205;
t180 = -t194 * t227 + t226;
t127 = t161 * t199 - t180 * t197;
t128 = t161 * t197 + t180 * t199;
t73 = Icges(7,5) * t128 + Icges(7,6) * t127 + Icges(7,3) * t162;
t75 = Icges(7,4) * t128 + Icges(7,2) * t127 + Icges(7,6) * t162;
t77 = Icges(7,1) * t128 + Icges(7,4) * t127 + Icges(7,5) * t162;
t29 = t129 * t75 + t130 * t77 + t164 * t73;
t74 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t164;
t76 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t164;
t78 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t164;
t30 = t129 * t76 + t130 * t78 + t164 * t74;
t175 = -t195 * t210 + t198 * t204;
t165 = t175 * t199 + t197 * t229;
t166 = t175 * t197 - t199 * t229;
t85 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t176;
t86 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t176;
t87 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t176;
t39 = t129 * t86 + t130 * t87 + t164 * t85;
t2 = t162 * t29 + t164 * t30 + t176 * t39;
t242 = t2 / 0.2e1;
t241 = t162 / 0.2e1;
t240 = t164 / 0.2e1;
t239 = t176 / 0.2e1;
t238 = pkin(3) * t193;
t79 = rSges(7,1) * t128 + rSges(7,2) * t127 + rSges(7,3) * t162;
t237 = pkin(5) * t180 + pkin(9) * t162 + t79;
t80 = rSges(7,1) * t130 + rSges(7,2) * t129 + rSges(7,3) * t164;
t236 = pkin(5) * t182 + pkin(9) * t164 + t80;
t232 = t191 * t192;
t214 = t190 * t232;
t108 = pkin(3) * t214 + pkin(8) * t182 + t183 * t238;
t160 = pkin(2) * t183 + qJ(3) * t182;
t158 = t195 * t160;
t235 = t195 * t108 + t158;
t89 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t176;
t234 = -pkin(5) * t229 + t176 * pkin(9) + t89;
t231 = t192 * t194;
t104 = rSges(6,1) * t182 - rSges(6,2) * t164 + rSges(6,3) * t163;
t125 = pkin(4) * t164 + qJ(5) * t163;
t224 = -t104 - t125;
t213 = t190 * t231;
t107 = -pkin(3) * t213 + pkin(8) * t180 + t181 * t238;
t159 = pkin(2) * t181 + qJ(3) * t180;
t223 = -t107 - t159;
t124 = pkin(4) * t162 + qJ(5) * t161;
t144 = pkin(4) * t176 + qJ(5) * t175;
t222 = t124 * t229 + t180 * t144;
t184 = (pkin(2) * t198 - qJ(3) * t200) * t192;
t221 = -pkin(3) * t233 - (-pkin(8) * t200 + t198 * t238) * t192 - t184;
t220 = t159 * t232 + t160 * t231;
t217 = t195 * t125 + t235;
t216 = -t125 - t236;
t215 = -m(4) - m(5) - m(6) - m(7);
t212 = -t124 + t223;
t211 = -t144 + t221;
t208 = (-t179 * rSges(4,1) - t178 * rSges(4,2) + rSges(4,3) * t229 - t184) * t192;
t207 = t107 * t232 + t108 * t231 + t220;
t140 = t176 * rSges(5,1) - t175 * rSges(5,2) - rSges(5,3) * t229;
t206 = (-t140 + t221) * t192;
t139 = -rSges(6,1) * t229 - t176 * rSges(6,2) + t175 * rSges(6,3);
t203 = (-t139 + t211) * t192;
t202 = t124 * t232 + t125 * t231 + t207;
t201 = (t211 - t234) * t192;
t177 = t195 * rSges(3,3) + (rSges(3,1) * t198 + rSges(3,2) * t200) * t192;
t174 = Icges(3,5) * t195 + (Icges(3,1) * t198 + Icges(3,4) * t200) * t192;
t172 = Icges(3,3) * t195 + (Icges(3,5) * t198 + Icges(3,6) * t200) * t192;
t170 = t183 * t193 + t214;
t169 = -t183 * t190 + t193 * t232;
t168 = t181 * t193 - t213;
t167 = -t181 * t190 - t193 * t231;
t154 = rSges(3,1) * t183 - rSges(3,2) * t182 + rSges(3,3) * t232;
t153 = rSges(3,1) * t181 - rSges(3,2) * t180 - rSges(3,3) * t231;
t150 = Icges(3,1) * t183 - Icges(3,4) * t182 + Icges(3,5) * t232;
t149 = Icges(3,1) * t181 - Icges(3,4) * t180 - Icges(3,5) * t231;
t148 = Icges(3,4) * t183 - Icges(3,2) * t182 + Icges(3,6) * t232;
t147 = Icges(3,4) * t181 - Icges(3,2) * t180 - Icges(3,6) * t231;
t146 = Icges(3,5) * t183 - Icges(3,6) * t182 + Icges(3,3) * t232;
t145 = Icges(3,5) * t181 - Icges(3,6) * t180 - Icges(3,3) * t231;
t143 = Icges(4,1) * t179 + Icges(4,4) * t178 - Icges(4,5) * t229;
t142 = Icges(4,4) * t179 + Icges(4,2) * t178 - Icges(4,6) * t229;
t138 = Icges(5,1) * t176 - Icges(5,4) * t175 - Icges(5,5) * t229;
t137 = Icges(5,4) * t176 - Icges(5,2) * t175 - Icges(5,6) * t229;
t136 = Icges(5,5) * t176 - Icges(5,6) * t175 - Icges(5,3) * t229;
t135 = -Icges(6,1) * t229 - Icges(6,4) * t176 + Icges(6,5) * t175;
t134 = -Icges(6,4) * t229 - Icges(6,2) * t176 + Icges(6,6) * t175;
t133 = -Icges(6,5) * t229 - Icges(6,6) * t176 + Icges(6,3) * t175;
t123 = -t153 * t195 - t177 * t231;
t122 = t154 * t195 - t177 * t232;
t117 = rSges(4,1) * t170 + rSges(4,2) * t169 + rSges(4,3) * t182;
t116 = rSges(4,1) * t168 + rSges(4,2) * t167 + rSges(4,3) * t180;
t115 = Icges(4,1) * t170 + Icges(4,4) * t169 + Icges(4,5) * t182;
t114 = Icges(4,1) * t168 + Icges(4,4) * t167 + Icges(4,5) * t180;
t113 = Icges(4,4) * t170 + Icges(4,2) * t169 + Icges(4,6) * t182;
t112 = Icges(4,4) * t168 + Icges(4,2) * t167 + Icges(4,6) * t180;
t111 = Icges(4,5) * t170 + Icges(4,6) * t169 + Icges(4,3) * t182;
t110 = Icges(4,5) * t168 + Icges(4,6) * t167 + Icges(4,3) * t180;
t109 = t182 * t124;
t106 = rSges(5,1) * t164 - rSges(5,2) * t163 + rSges(5,3) * t182;
t105 = rSges(5,1) * t162 - rSges(5,2) * t161 + rSges(5,3) * t180;
t103 = rSges(6,1) * t180 - rSges(6,2) * t162 + rSges(6,3) * t161;
t102 = Icges(5,1) * t164 - Icges(5,4) * t163 + Icges(5,5) * t182;
t101 = Icges(5,1) * t162 - Icges(5,4) * t161 + Icges(5,5) * t180;
t100 = Icges(6,1) * t182 - Icges(6,4) * t164 + Icges(6,5) * t163;
t99 = Icges(6,1) * t180 - Icges(6,4) * t162 + Icges(6,5) * t161;
t98 = Icges(5,4) * t164 - Icges(5,2) * t163 + Icges(5,6) * t182;
t97 = Icges(5,4) * t162 - Icges(5,2) * t161 + Icges(5,6) * t180;
t96 = Icges(6,4) * t182 - Icges(6,2) * t164 + Icges(6,6) * t163;
t95 = Icges(6,4) * t180 - Icges(6,2) * t162 + Icges(6,6) * t161;
t94 = Icges(5,5) * t164 - Icges(5,6) * t163 + Icges(5,3) * t182;
t93 = Icges(5,5) * t162 - Icges(5,6) * t161 + Icges(5,3) * t180;
t92 = Icges(6,5) * t182 - Icges(6,6) * t164 + Icges(6,3) * t163;
t91 = Icges(6,5) * t180 - Icges(6,6) * t162 + Icges(6,3) * t161;
t82 = (t153 * t191 + t154 * t194) * t192;
t72 = -t106 * t229 - t182 * t140;
t71 = t105 * t229 + t180 * t140;
t70 = (-t116 - t159) * t195 + t194 * t208;
t69 = t117 * t195 + t191 * t208 + t158;
t68 = -t136 * t229 - t175 * t137 + t176 * t138;
t67 = t175 * t133 - t176 * t134 - t135 * t229;
t66 = t105 * t182 - t106 * t180;
t65 = t136 * t182 - t137 * t163 + t138 * t164;
t64 = t136 * t180 - t137 * t161 + t138 * t162;
t63 = t133 * t163 - t134 * t164 + t135 * t182;
t62 = t133 * t161 - t134 * t162 + t135 * t180;
t61 = (t116 * t191 + t117 * t194) * t192 + t220;
t60 = -t164 * t89 + t176 * t80;
t59 = t162 * t89 - t176 * t79;
t58 = t224 * t229 + (-t139 - t144) * t182;
t57 = t103 * t229 + t180 * t139 + t222;
t56 = t176 * t102 - t175 * t98 - t229 * t94;
t55 = t176 * t101 - t175 * t97 - t229 * t93;
t54 = -t100 * t229 + t175 * t92 - t176 * t96;
t53 = t175 * t91 - t176 * t95 - t229 * t99;
t52 = t102 * t164 - t163 * t98 + t182 * t94;
t51 = t101 * t164 - t163 * t97 + t182 * t93;
t50 = t102 * t162 - t161 * t98 + t180 * t94;
t49 = t101 * t162 - t161 * t97 + t180 * t93;
t48 = t100 * t182 + t163 * t92 - t164 * t96;
t47 = t163 * t91 - t164 * t95 + t182 * t99;
t46 = t100 * t180 + t161 * t92 - t162 * t96;
t45 = t161 * t91 - t162 * t95 + t180 * t99;
t44 = (-t105 + t223) * t195 + t194 * t206;
t43 = t106 * t195 + t191 * t206 + t235;
t42 = t165 * t86 + t166 * t87 + t176 * t85;
t41 = -t162 * t80 + t164 * t79;
t40 = t103 * t182 + t180 * t224 + t109;
t38 = t127 * t86 + t128 * t87 + t162 * t85;
t37 = (t105 * t191 + t106 * t194) * t192 + t207;
t36 = (-t103 + t212) * t195 + t194 * t203;
t35 = t104 * t195 + t191 * t203 + t217;
t34 = t216 * t229 + (-t144 - t234) * t182;
t33 = t180 * t234 + t229 * t237 + t222;
t32 = t165 * t76 + t166 * t78 + t176 * t74;
t31 = t165 * t75 + t166 * t77 + t176 * t73;
t28 = t127 * t76 + t128 * t78 + t162 * t74;
t27 = t127 * t75 + t128 * t77 + t162 * t73;
t26 = (t103 * t191 + t104 * t194) * t192 + t202;
t25 = t180 * t216 + t182 * t237 + t109;
t24 = (t212 - t237) * t195 + t194 * t201;
t23 = t191 * t201 + t195 * t236 + t217;
t22 = t195 * t68 + (t191 * t56 - t194 * t55) * t192;
t21 = t195 * t67 + (t191 * t54 - t194 * t53) * t192;
t20 = (t191 * t237 + t194 * t236) * t192 + t202;
t19 = t55 * t180 + t56 * t182 - t229 * t68;
t18 = t53 * t180 + t54 * t182 - t229 * t67;
t17 = t195 * t65 + (t191 * t52 - t194 * t51) * t192;
t16 = t195 * t64 + (t191 * t50 - t194 * t49) * t192;
t15 = t195 * t63 + (t191 * t48 - t194 * t47) * t192;
t14 = t195 * t62 + (t191 * t46 - t194 * t45) * t192;
t13 = t51 * t180 + t52 * t182 - t229 * t65;
t12 = t49 * t180 + t50 * t182 - t229 * t64;
t11 = t47 * t180 + t48 * t182 - t229 * t63;
t10 = t45 * t180 + t46 * t182 - t229 * t62;
t9 = t195 * t42 + (t191 * t32 - t194 * t31) * t192;
t8 = t31 * t180 + t32 * t182 - t229 * t42;
t7 = t162 * t31 + t164 * t32 + t176 * t42;
t6 = t195 * t39 + (t191 * t30 - t194 * t29) * t192;
t5 = t195 * t38 + (t191 * t28 - t194 * t27) * t192;
t4 = t29 * t180 + t30 * t182 - t229 * t39;
t3 = t27 * t180 + t28 * t182 - t229 * t38;
t1 = t162 * t27 + t164 * t28 + t176 * t38;
t81 = [m(2) + m(3) - t215; m(3) * t82 + m(4) * t61 + m(5) * t37 + m(6) * t26 + m(7) * t20; m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t26 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t37 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(4) * (t61 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(3) * (t122 ^ 2 + t123 ^ 2 + t82 ^ 2) + (t6 + t17 + t15 + ((t111 * t182 + t113 * t169 + t115 * t170) * t191 - (t110 * t182 + t112 * t169 + t114 * t170) * t194) * t192 + (t146 * t232 - t148 * t182 + t150 * t183) * t232) * t232 + (-t5 - t16 - t14 - ((t111 * t180 + t113 * t167 + t115 * t168) * t191 - (t110 * t180 + t112 * t167 + t114 * t168) * t194) * t192 + (-t145 * t231 - t147 * t180 + t149 * t181) * t231 + (-t145 * t232 + t146 * t231 + t147 * t182 + t148 * t180 - t149 * t183 - t150 * t181) * t232) * t231 + (t9 + t22 + t21 + ((t148 * t200 + t150 * t198) * t191 - (t147 * t200 + t149 * t198) * t194) * t245 + ((-t145 * t194 + t146 * t191 + t173 * t200 + t174 * t198) * t192 - t141 * t229 + t178 * t142 + t179 * t143 + t195 * t172) * t195 + (-t111 * t229 + t178 * t113 + t179 * t115 + t142 * t169 + t143 * t170 + t172 * t232 + t174 * t183 + t182 * t243) * t232 + (t110 * t229 - t178 * t112 - t179 * t114 - t142 * t167 - t143 * t168 + t172 * t231 - t174 * t181 - t180 * t243) * t231) * t195; t215 * t229; m(7) * (t180 * t23 + t182 * t24 - t20 * t229) + m(6) * (t180 * t35 + t182 * t36 - t229 * t26) + m(5) * (t180 * t43 + t182 * t44 - t229 * t37) + m(4) * (t180 * t69 + t182 * t70 - t229 * t61); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t218) * (t200 ^ 2 * t245 + t180 ^ 2 + t182 ^ 2); m(5) * t66 + m(6) * t40 + m(7) * t25; (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t195 + (t6 / 0.2e1 + t17 / 0.2e1 + t15 / 0.2e1) * t182 + (t5 / 0.2e1 + t16 / 0.2e1 + t14 / 0.2e1) * t180 + m(7) * (t20 * t25 + t23 * t34 + t24 * t33) + m(6) * (t26 * t40 + t35 * t58 + t36 * t57) + m(5) * (t37 * t66 + t43 * t72 + t44 * t71) + ((-t9 / 0.2e1 - t22 / 0.2e1 - t21 / 0.2e1) * t200 + (-t3 / 0.2e1 - t12 / 0.2e1 - t10 / 0.2e1) * t194 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t191) * t192; m(5) * (t72 * t180 + t71 * t182 - t229 * t66) + m(6) * (t58 * t180 + t57 * t182 - t229 * t40) + m(7) * (t34 * t180 + t33 * t182 - t229 * t25); (-t18 - t19 - t8) * t229 + (t4 + t11 + t13) * t182 + (t3 + t10 + t12) * t180 + m(7) * (t25 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t40 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t66 ^ 2 + t71 ^ 2 + t72 ^ 2); t175 * t244; m(7) * (t161 * t23 + t163 * t24 + t175 * t20) + m(6) * (t161 * t35 + t163 * t36 + t175 * t26); (t161 * t180 + t163 * t182 - t175 * t229) * t244; m(7) * (t161 * t34 + t163 * t33 + t175 * t25) + m(6) * (t161 * t58 + t163 * t57 + t175 * t40); (t161 ^ 2 + t163 ^ 2 + t175 ^ 2) * t244; m(7) * t41; t6 * t240 + t195 * t7 / 0.2e1 + t9 * t239 + t5 * t241 + m(7) * (t20 * t41 + t23 * t60 + t24 * t59) + (-t194 * t1 / 0.2e1 + t191 * t242) * t192; m(7) * (t60 * t180 + t59 * t182 - t229 * t41); m(7) * (t25 * t41 + t33 * t59 + t34 * t60) + t182 * t242 + t3 * t241 + t8 * t239 + t4 * t240 + t180 * t1 / 0.2e1 - t7 * t229 / 0.2e1; m(7) * (t161 * t60 + t163 * t59 + t175 * t41); t164 * t2 + t162 * t1 + t176 * t7 + m(7) * (t41 ^ 2 + t59 ^ 2 + t60 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t81(1) t81(2) t81(4) t81(7) t81(11) t81(16); t81(2) t81(3) t81(5) t81(8) t81(12) t81(17); t81(4) t81(5) t81(6) t81(9) t81(13) t81(18); t81(7) t81(8) t81(9) t81(10) t81(14) t81(19); t81(11) t81(12) t81(13) t81(14) t81(15) t81(20); t81(16) t81(17) t81(18) t81(19) t81(20) t81(21);];
Mq  = res;
