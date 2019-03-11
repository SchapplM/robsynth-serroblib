% Calculate joint inertia matrix for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:00:10
% EndTime: 2019-03-08 20:00:20
% DurationCPUTime: 5.36s
% Computational Cost: add. (23147->526), mult. (59574->771), div. (0->0), fcn. (78654->12), ass. (0->225)
t256 = rSges(7,1) + pkin(5);
t255 = rSges(7,3) + qJ(6);
t202 = sin(pkin(10));
t204 = cos(pkin(10));
t243 = sin(pkin(11));
t244 = cos(pkin(11));
t248 = sin(qJ(2));
t251 = cos(qJ(2));
t195 = -t248 * t243 + t251 * t244;
t205 = cos(pkin(6));
t208 = t205 * t195;
t209 = t243 * t251 + t248 * t244;
t173 = -t202 * t209 + t204 * t208;
t187 = t209 * t205;
t174 = t187 * t204 + t195 * t202;
t203 = sin(pkin(6));
t241 = t203 * t204;
t127 = Icges(4,5) * t174 + Icges(4,6) * t173 - Icges(4,3) * t241;
t220 = t205 * t251;
t190 = -t202 * t248 + t204 * t220;
t219 = t205 * t248;
t191 = t202 * t251 + t204 * t219;
t160 = Icges(3,5) * t191 + Icges(3,6) * t190 - Icges(3,3) * t241;
t254 = -t160 - t127;
t175 = -t202 * t208 - t204 * t209;
t176 = -t187 * t202 + t195 * t204;
t242 = t202 * t203;
t128 = Icges(4,5) * t176 + Icges(4,6) * t175 + Icges(4,3) * t242;
t192 = -t202 * t220 - t204 * t248;
t193 = -t202 * t219 + t204 * t251;
t161 = Icges(3,5) * t193 + Icges(3,6) * t192 + Icges(3,3) * t242;
t253 = t161 + t128;
t185 = t195 * t203;
t186 = t209 * t203;
t252 = Icges(4,5) * t186 + Icges(4,6) * t185 + (t248 * Icges(3,5) + Icges(3,6) * t251) * t203 + (Icges(4,3) + Icges(3,3)) * t205;
t250 = cos(qJ(4));
t249 = cos(qJ(5));
t247 = pkin(2) * t251;
t207 = sin(qJ(4));
t240 = t203 * t207;
t150 = t174 * t250 - t204 * t240;
t206 = sin(qJ(5));
t123 = t150 * t206 + t173 * t249;
t124 = t150 * t249 - t173 * t206;
t221 = t203 * t250;
t149 = t174 * t207 + t204 * t221;
t246 = rSges(7,2) * t149 + t255 * t123 + t256 * t124;
t152 = t176 * t250 + t202 * t240;
t125 = t152 * t206 + t175 * t249;
t126 = t152 * t249 - t175 * t206;
t151 = t176 * t207 - t202 * t221;
t245 = rSges(7,2) * t151 + t255 * t125 + t256 * t126;
t178 = t186 * t250 + t205 * t207;
t147 = t178 * t206 + t185 * t249;
t148 = t178 * t249 - t185 * t206;
t177 = t186 * t207 - t205 * t250;
t239 = rSges(7,2) * t177 + t255 * t147 + t256 * t148;
t139 = pkin(3) * t176 - pkin(8) * t175;
t218 = pkin(2) * t219 - qJ(3) * t203;
t171 = -t202 * t218 + t204 * t247;
t159 = t205 * t171;
t238 = t205 * t139 + t159;
t138 = pkin(3) * t174 - pkin(8) * t173;
t170 = t202 * t247 + t204 * t218;
t237 = -t138 - t170;
t236 = t170 * t242 + t171 * t241;
t196 = pkin(2) * t203 * t248 + t205 * qJ(3);
t235 = -pkin(3) * t186 + pkin(8) * t185 - t196;
t77 = Icges(7,5) * t124 + Icges(7,6) * t149 + Icges(7,3) * t123;
t81 = Icges(7,4) * t124 + Icges(7,2) * t149 + Icges(7,6) * t123;
t85 = Icges(7,1) * t124 + Icges(7,4) * t149 + Icges(7,5) * t123;
t33 = t123 * t77 + t124 * t85 + t149 * t81;
t78 = Icges(7,5) * t126 + Icges(7,6) * t151 + Icges(7,3) * t125;
t82 = Icges(7,4) * t126 + Icges(7,2) * t151 + Icges(7,6) * t125;
t86 = Icges(7,1) * t126 + Icges(7,4) * t151 + Icges(7,5) * t125;
t34 = t123 * t78 + t124 * t86 + t149 * t82;
t106 = Icges(7,5) * t148 + Icges(7,6) * t177 + Icges(7,3) * t147;
t108 = Icges(7,4) * t148 + Icges(7,2) * t177 + Icges(7,6) * t147;
t110 = Icges(7,1) * t148 + Icges(7,4) * t177 + Icges(7,5) * t147;
t52 = t106 * t123 + t108 * t149 + t110 * t124;
t1 = t149 * t33 + t151 * t34 + t177 * t52;
t79 = Icges(6,5) * t124 - Icges(6,6) * t123 + Icges(6,3) * t149;
t83 = Icges(6,4) * t124 - Icges(6,2) * t123 + Icges(6,6) * t149;
t87 = Icges(6,1) * t124 - Icges(6,4) * t123 + Icges(6,5) * t149;
t35 = -t123 * t83 + t124 * t87 + t149 * t79;
t80 = Icges(6,5) * t126 - Icges(6,6) * t125 + Icges(6,3) * t151;
t84 = Icges(6,4) * t126 - Icges(6,2) * t125 + Icges(6,6) * t151;
t88 = Icges(6,1) * t126 - Icges(6,4) * t125 + Icges(6,5) * t151;
t36 = -t123 * t84 + t124 * t88 + t149 * t80;
t107 = Icges(6,5) * t148 - Icges(6,6) * t147 + Icges(6,3) * t177;
t109 = Icges(6,4) * t148 - Icges(6,2) * t147 + Icges(6,6) * t177;
t111 = Icges(6,1) * t148 - Icges(6,4) * t147 + Icges(6,5) * t177;
t53 = t107 * t149 - t109 * t123 + t111 * t124;
t2 = t149 * t35 + t151 * t36 + t177 * t53;
t234 = -t2 / 0.2e1 - t1 / 0.2e1;
t37 = t125 * t77 + t126 * t85 + t151 * t81;
t38 = t125 * t78 + t126 * t86 + t151 * t82;
t54 = t106 * t125 + t108 * t151 + t110 * t126;
t3 = t149 * t37 + t151 * t38 + t177 * t54;
t39 = -t125 * t83 + t126 * t87 + t151 * t79;
t40 = -t125 * t84 + t126 * t88 + t151 * t80;
t55 = t107 * t151 - t109 * t125 + t111 * t126;
t4 = t149 * t39 + t151 * t40 + t177 * t55;
t233 = -t4 / 0.2e1 - t3 / 0.2e1;
t5 = -t173 * t33 - t175 * t34 - t185 * t52;
t6 = -t173 * t35 - t175 * t36 - t185 * t53;
t232 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = -t173 * t37 - t175 * t38 - t185 * t54;
t8 = -t173 * t39 - t175 * t40 - t185 * t55;
t231 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t205 * t53 + (t202 * t36 - t204 * t35) * t203;
t9 = t205 * t52 + (t202 * t34 - t204 * t33) * t203;
t230 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t205 * t54 + (t202 * t38 - t204 * t37) * t203;
t12 = t205 * t55 + (t202 * t40 - t204 * t39) * t203;
t229 = t12 / 0.2e1 + t11 / 0.2e1;
t41 = t147 * t77 + t148 * t85 + t177 * t81;
t42 = t147 * t78 + t148 * t86 + t177 * t82;
t62 = t106 * t147 + t108 * t177 + t110 * t148;
t13 = t149 * t41 + t151 * t42 + t177 * t62;
t43 = -t147 * t83 + t148 * t87 + t177 * t79;
t44 = -t147 * t84 + t148 * t88 + t177 * t80;
t63 = t107 * t177 - t109 * t147 + t111 * t148;
t14 = t149 * t43 + t151 * t44 + t177 * t63;
t228 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = -t173 * t41 - t175 * t42 - t185 * t62;
t16 = -t173 * t43 - t175 * t44 - t185 * t63;
t227 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t205 * t62 + (t202 * t42 - t204 * t41) * t203;
t18 = t205 * t63 + (t202 * t44 - t204 * t43) * t203;
t226 = t18 / 0.2e1 + t17 / 0.2e1;
t225 = m(4) + m(5) + m(6) + m(7);
t120 = pkin(4) * t152 + pkin(9) * t151;
t224 = t205 * t120 + t238;
t119 = pkin(4) * t150 + pkin(9) * t149;
t223 = -t119 + t237;
t146 = pkin(4) * t178 + pkin(9) * t177;
t222 = -t146 + t235;
t217 = (-rSges(4,1) * t186 - rSges(4,2) * t185 - rSges(4,3) * t205 - t196) * t203;
t216 = t138 * t242 + t139 * t241 + t236;
t143 = rSges(5,1) * t178 - rSges(5,2) * t177 - rSges(5,3) * t185;
t215 = (-t143 + t235) * t203;
t113 = rSges(6,1) * t148 - rSges(6,2) * t147 + rSges(6,3) * t177;
t212 = (-t113 + t222) * t203;
t211 = t119 * t242 + t120 * t241 + t216;
t210 = (t222 - t239) * t203;
t184 = t205 * rSges(3,3) + (t248 * rSges(3,1) + rSges(3,2) * t251) * t203;
t183 = Icges(3,5) * t205 + (t248 * Icges(3,1) + Icges(3,4) * t251) * t203;
t182 = Icges(3,6) * t205 + (t248 * Icges(3,4) + Icges(3,2) * t251) * t203;
t167 = rSges(3,1) * t193 + rSges(3,2) * t192 + rSges(3,3) * t242;
t166 = rSges(3,1) * t191 + rSges(3,2) * t190 - rSges(3,3) * t241;
t165 = Icges(3,1) * t193 + Icges(3,4) * t192 + Icges(3,5) * t242;
t164 = Icges(3,1) * t191 + Icges(3,4) * t190 - Icges(3,5) * t241;
t163 = Icges(3,4) * t193 + Icges(3,2) * t192 + Icges(3,6) * t242;
t162 = Icges(3,4) * t191 + Icges(3,2) * t190 - Icges(3,6) * t241;
t157 = Icges(4,1) * t186 + Icges(4,4) * t185 + Icges(4,5) * t205;
t156 = Icges(4,4) * t186 + Icges(4,2) * t185 + Icges(4,6) * t205;
t145 = -t166 * t205 - t184 * t241;
t144 = t167 * t205 - t184 * t242;
t142 = Icges(5,1) * t178 - Icges(5,4) * t177 - Icges(5,5) * t185;
t141 = Icges(5,4) * t178 - Icges(5,2) * t177 - Icges(5,6) * t185;
t140 = Icges(5,5) * t178 - Icges(5,6) * t177 - Icges(5,3) * t185;
t134 = rSges(4,1) * t176 + rSges(4,2) * t175 + rSges(4,3) * t242;
t133 = rSges(4,1) * t174 + rSges(4,2) * t173 - rSges(4,3) * t241;
t132 = Icges(4,1) * t176 + Icges(4,4) * t175 + Icges(4,5) * t242;
t131 = Icges(4,1) * t174 + Icges(4,4) * t173 - Icges(4,5) * t241;
t130 = Icges(4,4) * t176 + Icges(4,2) * t175 + Icges(4,6) * t242;
t129 = Icges(4,4) * t174 + Icges(4,2) * t173 - Icges(4,6) * t241;
t122 = (t166 * t202 + t167 * t204) * t203;
t121 = t173 * t146;
t114 = t185 * t120;
t105 = t175 * t119;
t104 = rSges(5,1) * t152 - rSges(5,2) * t151 - rSges(5,3) * t175;
t103 = rSges(5,1) * t150 - rSges(5,2) * t149 - rSges(5,3) * t173;
t102 = Icges(5,1) * t152 - Icges(5,4) * t151 - Icges(5,5) * t175;
t101 = Icges(5,1) * t150 - Icges(5,4) * t149 - Icges(5,5) * t173;
t100 = Icges(5,4) * t152 - Icges(5,2) * t151 - Icges(5,6) * t175;
t99 = Icges(5,4) * t150 - Icges(5,2) * t149 - Icges(5,6) * t173;
t98 = Icges(5,5) * t152 - Icges(5,6) * t151 - Icges(5,3) * t175;
t97 = Icges(5,5) * t150 - Icges(5,6) * t149 - Icges(5,3) * t173;
t94 = (-t133 - t170) * t205 + t204 * t217;
t93 = t134 * t205 + t202 * t217 + t159;
t92 = rSges(6,1) * t126 - rSges(6,2) * t125 + rSges(6,3) * t151;
t90 = rSges(6,1) * t124 - rSges(6,2) * t123 + rSges(6,3) * t149;
t76 = (t133 * t202 + t134 * t204) * t203 + t236;
t75 = -t104 * t185 + t143 * t175;
t74 = t103 * t185 - t143 * t173;
t73 = -t140 * t185 - t141 * t177 + t142 * t178;
t72 = -t103 * t175 + t104 * t173;
t71 = -t140 * t175 - t141 * t151 + t142 * t152;
t70 = -t140 * t173 - t141 * t149 + t142 * t150;
t69 = (-t103 + t237) * t205 + t204 * t215;
t68 = t104 * t205 + t202 * t215 + t238;
t67 = -t113 * t151 + t177 * t92;
t66 = t113 * t149 - t177 * t90;
t65 = -t100 * t177 + t102 * t178 - t185 * t98;
t64 = t101 * t178 - t177 * t99 - t185 * t97;
t61 = (t103 * t202 + t104 * t204) * t203 + t216;
t60 = -t100 * t151 + t102 * t152 - t175 * t98;
t59 = t101 * t152 - t151 * t99 - t175 * t97;
t58 = -t100 * t149 + t102 * t150 - t173 * t98;
t57 = t101 * t150 - t149 * t99 - t173 * t97;
t56 = -t149 * t92 + t151 * t90;
t51 = -t185 * t92 - t114 + (t113 + t146) * t175;
t50 = -t113 * t173 - t121 - (-t119 - t90) * t185;
t49 = (-t90 + t223) * t205 + t204 * t212;
t48 = t202 * t212 + t205 * t92 + t224;
t47 = -t151 * t239 + t177 * t245;
t46 = t149 * t239 - t177 * t246;
t45 = -t175 * t90 - t105 + (t120 + t92) * t173;
t32 = (t202 * t90 + t204 * t92) * t203 + t211;
t31 = -t114 - t245 * t185 + (t146 + t239) * t175;
t30 = -t121 - t239 * t173 - (-t119 - t246) * t185;
t29 = (t223 - t246) * t205 + t204 * t210;
t28 = t202 * t210 + t205 * t245 + t224;
t27 = -t149 * t245 + t151 * t246;
t26 = -t105 - t246 * t175 + (t120 + t245) * t173;
t25 = t205 * t73 + (t202 * t65 - t204 * t64) * t203;
t24 = (t202 * t246 + t204 * t245) * t203 + t211;
t23 = -t173 * t64 - t175 * t65 - t185 * t73;
t22 = t205 * t71 + (t202 * t60 - t204 * t59) * t203;
t21 = t205 * t70 + (t202 * t58 - t204 * t57) * t203;
t20 = -t173 * t59 - t175 * t60 - t185 * t71;
t19 = -t173 * t57 - t175 * t58 - t185 * t70;
t89 = [m(3) + m(2) + t225; m(3) * t122 + m(4) * t76 + m(5) * t61 + m(6) * t32 + m(7) * t24; m(7) * (t24 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t61 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t76 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t122 ^ 2 + t144 ^ 2 + t145 ^ 2) + (t17 + t18 + t25 + (t185 * t156 + t186 * t157 + t252 * t205) * t205 + ((t185 * t130 + t186 * t132) * t202 - (t185 * t129 + t186 * t131) * t204 + (-t127 * t204 + t128 * t202 + t182 * t251 + t248 * t183) * t205) * t203) * t205 + (t12 + t11 + t22 + (t175 * t130 + t176 * t132 + t192 * t163 + t193 * t165 + t253 * t242) * t242 + ((t163 * t251 + t248 * t165) * t203 + t192 * t182 + t193 * t183 + t175 * t156 + t176 * t157 + t252 * t242 + t205 * t161) * t205) * t242 + (-t10 - t9 - t21 + (t173 * t129 + t174 * t131 + t190 * t162 + t191 * t164 + t254 * t241) * t241 + (-(t162 * t251 + t248 * t164) * t203 - t190 * t182 - t191 * t183 - t173 * t156 - t174 * t157 + t252 * t241 - t205 * t160) * t205 + (-t175 * t129 - t173 * t130 - t176 * t131 - t174 * t132 - t192 * t162 - t190 * t163 - t193 * t164 - t191 * t165 + t253 * t241 + t254 * t242) * t242) * t241; t225 * t205; m(7) * (t205 * t24 + (t202 * t29 - t204 * t28) * t203) + m(6) * (t205 * t32 + (t202 * t49 - t204 * t48) * t203) + m(5) * (t205 * t61 + (t202 * t69 - t204 * t68) * t203) + m(4) * (t205 * t76 + (t202 * t94 - t204 * t93) * t203); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t205 ^ 2 + (t202 ^ 2 + t204 ^ 2) * t203 ^ 2); m(5) * t72 + m(6) * t45 + m(7) * t26; (t23 / 0.2e1 + t227) * t205 - (t25 / 0.2e1 + t226) * t185 + (-t22 / 0.2e1 - t229) * t175 + (-t21 / 0.2e1 - t230) * t173 + m(7) * (t24 * t26 + t28 * t31 + t29 * t30) + m(6) * (t32 * t45 + t48 * t51 + t49 * t50) + m(5) * (t61 * t72 + t68 * t75 + t69 * t74) + ((-t19 / 0.2e1 - t232) * t204 + (t20 / 0.2e1 + t231) * t202) * t203; m(5) * (t205 * t72 + (t202 * t74 - t204 * t75) * t203) + m(6) * (t205 * t45 + (t202 * t50 - t204 * t51) * t203) + m(7) * (t205 * t26 + (t202 * t30 - t204 * t31) * t203); -(t16 + t15 + t23) * t185 + (-t7 - t8 - t20) * t175 + (-t5 - t6 - t19) * t173 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t45 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t72 ^ 2 + t74 ^ 2 + t75 ^ 2); m(6) * t56 + m(7) * t27; t228 * t205 + t226 * t177 + t229 * t151 + t230 * t149 + m(7) * (t24 * t27 + t28 * t47 + t29 * t46) + m(6) * (t32 * t56 + t48 * t67 + t49 * t66) + (-t202 * t233 + t204 * t234) * t203; m(6) * (t205 * t56 + (t202 * t66 - t204 * t67) * t203) + m(7) * (t205 * t27 + (t202 * t46 - t204 * t47) * t203); -t228 * t185 + t227 * t177 + t233 * t175 + t234 * t173 + t231 * t151 + t232 * t149 + m(7) * (t26 * t27 + t30 * t46 + t31 * t47) + m(6) * (t45 * t56 + t50 * t66 + t51 * t67); (t13 + t14) * t177 + (t3 + t4) * t151 + (t1 + t2) * t149 + m(7) * (t27 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t56 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * t147; m(7) * (t123 * t28 + t125 * t29 + t147 * t24); m(7) * (t147 * t205 + (-t123 * t204 + t125 * t202) * t203); m(7) * (t123 * t31 + t125 * t30 + t147 * t26); m(7) * (t123 * t47 + t125 * t46 + t147 * t27); m(7) * (t123 ^ 2 + t125 ^ 2 + t147 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t89(1) t89(2) t89(4) t89(7) t89(11) t89(16); t89(2) t89(3) t89(5) t89(8) t89(12) t89(17); t89(4) t89(5) t89(6) t89(9) t89(13) t89(18); t89(7) t89(8) t89(9) t89(10) t89(14) t89(19); t89(11) t89(12) t89(13) t89(14) t89(15) t89(20); t89(16) t89(17) t89(18) t89(19) t89(20) t89(21);];
Mq  = res;
