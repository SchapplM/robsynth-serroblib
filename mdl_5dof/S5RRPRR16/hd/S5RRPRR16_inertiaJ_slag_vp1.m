% Calculate joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR16_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR16_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:44
% EndTime: 2019-12-31 20:44:53
% DurationCPUTime: 3.41s
% Computational Cost: add. (8266->441), mult. (21140->631), div. (0->0), fcn. (26545->10), ass. (0->207)
t194 = sin(pkin(5));
t195 = cos(pkin(5));
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t158 = Icges(3,3) * t195 + (Icges(3,5) * t198 + Icges(3,6) * t201) * t194;
t159 = Icges(3,6) * t195 + (Icges(3,4) * t198 + Icges(3,2) * t201) * t194;
t160 = Icges(3,5) * t195 + (Icges(3,1) * t198 + Icges(3,4) * t201) * t194;
t161 = Icges(4,5) * t195 + (-Icges(4,6) * t198 - Icges(4,3) * t201) * t194;
t162 = Icges(4,4) * t195 + (-Icges(4,2) * t198 - Icges(4,6) * t201) * t194;
t163 = Icges(4,1) * t195 + (-Icges(4,4) * t198 - Icges(4,5) * t201) * t194;
t235 = t194 * t201;
t236 = t194 * t198;
t254 = (-t201 * t161 - t198 * t162) * t194 + t159 * t235 + t160 * t236 + (t163 + t158) * t195;
t202 = cos(qJ(1));
t230 = t201 * t202;
t199 = sin(qJ(1));
t234 = t198 * t199;
t177 = -t195 * t230 + t234;
t231 = t199 * t201;
t233 = t198 * t202;
t178 = t195 * t233 + t231;
t229 = t202 * t194;
t117 = -Icges(4,5) * t229 - Icges(4,6) * t178 + Icges(4,3) * t177;
t124 = Icges(3,4) * t178 - Icges(3,2) * t177 - Icges(3,6) * t229;
t253 = t117 - t124;
t119 = -Icges(4,4) * t229 - Icges(4,2) * t178 + Icges(4,6) * t177;
t126 = Icges(3,1) * t178 - Icges(3,4) * t177 - Icges(3,5) * t229;
t252 = t119 - t126;
t179 = t195 * t231 + t233;
t180 = -t195 * t234 + t230;
t232 = t199 * t194;
t116 = Icges(4,5) * t232 - Icges(4,6) * t180 + Icges(4,3) * t179;
t125 = Icges(3,4) * t180 - Icges(3,2) * t179 + Icges(3,6) * t232;
t251 = -t125 + t116;
t118 = Icges(4,4) * t232 - Icges(4,2) * t180 + Icges(4,6) * t179;
t127 = Icges(3,1) * t180 - Icges(3,4) * t179 + Icges(3,5) * t232;
t250 = t127 - t118;
t249 = t194 ^ 2;
t197 = sin(qJ(4));
t243 = cos(qJ(4));
t148 = -t179 * t243 + t197 * t232;
t150 = t177 * t243 + t197 * t229;
t218 = t194 * t243;
t149 = t179 * t197 + t199 * t218;
t196 = sin(qJ(5));
t200 = cos(qJ(5));
t107 = -t149 * t196 + t180 * t200;
t108 = t149 * t200 + t180 * t196;
t63 = Icges(6,5) * t108 + Icges(6,6) * t107 + Icges(6,3) * t148;
t65 = Icges(6,4) * t108 + Icges(6,2) * t107 + Icges(6,6) * t148;
t67 = Icges(6,1) * t108 + Icges(6,4) * t107 + Icges(6,5) * t148;
t17 = t107 * t65 + t108 * t67 + t148 * t63;
t175 = t195 * t197 + t201 * t218;
t151 = t177 * t197 - t202 * t218;
t109 = -t151 * t196 + t178 * t200;
t110 = t151 * t200 + t178 * t196;
t64 = Icges(6,5) * t110 + Icges(6,6) * t109 - Icges(6,3) * t150;
t66 = Icges(6,4) * t110 + Icges(6,2) * t109 - Icges(6,6) * t150;
t68 = Icges(6,1) * t110 + Icges(6,4) * t109 - Icges(6,5) * t150;
t18 = t107 * t66 + t108 * t68 + t148 * t64;
t176 = t195 * t243 - t197 * t235;
t145 = -t176 * t196 + t200 * t236;
t146 = t176 * t200 + t196 * t236;
t85 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t175;
t86 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t175;
t87 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t175;
t26 = t107 * t86 + t108 * t87 + t148 * t85;
t1 = t148 * t17 - t150 * t18 + t175 * t26;
t248 = t1 / 0.2e1;
t21 = t145 * t65 + t146 * t67 + t175 * t63;
t22 = t145 * t66 + t146 * t68 + t175 * t64;
t34 = t145 * t86 + t146 * t87 + t175 * t85;
t31 = t34 * t175;
t7 = t21 * t148 - t22 * t150 + t31;
t247 = t7 / 0.2e1;
t246 = t148 / 0.2e1;
t245 = -t150 / 0.2e1;
t244 = t175 / 0.2e1;
t242 = pkin(4) * t151;
t241 = t178 * pkin(2);
t69 = t108 * rSges(6,1) + t107 * rSges(6,2) + t148 * rSges(6,3);
t240 = t149 * pkin(4) + pkin(9) * t148 + t69;
t239 = t254 * t195;
t211 = -rSges(6,1) * t110 - rSges(6,2) * t109;
t70 = -rSges(6,3) * t150 - t211;
t238 = -pkin(9) * t150 + t242 + t70;
t88 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t175;
t237 = pkin(4) * t176 + pkin(9) * t175 + t88;
t121 = -Icges(4,1) * t229 - Icges(4,4) * t178 + Icges(4,5) * t177;
t122 = Icges(3,5) * t178 - Icges(3,6) * t177 - Icges(3,3) * t229;
t228 = -t122 - t121;
t120 = Icges(4,1) * t232 - Icges(4,4) * t180 + Icges(4,5) * t179;
t123 = Icges(3,5) * t180 - Icges(3,6) * t179 + Icges(3,3) * t232;
t227 = t123 + t120;
t167 = t177 * qJ(3);
t137 = t167 + t241;
t138 = t180 * pkin(2) + qJ(3) * t179;
t226 = t137 * t232 + t138 * t229;
t135 = t195 * t138;
t154 = pkin(3) * t232 + pkin(8) * t180;
t225 = t195 * t154 + t135;
t155 = -pkin(3) * t229 + t178 * pkin(8);
t224 = -t137 - t155;
t181 = (pkin(2) * t198 - qJ(3) * t201) * t194;
t223 = -pkin(3) * t195 - pkin(8) * t236 - t181;
t222 = t202 * pkin(1) + pkin(7) * t232;
t221 = t26 / 0.2e1 + t21 / 0.2e1;
t27 = t109 * t86 + t110 * t87 - t150 * t85;
t220 = -t27 / 0.2e1 - t22 / 0.2e1;
t112 = Icges(5,5) * t176 - Icges(5,6) * t175 + Icges(5,3) * t236;
t113 = Icges(5,4) * t176 - Icges(5,2) * t175 + Icges(5,6) * t236;
t114 = Icges(5,1) * t176 - Icges(5,4) * t175 + Icges(5,5) * t236;
t54 = t112 * t236 - t175 * t113 + t176 * t114;
t95 = t149 * rSges(5,1) - t148 * rSges(5,2) + t180 * rSges(5,3);
t131 = t180 * rSges(3,1) - t179 * rSges(3,2) + rSges(3,3) * t232;
t128 = rSges(4,1) * t232 - t180 * rSges(4,2) + t179 * rSges(4,3);
t217 = -t199 * pkin(1) + pkin(7) * t229;
t216 = t194 * (-rSges(4,1) * t195 - (-rSges(4,2) * t198 - rSges(4,3) * t201) * t194 - t181);
t215 = t154 * t229 + t155 * t232 + t226;
t214 = -t167 + t217;
t115 = rSges(5,1) * t176 - rSges(5,2) * t175 + rSges(5,3) * t236;
t213 = t194 * (-t115 + t223);
t212 = -rSges(5,1) * t151 - rSges(5,2) * t150;
t210 = t194 * (t223 - t237);
t209 = t138 + t222;
t89 = Icges(5,5) * t149 - Icges(5,6) * t148 + Icges(5,3) * t180;
t91 = Icges(5,4) * t149 - Icges(5,2) * t148 + Icges(5,6) * t180;
t93 = Icges(5,1) * t149 - Icges(5,4) * t148 + Icges(5,5) * t180;
t40 = -t175 * t91 + t176 * t93 + t236 * t89;
t48 = t112 * t180 - t113 * t148 + t114 * t149;
t208 = t40 / 0.2e1 + t48 / 0.2e1 + t221;
t90 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t178;
t92 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t178;
t94 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t178;
t41 = -t175 * t92 + t176 * t94 + t236 * t90;
t49 = t112 * t178 + t113 * t150 + t114 * t151;
t207 = t49 / 0.2e1 + t41 / 0.2e1 - t220;
t206 = rSges(4,1) * t229 - t177 * rSges(4,3);
t205 = t214 - t155;
t130 = t178 * rSges(3,1) - t177 * rSges(3,2) - rSges(3,3) * t229;
t203 = t154 + t209;
t184 = rSges(2,1) * t202 - t199 * rSges(2,2);
t183 = -t199 * rSges(2,1) - rSges(2,2) * t202;
t164 = rSges(3,3) * t195 + (rSges(3,1) * t198 + rSges(3,2) * t201) * t194;
t129 = -t178 * rSges(4,2) - t206;
t104 = t131 + t222;
t103 = -t130 + t217;
t98 = -t195 * t130 - t164 * t229;
t97 = t131 * t195 - t164 * t232;
t96 = rSges(5,3) * t178 - t212;
t82 = t209 + t128;
t81 = (rSges(4,2) - pkin(2)) * t178 + t206 + t214;
t79 = (t130 * t199 + t131 * t202) * t194;
t78 = t158 * t232 - t159 * t179 + t160 * t180;
t77 = -t158 * t229 - t177 * t159 + t178 * t160;
t76 = t177 * t161 - t178 * t162 - t163 * t229;
t75 = t161 * t179 - t162 * t180 + t163 * t232;
t72 = (-t129 - t137) * t195 + t202 * t216;
t71 = t128 * t195 + t199 * t216 + t135;
t62 = -t115 * t180 + t236 * t95;
t61 = t115 * t178 - t236 * t96;
t60 = t121 * t195 + (-t117 * t201 - t119 * t198) * t194;
t59 = t120 * t195 + (-t116 * t201 - t118 * t198) * t194;
t58 = t123 * t195 + (t125 * t201 + t127 * t198) * t194;
t57 = t122 * t195 + (t124 * t201 + t126 * t198) * t194;
t56 = t203 + t95;
t55 = (-rSges(5,3) - pkin(2)) * t178 + t205 + t212;
t53 = t54 * t195;
t52 = t54 * t236;
t51 = (t128 * t202 + t129 * t199) * t194 + t226;
t50 = -t178 * t95 + t180 * t96;
t47 = (-t96 + t224) * t195 + t202 * t213;
t46 = t195 * t95 + t199 * t213 + t225;
t45 = -t150 * t88 - t175 * t70;
t44 = -t148 * t88 + t175 * t69;
t43 = t203 + t240;
t42 = -t241 - t242 + (rSges(6,3) + pkin(9)) * t150 + t205 + t211;
t39 = t150 * t92 + t151 * t94 + t178 * t90;
t38 = t150 * t91 + t151 * t93 + t178 * t89;
t37 = -t148 * t92 + t149 * t94 + t180 * t90;
t36 = -t148 * t91 + t149 * t93 + t180 * t89;
t35 = (t199 * t96 + t202 * t95) * t194 + t215;
t33 = t34 * t195;
t32 = t34 * t236;
t30 = t148 * t70 + t150 * t69;
t29 = -t180 * t237 + t236 * t240;
t28 = t178 * t237 - t236 * t238;
t25 = (t224 - t238) * t195 + t202 * t210;
t24 = t195 * t240 + t199 * t210 + t225;
t23 = -t178 * t240 + t180 * t238;
t20 = t109 * t66 + t110 * t68 - t150 * t64;
t19 = t109 * t65 + t110 * t67 - t150 * t63;
t16 = (t199 * t238 + t202 * t240) * t194 + t215;
t15 = t53 + (t40 * t199 - t41 * t202) * t194;
t14 = t41 * t178 + t40 * t180 + t52;
t13 = t49 * t195 + (t199 * t38 - t202 * t39) * t194;
t12 = t48 * t195 + (t199 * t36 - t202 * t37) * t194;
t11 = t178 * t39 + t180 * t38 + t236 * t49;
t10 = t178 * t37 + t180 * t36 + t236 * t48;
t9 = t33 + (t21 * t199 - t22 * t202) * t194;
t8 = t22 * t178 + t21 * t180 + t32;
t6 = t27 * t195 + (t19 * t199 - t20 * t202) * t194;
t5 = t26 * t195 + (t17 * t199 - t18 * t202) * t194;
t4 = t178 * t20 + t180 * t19 + t236 * t27;
t3 = t17 * t180 + t178 * t18 + t236 * t26;
t2 = t148 * t19 - t150 * t20 + t175 * t27;
t73 = [Icges(2,3) + m(6) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t55 ^ 2 + t56 ^ 2) + m(4) * (t81 ^ 2 + t82 ^ 2) + m(3) * (t103 ^ 2 + t104 ^ 2) + m(2) * (t183 ^ 2 + t184 ^ 2) + t54 + t34 + t254; t33 + t53 + m(6) * (t24 * t43 + t25 * t42) + m(5) * (t46 * t56 + t47 * t55) + m(4) * (t71 * t82 + t72 * t81) + m(3) * (t103 * t98 + t104 * t97) + ((-t57 / 0.2e1 - t77 / 0.2e1 - t76 / 0.2e1 - t60 / 0.2e1 - t207) * t202 + (t58 / 0.2e1 + t78 / 0.2e1 + t75 / 0.2e1 + t59 / 0.2e1 + t208) * t199) * t194 + t239; (t9 + t15 + t239) * t195 + m(6) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t35 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(4) * (t51 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (t79 ^ 2 + t97 ^ 2 + t98 ^ 2) + ((-t13 - t6 + ((t253 * t177 - t252 * t178) * t194 + t228 * t249 * t202) * t202 + (-t57 - t60 - t76 - t77) * t195) * t202 + (t5 + t12 + ((t251 * t179 + t250 * t180) * t194 + t227 * t249 * t199) * t199 + (t78 + t75 + t58 + t59) * t195 + ((t199 * t228 + t202 * t227) * t194 + t252 * t180 - t253 * t179 - t250 * t178 - t251 * t177) * t229) * t199) * t194; m(6) * (t177 * t43 + t179 * t42) + m(5) * (t177 * t56 + t179 * t55) + m(4) * (t177 * t82 + t179 * t81); m(6) * (-t16 * t235 + t177 * t24 + t179 * t25) + m(5) * (t177 * t46 + t179 * t47 - t235 * t35) + m(4) * (t177 * t71 + t179 * t72 - t235 * t51); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t201 ^ 2 * t249 + t177 ^ 2 + t179 ^ 2); t32 + t52 + m(6) * (t28 * t42 + t29 * t43) + m(5) * (t55 * t61 + t56 * t62) + t208 * t180 + t207 * t178; (t8 / 0.2e1 + t14 / 0.2e1) * t195 + (t5 / 0.2e1 + t12 / 0.2e1) * t180 + (t6 / 0.2e1 + t13 / 0.2e1) * t178 + m(6) * (t16 * t23 + t24 * t29 + t25 * t28) + m(5) * (t35 * t50 + t46 * t62 + t47 * t61) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t202 + (t3 / 0.2e1 + t10 / 0.2e1) * t199 + (t9 / 0.2e1 + t15 / 0.2e1) * t198) * t194; m(5) * (t177 * t62 + t179 * t61 - t235 * t50) + m(6) * (t177 * t29 + t179 * t28 - t23 * t235); (t14 + t8) * t236 + (t3 + t10) * t180 + (t4 + t11) * t178 + m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t50 ^ 2 + t61 ^ 2 + t62 ^ 2); m(6) * (t42 * t45 + t43 * t44) + t31 + t220 * t150 + t221 * t148; m(6) * (t16 * t30 + t24 * t44 + t25 * t45) + t5 * t246 + t6 * t245 + t195 * t247 + t9 * t244 + (t199 * t248 - t202 * t2 / 0.2e1) * t194; m(6) * (t177 * t44 + t179 * t45 - t235 * t30); m(6) * (t23 * t30 + t28 * t45 + t29 * t44) + t8 * t244 + t4 * t245 + t3 * t246 + t180 * t248 + t178 * t2 / 0.2e1 + t236 * t247; m(6) * (t30 ^ 2 + t44 ^ 2 + t45 ^ 2) + t148 * t1 - t150 * t2 + t175 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t73(1), t73(2), t73(4), t73(7), t73(11); t73(2), t73(3), t73(5), t73(8), t73(12); t73(4), t73(5), t73(6), t73(9), t73(13); t73(7), t73(8), t73(9), t73(10), t73(14); t73(11), t73(12), t73(13), t73(14), t73(15);];
Mq = res;
