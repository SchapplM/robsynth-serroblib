% Calculate joint inertia matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:42
% EndTime: 2019-03-09 05:04:48
% DurationCPUTime: 2.63s
% Computational Cost: add. (9969->420), mult. (12735->601), div. (0->0), fcn. (15006->10), ass. (0->201)
t180 = sin(qJ(4));
t181 = sin(qJ(3));
t184 = cos(qJ(4));
t185 = cos(qJ(3));
t144 = -Icges(5,3) * t185 + (Icges(5,5) * t184 - Icges(5,6) * t180) * t181;
t145 = -Icges(6,2) * t185 + (Icges(6,4) * t184 + Icges(6,6) * t180) * t181;
t256 = -t144 - t145;
t178 = qJ(1) + pkin(10);
t175 = sin(t178);
t176 = cos(t178);
t222 = t180 * t185;
t139 = t175 * t222 + t176 * t184;
t220 = t184 * t185;
t140 = t175 * t220 - t176 * t180;
t226 = t175 * t181;
t80 = Icges(6,5) * t140 + Icges(6,6) * t226 + Icges(6,3) * t139;
t84 = Icges(6,4) * t140 + Icges(6,2) * t226 + Icges(6,6) * t139;
t88 = Icges(6,1) * t140 + Icges(6,4) * t226 + Icges(6,5) * t139;
t43 = -t185 * t84 + (t180 * t80 + t184 * t88) * t181;
t82 = Icges(5,5) * t140 - Icges(5,6) * t139 + Icges(5,3) * t226;
t86 = Icges(5,4) * t140 - Icges(5,2) * t139 + Icges(5,6) * t226;
t90 = Icges(5,1) * t140 - Icges(5,4) * t139 + Icges(5,5) * t226;
t45 = -t185 * t82 + (-t180 * t86 + t184 * t90) * t181;
t255 = -t43 - t45;
t141 = -t175 * t184 + t176 * t222;
t142 = t175 * t180 + t176 * t220;
t225 = t176 * t181;
t81 = Icges(6,5) * t142 + Icges(6,6) * t225 + Icges(6,3) * t141;
t85 = Icges(6,4) * t142 + Icges(6,2) * t225 + Icges(6,6) * t141;
t89 = Icges(6,1) * t142 + Icges(6,4) * t225 + Icges(6,5) * t141;
t44 = -t185 * t85 + (t180 * t81 + t184 * t89) * t181;
t83 = Icges(5,5) * t142 - Icges(5,6) * t141 + Icges(5,3) * t225;
t87 = Icges(5,4) * t142 - Icges(5,2) * t141 + Icges(5,6) * t225;
t91 = Icges(5,1) * t142 - Icges(5,4) * t141 + Icges(5,5) * t225;
t46 = -t185 * t83 + (-t180 * t87 + t184 * t91) * t181;
t254 = t44 + t46;
t146 = -Icges(5,6) * t185 + (Icges(5,4) * t184 - Icges(5,2) * t180) * t181;
t223 = t180 * t181;
t143 = -Icges(6,6) * t185 + (Icges(6,5) * t184 + Icges(6,3) * t180) * t181;
t147 = -Icges(6,4) * t185 + (Icges(6,1) * t184 + Icges(6,5) * t180) * t181;
t148 = -Icges(5,5) * t185 + (Icges(5,1) * t184 - Icges(5,4) * t180) * t181;
t221 = t181 * t184;
t252 = t143 * t223 + (t147 + t148) * t221;
t253 = (-t146 * t223 + t185 * t256 + t252) * t185;
t179 = sin(qJ(6));
t183 = cos(qJ(6));
t104 = t139 * t183 - t140 * t179;
t105 = t139 * t179 + t140 * t183;
t59 = Icges(7,5) * t105 + Icges(7,6) * t104 - Icges(7,3) * t226;
t61 = Icges(7,4) * t105 + Icges(7,2) * t104 - Icges(7,6) * t226;
t63 = Icges(7,1) * t105 + Icges(7,4) * t104 - Icges(7,5) * t226;
t16 = t104 * t61 + t105 * t63 - t59 * t226;
t106 = t141 * t183 - t142 * t179;
t107 = t141 * t179 + t142 * t183;
t60 = Icges(7,5) * t107 + Icges(7,6) * t106 - Icges(7,3) * t225;
t62 = Icges(7,4) * t107 + Icges(7,2) * t106 - Icges(7,6) * t225;
t64 = Icges(7,1) * t107 + Icges(7,4) * t106 - Icges(7,5) * t225;
t17 = t104 * t62 + t105 * t64 - t60 * t226;
t153 = (-t179 * t184 + t180 * t183) * t181;
t154 = (t179 * t180 + t183 * t184) * t181;
t108 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t185;
t109 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t185;
t110 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t185;
t28 = t104 * t109 + t105 * t110 - t108 * t226;
t1 = -t28 * t185 + (t16 * t175 + t17 * t176) * t181;
t34 = t139 * t80 + t140 * t88 + t84 * t226;
t35 = t139 * t81 + t140 * t89 + t85 * t226;
t36 = -t139 * t86 + t140 * t90 + t82 * t226;
t37 = -t139 * t87 + t140 * t91 + t83 * t226;
t55 = t139 * t143 + t140 * t147 + t145 * t226;
t56 = -t139 * t146 + t140 * t148 + t144 * t226;
t251 = t1 + (-t55 - t56) * t185 + ((t35 + t37) * t176 + (t34 + t36) * t175) * t181;
t18 = t106 * t61 + t107 * t63 - t59 * t225;
t19 = t106 * t62 + t107 * t64 - t60 * t225;
t29 = t106 * t109 + t107 * t110 - t108 * t225;
t2 = -t29 * t185 + (t175 * t18 + t176 * t19) * t181;
t38 = t141 * t80 + t142 * t88 + t84 * t225;
t39 = t141 * t81 + t142 * t89 + t85 * t225;
t40 = -t141 * t86 + t142 * t90 + t82 * t225;
t41 = -t141 * t87 + t142 * t91 + t83 * t225;
t57 = t141 * t143 + t142 * t147 + t145 * t225;
t58 = -t141 * t146 + t142 * t148 + t144 * t225;
t250 = t2 + (-t58 - t57) * t185 + ((t39 + t41) * t176 + (t38 + t40) * t175) * t181;
t21 = t153 * t61 + t154 * t63 + t185 * t59;
t22 = t153 * t62 + t154 * t64 + t185 * t60;
t212 = t185 * t108 + t153 * t109 + t154 * t110;
t42 = t212 * t185;
t3 = -t42 + (t21 * t175 + t22 * t176) * t181;
t249 = (t175 * t1 + t176 * t2) * t181 - t185 * t3;
t248 = t175 ^ 2;
t247 = t176 ^ 2;
t246 = m(6) + m(7);
t245 = t175 / 0.2e1;
t243 = -t181 / 0.2e1;
t242 = -t185 / 0.2e1;
t241 = -rSges(7,3) - pkin(9);
t161 = t181 * rSges(4,1) + t185 * rSges(4,2);
t240 = m(4) * t161;
t239 = pkin(3) * t185;
t182 = sin(qJ(1));
t238 = t182 * pkin(1);
t236 = t185 * (t22 * t175 - t21 * t176);
t234 = t139 * rSges(6,3);
t233 = t176 * rSges(4,3);
t232 = t107 * rSges(7,1) + t106 * rSges(7,2);
t113 = t142 * pkin(4) + t141 * qJ(5);
t94 = t142 * rSges(6,1) + rSges(6,2) * t225 + t141 * rSges(6,3);
t231 = -t113 - t94;
t197 = -t105 * rSges(7,1) - t104 * rSges(7,2);
t65 = -rSges(7,3) * t226 - t197;
t230 = t140 * pkin(5) - pkin(9) * t226 + t65;
t137 = t142 * pkin(5);
t66 = -rSges(7,3) * t225 + t232;
t229 = -pkin(9) * t225 + t137 + t66;
t228 = Icges(4,4) * t181;
t227 = Icges(4,4) * t185;
t224 = t176 * t185;
t132 = t139 * qJ(5);
t112 = t140 * pkin(4) + t132;
t155 = (pkin(4) * t184 + qJ(5) * t180) * t181;
t219 = t185 * t112 + t155 * t226;
t111 = t154 * rSges(7,1) + t153 * rSges(7,2) + t185 * rSges(7,3);
t218 = pkin(5) * t221 + t185 * pkin(9) + t111;
t213 = pkin(3) * t224 + pkin(8) * t225;
t216 = t248 * (pkin(8) * t181 + t239) + t176 * t213;
t149 = -t185 * rSges(6,2) + (rSges(6,1) * t184 + rSges(6,3) * t180) * t181;
t215 = -t149 - t155;
t150 = -t185 * rSges(5,3) + (rSges(5,1) * t184 - rSges(5,2) * t180) * t181;
t167 = t181 * pkin(3) - t185 * pkin(8);
t214 = -t150 - t167;
t211 = -t22 / 0.2e1 - t29 / 0.2e1;
t210 = -t28 / 0.2e1 - t21 / 0.2e1;
t209 = -t113 - t229;
t208 = -t155 - t218;
t95 = t142 * rSges(5,1) - t141 * rSges(5,2) + rSges(5,3) * t225;
t207 = -t167 + t215;
t186 = cos(qJ(1));
t177 = t186 * pkin(1);
t206 = t176 * pkin(2) + t175 * pkin(7) + t177;
t205 = -pkin(2) - t239;
t204 = t176 * pkin(7) - t238;
t203 = t175 * t112 + t176 * t113 + t216;
t202 = -t167 + t208;
t201 = -t132 + t204;
t199 = rSges(4,1) * t185 - rSges(4,2) * t181;
t198 = -t140 * rSges(5,1) + t139 * rSges(5,2);
t196 = t206 + t213;
t195 = Icges(4,1) * t185 - t228;
t194 = -Icges(4,2) * t181 + t227;
t193 = Icges(4,5) * t185 - Icges(4,6) * t181;
t190 = rSges(4,1) * t224 - rSges(4,2) * t225 + t175 * rSges(4,3);
t189 = t44 / 0.2e1 + t46 / 0.2e1 + t58 / 0.2e1 + t57 / 0.2e1 - t211;
t188 = t43 / 0.2e1 + t45 / 0.2e1 + t56 / 0.2e1 + t55 / 0.2e1 - t210;
t187 = t113 + t196;
t163 = t186 * rSges(2,1) - t182 * rSges(2,2);
t162 = -t182 * rSges(2,1) - t186 * rSges(2,2);
t158 = Icges(4,5) * t181 + Icges(4,6) * t185;
t152 = t176 * rSges(3,1) - t175 * rSges(3,2) + t177;
t151 = -t175 * rSges(3,1) - t176 * rSges(3,2) - t238;
t121 = Icges(4,3) * t175 + t193 * t176;
t120 = -Icges(4,3) * t176 + t193 * t175;
t117 = t214 * t176;
t116 = t214 * t175;
t115 = t190 + t206;
t114 = t233 + (-pkin(2) - t199) * t175 + t204;
t96 = t112 * t225;
t93 = rSges(5,3) * t226 - t198;
t92 = t140 * rSges(6,1) + rSges(6,2) * t226 + t234;
t79 = t207 * t176;
t78 = t207 * t175;
t75 = t176 * t190 + (t199 * t175 - t233) * t175;
t72 = -t150 * t225 - t185 * t95;
t71 = t150 * t226 + t185 * t93;
t70 = t196 + t95;
t69 = ((-rSges(5,3) - pkin(8)) * t181 + t205) * t175 + t198 + t204;
t68 = t202 * t176;
t67 = t202 * t175;
t54 = (-t175 * t95 + t176 * t93) * t181;
t53 = t187 + t94;
t52 = -t234 + (-rSges(6,1) - pkin(4)) * t140 + ((-rSges(6,2) - pkin(8)) * t181 + t205) * t175 + t201;
t51 = t231 * t185 + t215 * t225;
t50 = t149 * t226 + t185 * t92 + t219;
t49 = t111 * t225 + t185 * t66;
t48 = -t111 * t226 - t185 * t65;
t47 = t175 * t93 + t176 * t95 + t216;
t33 = t241 * t225 + t137 + t187 + t232;
t32 = (-pkin(4) - pkin(5)) * t140 + ((-pkin(8) - t241) * t181 + t205) * t175 + t197 + t201;
t31 = t96 + (t231 * t175 + t176 * t92) * t181;
t30 = (t175 * t66 - t176 * t65) * t181;
t25 = t175 * t92 + t176 * t94 + t203;
t24 = t209 * t185 + t208 * t225;
t23 = t230 * t185 + t218 * t226 + t219;
t20 = t96 + (t209 * t175 + t230 * t176) * t181;
t15 = t41 * t175 - t40 * t176;
t14 = t39 * t175 - t38 * t176;
t13 = t37 * t175 - t36 * t176;
t12 = t35 * t175 - t34 * t176;
t11 = t230 * t175 + t229 * t176 + t203;
t5 = t19 * t175 - t18 * t176;
t4 = -t16 * t176 + t17 * t175;
t6 = [Icges(2,3) + Icges(3,3) + (Icges(4,1) * t181 - t180 * t146 + t227) * t181 + (Icges(4,2) * t185 + t228 + t256) * t185 + m(7) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t69 ^ 2 + t70 ^ 2) + m(6) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t114 ^ 2 + t115 ^ 2) + m(3) * (t151 ^ 2 + t152 ^ 2) + m(2) * (t162 ^ 2 + t163 ^ 2) + t212 + t252; 0; m(3) + m(4) + m(5) + t246; m(7) * (t32 * t68 + t33 * t67) + m(5) * (t116 * t70 + t117 * t69) + m(6) * (t52 * t79 + t53 * t78) + ((-Icges(4,6) * t176 + t194 * t175) * t242 + (-Icges(4,5) * t176 + t195 * t175) * t243 - t114 * t240 + t176 * t158 / 0.2e1 - t188) * t176 + (t185 * (Icges(4,6) * t175 + t194 * t176) / 0.2e1 + t181 * (Icges(4,5) * t175 + t195 * t176) / 0.2e1 - t115 * t240 + t158 * t245 + t189) * t175; m(4) * t75 + m(5) * t47 + m(6) * t25 + m(7) * t11; m(7) * (t11 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t25 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2 + t47 ^ 2) + m(4) * (t75 ^ 2 + (t247 + t248) * t161 ^ 2) + (-t247 * t120 - t12 - t13 - t4) * t176 + (t248 * t121 + t14 + t15 + t5 + (-t175 * t120 + t176 * t121) * t176) * t175; -t42 - t253 + m(7) * (t23 * t32 + t24 * t33) + m(5) * (t69 * t71 + t70 * t72) + m(6) * (t50 * t52 + t51 * t53) + (t188 * t175 + t189 * t176) * t181; m(5) * t54 + m(6) * t31 + m(7) * t20; -t236 / 0.2e1 + m(7) * (t11 * t20 + t23 * t68 + t24 * t67) + m(6) * (t25 * t31 + t50 * t79 + t51 * t78) + m(5) * (t116 * t72 + t117 * t71 + t47 * t54) + ((t5 / 0.2e1 + t14 / 0.2e1 + t15 / 0.2e1) * t176 + (t4 / 0.2e1 + t13 / 0.2e1 + t12 / 0.2e1) * t175) * t181 + (t254 * t175 + t255 * t176) * t242 + t250 * t245 - t251 * t176 / 0.2e1; m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t31 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t54 ^ 2 + t71 ^ 2 + t72 ^ 2) + (-t3 + t253) * t185 + ((-t254 * t185 + t250) * t176 + (t255 * t185 + t251) * t175) * t181; m(7) * (t139 * t33 + t141 * t32) + m(6) * (t139 * t53 + t141 * t52); t246 * t223; m(7) * (t11 * t223 + t139 * t67 + t141 * t68) + m(6) * (t139 * t78 + t141 * t79 + t25 * t223); m(7) * (t139 * t24 + t141 * t23 + t20 * t223) + m(6) * (t139 * t51 + t141 * t50 + t31 * t223); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t181 ^ 2 * t180 ^ 2 + t139 ^ 2 + t141 ^ 2); m(7) * (t32 * t48 + t33 * t49) + t42 + (t210 * t175 + t211 * t176) * t181; m(7) * t30; m(7) * (t11 * t30 + t48 * t68 + t49 * t67) + t236 / 0.2e1 + (t5 * t243 + t1 / 0.2e1) * t176 + (-t2 / 0.2e1 + t4 * t243) * t175; m(7) * (t20 * t30 + t23 * t48 + t24 * t49) - t249; m(7) * (t49 * t139 + t48 * t141 + t30 * t223); m(7) * (t30 ^ 2 + t48 ^ 2 + t49 ^ 2) + t249;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
