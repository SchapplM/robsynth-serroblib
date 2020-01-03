% Calculate joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:20
% EndTime: 2019-12-31 21:43:27
% DurationCPUTime: 2.58s
% Computational Cost: add. (9829->478), mult. (25147->687), div. (0->0), fcn. (31672->10), ass. (0->220)
t210 = cos(pkin(5));
t212 = sin(qJ(2));
t216 = cos(qJ(1));
t248 = t216 * t212;
t213 = sin(qJ(1));
t215 = cos(qJ(2));
t250 = t213 * t215;
t195 = t210 * t248 + t250;
t209 = sin(pkin(5));
t258 = sin(qJ(3));
t230 = t209 * t258;
t259 = cos(qJ(3));
t174 = t195 * t259 - t216 * t230;
t247 = t216 * t215;
t251 = t213 * t212;
t197 = -t210 * t251 + t247;
t176 = t197 * t259 + t213 * t230;
t231 = t209 * t259;
t193 = t210 * t258 + t212 * t231;
t175 = t197 * t258 - t213 * t231;
t196 = t210 * t250 + t248;
t211 = sin(qJ(5));
t214 = cos(qJ(5));
t138 = t175 * t214 - t196 * t211;
t139 = t175 * t211 + t196 * t214;
t173 = t195 * t258 + t216 * t231;
t194 = -t210 * t247 + t251;
t136 = t173 * t214 - t194 * t211;
t137 = t173 * t211 + t194 * t214;
t80 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t174;
t82 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t174;
t84 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t174;
t26 = t138 * t82 + t139 * t84 + t176 * t80;
t81 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t176;
t83 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t176;
t85 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t176;
t27 = t138 * t83 + t139 * t85 + t176 * t81;
t192 = -t210 * t259 + t212 * t230;
t253 = t209 * t215;
t171 = t192 * t214 + t211 * t253;
t172 = t192 * t211 - t214 * t253;
t97 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t193;
t98 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t193;
t99 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t193;
t36 = t138 * t98 + t139 * t99 + t176 * t97;
t2 = t26 * t174 + t27 * t176 + t36 * t193;
t263 = t2 / 0.2e1;
t262 = t174 / 0.2e1;
t261 = t176 / 0.2e1;
t260 = t193 / 0.2e1;
t257 = t194 * pkin(4);
t224 = -t137 * rSges(6,1) - t136 * rSges(6,2);
t86 = t174 * rSges(6,3) - t224;
t256 = t174 * pkin(9) + t257 + t86;
t87 = t139 * rSges(6,1) + t138 * rSges(6,2) + t176 * rSges(6,3);
t255 = t196 * pkin(4) + t176 * pkin(9) + t87;
t254 = t209 * t213;
t252 = t209 * t216;
t249 = t213 * t216;
t100 = t172 * rSges(6,1) + t171 * rSges(6,2) + t193 * rSges(6,3);
t246 = -pkin(4) * t253 + t193 * pkin(9) + t100;
t114 = t196 * rSges(5,1) - t176 * rSges(5,2) + t175 * rSges(5,3);
t126 = t176 * pkin(3) + t175 * qJ(4);
t245 = -t114 - t126;
t164 = t173 * qJ(4);
t125 = t174 * pkin(3) + t164;
t161 = t193 * pkin(3) + t192 * qJ(4);
t244 = t125 * t253 + t194 * t161;
t163 = t197 * pkin(2) + t196 * pkin(8);
t160 = t210 * t163;
t243 = t210 * t126 + t160;
t162 = t195 * pkin(2) + t194 * pkin(8);
t242 = -t125 - t162;
t142 = -Icges(5,5) * t253 - Icges(5,6) * t193 + Icges(5,3) * t192;
t143 = -Icges(5,4) * t253 - Icges(5,2) * t193 + Icges(5,6) * t192;
t241 = t192 * t142 - t193 * t143;
t146 = Icges(4,4) * t193 - Icges(4,2) * t192 - Icges(4,6) * t253;
t147 = Icges(4,1) * t193 - Icges(4,4) * t192 - Icges(4,5) * t253;
t240 = -t192 * t146 + t193 * t147;
t148 = -rSges(5,1) * t253 - t193 * rSges(5,2) + t192 * rSges(5,3);
t239 = -t148 - t161;
t238 = t162 * t254 + t163 * t252;
t237 = t216 * pkin(1) + pkin(7) * t254;
t41 = t171 * t98 + t172 * t99 + t193 * t97;
t28 = t171 * t82 + t172 * t84 + t193 * t80;
t35 = t136 * t98 + t137 * t99 + t174 * t97;
t236 = t35 / 0.2e1 + t28 / 0.2e1;
t29 = t171 * t83 + t172 * t85 + t193 * t81;
t235 = t36 / 0.2e1 + t29 / 0.2e1;
t234 = -t126 - t255;
t233 = -t161 - t246;
t116 = t176 * rSges(4,1) - t175 * rSges(4,2) + t196 * rSges(4,3);
t181 = Icges(3,3) * t210 + (Icges(3,5) * t212 + Icges(3,6) * t215) * t209;
t182 = Icges(3,6) * t210 + (Icges(3,4) * t212 + Icges(3,2) * t215) * t209;
t183 = Icges(3,5) * t210 + (Icges(3,1) * t212 + Icges(3,4) * t215) * t209;
t232 = t209 * t212 * t183 + t210 * t181 + t182 * t253;
t157 = t197 * rSges(3,1) - t196 * rSges(3,2) + rSges(3,3) * t254;
t229 = -t213 * pkin(1) + pkin(7) * t252;
t149 = t193 * rSges(4,1) - t192 * rSges(4,2) - rSges(4,3) * t253;
t198 = (pkin(2) * t212 - pkin(8) * t215) * t209;
t228 = t209 * (-t149 - t198);
t227 = t125 * t254 + t126 * t252 + t238;
t226 = t209 * (-t198 + t239);
t225 = -t194 * rSges(5,1) - t173 * rSges(5,3);
t223 = t163 + t237;
t222 = t209 * (-t198 + t233);
t221 = -t162 + t229;
t115 = t174 * rSges(4,1) - t173 * rSges(4,2) + t194 * rSges(4,3);
t220 = -t164 + t221;
t156 = t195 * rSges(3,1) - t194 * rSges(3,2) - rSges(3,3) * t252;
t101 = Icges(5,5) * t194 - Icges(5,6) * t174 + Icges(5,3) * t173;
t105 = Icges(5,4) * t194 - Icges(5,2) * t174 + Icges(5,6) * t173;
t109 = Icges(5,1) * t194 - Icges(5,4) * t174 + Icges(5,5) * t173;
t54 = t192 * t101 - t193 * t105 - t109 * t253;
t103 = Icges(4,5) * t174 - Icges(4,6) * t173 + Icges(4,3) * t194;
t107 = Icges(4,4) * t174 - Icges(4,2) * t173 + Icges(4,6) * t194;
t111 = Icges(4,1) * t174 - Icges(4,4) * t173 + Icges(4,5) * t194;
t56 = -t103 * t253 - t192 * t107 + t193 * t111;
t144 = -Icges(5,1) * t253 - Icges(5,4) * t193 + Icges(5,5) * t192;
t63 = t173 * t142 - t174 * t143 + t194 * t144;
t145 = Icges(4,5) * t193 - Icges(4,6) * t192 - Icges(4,3) * t253;
t65 = t194 * t145 - t173 * t146 + t174 * t147;
t219 = t54 / 0.2e1 + t56 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 + t236;
t102 = Icges(5,5) * t196 - Icges(5,6) * t176 + Icges(5,3) * t175;
t106 = Icges(5,4) * t196 - Icges(5,2) * t176 + Icges(5,6) * t175;
t110 = Icges(5,1) * t196 - Icges(5,4) * t176 + Icges(5,5) * t175;
t55 = t192 * t102 - t193 * t106 - t110 * t253;
t104 = Icges(4,5) * t176 - Icges(4,6) * t175 + Icges(4,3) * t196;
t108 = Icges(4,4) * t176 - Icges(4,2) * t175 + Icges(4,6) * t196;
t112 = Icges(4,1) * t176 - Icges(4,4) * t175 + Icges(4,5) * t196;
t57 = -t104 * t253 - t192 * t108 + t193 * t112;
t64 = t175 * t142 - t176 * t143 + t196 * t144;
t66 = t196 * t145 - t175 * t146 + t176 * t147;
t218 = t57 / 0.2e1 + t66 / 0.2e1 + t64 / 0.2e1 + t55 / 0.2e1 + t235;
t217 = t126 + t223;
t200 = t216 * rSges(2,1) - t213 * rSges(2,2);
t199 = -t213 * rSges(2,1) - t216 * rSges(2,2);
t184 = t210 * rSges(3,3) + (rSges(3,1) * t212 + rSges(3,2) * t215) * t209;
t155 = Icges(3,1) * t197 - Icges(3,4) * t196 + Icges(3,5) * t254;
t154 = Icges(3,1) * t195 - Icges(3,4) * t194 - Icges(3,5) * t252;
t153 = Icges(3,4) * t197 - Icges(3,2) * t196 + Icges(3,6) * t254;
t152 = Icges(3,4) * t195 - Icges(3,2) * t194 - Icges(3,6) * t252;
t151 = Icges(3,5) * t197 - Icges(3,6) * t196 + Icges(3,3) * t254;
t150 = Icges(3,5) * t195 - Icges(3,6) * t194 - Icges(3,3) * t252;
t132 = t157 + t237;
t131 = -t156 + t229;
t120 = -t210 * t156 - t184 * t252;
t119 = t210 * t157 - t184 * t254;
t117 = t196 * t125;
t113 = -t174 * rSges(5,2) - t225;
t96 = t232 * t210;
t94 = (t156 * t213 + t157 * t216) * t209;
t93 = t181 * t254 - t196 * t182 + t197 * t183;
t92 = -t181 * t252 - t194 * t182 + t195 * t183;
t89 = t223 + t116;
t88 = -t115 + t221;
t79 = -t116 * t253 - t196 * t149;
t78 = t115 * t253 + t194 * t149;
t77 = t210 * t151 + (t153 * t215 + t155 * t212) * t209;
t76 = t210 * t150 + (t152 * t215 + t154 * t212) * t209;
t75 = -t145 * t253 + t240;
t74 = -t144 * t253 + t241;
t73 = t75 * t210;
t72 = t74 * t210;
t71 = t217 + t114;
t70 = (rSges(5,2) - pkin(3)) * t174 + t220 + t225;
t69 = t196 * t115 - t194 * t116;
t68 = (-t115 - t162) * t210 + t216 * t228;
t67 = t210 * t116 + t213 * t228 + t160;
t62 = (t115 * t213 + t116 * t216) * t209 + t238;
t61 = -t176 * t100 + t193 * t87;
t60 = t174 * t100 - t193 * t86;
t59 = t239 * t196 + t245 * t253;
t58 = t113 * t253 + t194 * t148 + t244;
t53 = t217 + t255;
t52 = -t257 + (-rSges(6,3) - pkin(3) - pkin(9)) * t174 + t220 + t224;
t51 = (-t113 + t242) * t210 + t216 * t226;
t50 = t210 * t114 + t213 * t226 + t243;
t49 = t196 * t104 - t175 * t108 + t176 * t112;
t48 = t196 * t103 - t175 * t107 + t176 * t111;
t47 = t194 * t104 - t173 * t108 + t174 * t112;
t46 = t194 * t103 - t173 * t107 + t174 * t111;
t45 = t175 * t102 - t176 * t106 + t196 * t110;
t44 = t175 * t101 - t176 * t105 + t196 * t109;
t43 = t173 * t102 - t174 * t106 + t194 * t110;
t42 = t173 * t101 - t174 * t105 + t194 * t109;
t40 = t41 * t210;
t39 = t41 * t193;
t38 = -t174 * t87 + t176 * t86;
t37 = t196 * t113 + t245 * t194 + t117;
t34 = (t113 * t213 + t114 * t216) * t209 + t227;
t33 = t233 * t196 + t234 * t253;
t32 = t246 * t194 + t256 * t253 + t244;
t31 = (t242 - t256) * t210 + t216 * t222;
t30 = t255 * t210 + t213 * t222 + t243;
t25 = t136 * t83 + t137 * t85 + t174 * t81;
t24 = t136 * t82 + t137 * t84 + t174 * t80;
t23 = t234 * t194 + t256 * t196 + t117;
t22 = (t256 * t213 + t255 * t216) * t209 + t227;
t21 = t73 + (t57 * t213 - t56 * t216) * t209;
t20 = t72 + (t55 * t213 - t54 * t216) * t209;
t19 = t56 * t194 + t57 * t196 - t75 * t253;
t18 = t54 * t194 + t55 * t196 - t74 * t253;
t17 = t66 * t210 + (t213 * t49 - t216 * t48) * t209;
t16 = t65 * t210 + (t213 * t47 - t216 * t46) * t209;
t15 = t64 * t210 + (t213 * t45 - t216 * t44) * t209;
t14 = t63 * t210 + (t213 * t43 - t216 * t42) * t209;
t13 = t48 * t194 + t49 * t196 - t66 * t253;
t12 = t46 * t194 + t47 * t196 - t65 * t253;
t11 = t44 * t194 + t45 * t196 - t64 * t253;
t10 = t42 * t194 + t43 * t196 - t63 * t253;
t9 = t40 + (t29 * t213 - t28 * t216) * t209;
t8 = t28 * t194 + t29 * t196 - t41 * t253;
t7 = t28 * t174 + t29 * t176 + t39;
t6 = t36 * t210 + (t213 * t27 - t216 * t26) * t209;
t5 = t35 * t210 + (t213 * t25 - t216 * t24) * t209;
t4 = t26 * t194 + t27 * t196 - t36 * t253;
t3 = t24 * t194 + t25 * t196 - t35 * t253;
t1 = t24 * t174 + t25 * t176 + t35 * t193;
t90 = [Icges(2,3) + (-t144 - t145) * t253 + m(6) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t70 ^ 2 + t71 ^ 2) + m(4) * (t88 ^ 2 + t89 ^ 2) + m(3) * (t131 ^ 2 + t132 ^ 2) + m(2) * (t199 ^ 2 + t200 ^ 2) + t232 + t41 + t240 + t241; t40 + t72 + t73 + t96 + m(6) * (t30 * t53 + t31 * t52) + m(5) * (t50 * t71 + t51 * t70) + m(4) * (t67 * t89 + t68 * t88) + m(3) * (t119 * t132 + t120 * t131) + ((-t76 / 0.2e1 - t92 / 0.2e1 - t219) * t216 + (t77 / 0.2e1 + t93 / 0.2e1 + t218) * t213) * t209; (t21 + t20 + t96 + t9) * t210 + m(5) * (t34 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t62 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t119 ^ 2 + t120 ^ 2 + t94 ^ 2) + m(6) * (t22 ^ 2 + t30 ^ 2 + t31 ^ 2) + (t213 * t6 - t216 * t5 - t216 * t16 - t216 * t14 + t213 * t17 + t213 * t15 + (t213 * ((-t196 * t153 + t197 * t155) * t213 - (-t196 * t152 + t197 * t154) * t216) - t216 * ((-t194 * t153 + t195 * t155) * t213 - (-t194 * t152 + t195 * t154) * t216) + (t213 * (t151 * t213 ^ 2 - t150 * t249) - t216 * (t150 * t216 ^ 2 - t151 * t249)) * t209) * t209 + ((-t76 - t92) * t216 + (t77 + t93) * t213) * t210) * t209; (-t41 - t74 - t75) * t253 + m(6) * (t32 * t52 + t33 * t53) + m(5) * (t58 * t70 + t59 * t71) + m(4) * (t78 * t88 + t79 * t89) + t218 * t196 + t219 * t194; (t8 / 0.2e1 + t19 / 0.2e1 + t18 / 0.2e1) * t210 + (t6 / 0.2e1 + t17 / 0.2e1 + t15 / 0.2e1) * t196 + (t5 / 0.2e1 + t16 / 0.2e1 + t14 / 0.2e1) * t194 + m(6) * (t23 * t22 + t33 * t30 + t32 * t31) + m(5) * (t37 * t34 + t59 * t50 + t58 * t51) + m(4) * (t69 * t62 + t79 * t67 + t78 * t68) + ((-t3 / 0.2e1 - t12 / 0.2e1 - t10 / 0.2e1) * t216 + (-t9 / 0.2e1 - t21 / 0.2e1 - t20 / 0.2e1) * t215 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t213) * t209; (-t18 - t19 - t8) * t253 + (t4 + t11 + t13) * t196 + (t3 + t12 + t10) * t194 + m(6) * (t23 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t37 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t69 ^ 2 + t78 ^ 2 + t79 ^ 2); m(6) * (t173 * t53 + t175 * t52) + m(5) * (t173 * t71 + t175 * t70); m(6) * (t173 * t30 + t175 * t31 + t192 * t22) + m(5) * (t173 * t50 + t175 * t51 + t192 * t34); m(6) * (t173 * t33 + t175 * t32 + t192 * t23) + m(5) * (t173 * t59 + t175 * t58 + t192 * t37); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t173 ^ 2 + t175 ^ 2 + t192 ^ 2); m(6) * (t60 * t52 + t61 * t53) + t39 + t235 * t176 + t236 * t174; m(6) * (t38 * t22 + t61 * t30 + t60 * t31) + t6 * t261 + t5 * t262 + t9 * t260 + t210 * t7 / 0.2e1 + (t213 * t263 - t216 * t1 / 0.2e1) * t209; m(6) * (t38 * t23 + t60 * t32 + t61 * t33) + t194 * t1 / 0.2e1 + t8 * t260 + t3 * t262 + t196 * t263 + t4 * t261 - t7 * t253 / 0.2e1; m(6) * (t61 * t173 + t60 * t175 + t38 * t192); m(6) * (t38 ^ 2 + t60 ^ 2 + t61 ^ 2) + t176 * t2 + t174 * t1 + t193 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t90(1), t90(2), t90(4), t90(7), t90(11); t90(2), t90(3), t90(5), t90(8), t90(12); t90(4), t90(5), t90(6), t90(9), t90(13); t90(7), t90(8), t90(9), t90(10), t90(14); t90(11), t90(12), t90(13), t90(14), t90(15);];
Mq = res;
