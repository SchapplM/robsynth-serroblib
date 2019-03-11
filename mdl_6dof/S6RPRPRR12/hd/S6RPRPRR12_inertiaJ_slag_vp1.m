% Calculate joint inertia matrix for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:13
% EndTime: 2019-03-09 04:18:21
% DurationCPUTime: 3.22s
% Computational Cost: add. (4307->373), mult. (7097->538), div. (0->0), fcn. (7501->8), ass. (0->191)
t268 = -Icges(4,4) - Icges(5,6);
t267 = Icges(4,1) + Icges(5,2);
t266 = Icges(4,2) + Icges(5,3);
t175 = cos(qJ(3));
t265 = t268 * t175;
t172 = sin(qJ(3));
t264 = t268 * t172;
t263 = Icges(5,1) + Icges(4,3);
t262 = t266 * t175 - t264;
t261 = t267 * t172 - t265;
t260 = (-Icges(5,5) + Icges(4,6)) * t175 + (-Icges(5,4) + Icges(4,5)) * t172;
t173 = sin(qJ(1));
t259 = -t173 / 0.2e1;
t248 = t173 / 0.2e1;
t176 = cos(qJ(1));
t258 = -t176 / 0.2e1;
t246 = t176 / 0.2e1;
t247 = t175 / 0.2e1;
t242 = rSges(5,2) * t172;
t178 = -t242 + (-rSges(5,3) - qJ(4)) * t175;
t233 = t172 * t176;
t156 = pkin(3) * t233;
t161 = t176 * qJ(2);
t222 = t156 + t161;
t250 = -pkin(1) - pkin(7);
t67 = t178 * t176 + (-rSges(5,1) + t250) * t173 + t222;
t234 = t172 * t173;
t153 = pkin(3) * t234;
t221 = t176 * pkin(1) + t173 * qJ(2);
t214 = t176 * pkin(7) + t221;
t206 = t153 + t214;
t68 = t176 * rSges(5,1) + t178 * t173 + t206;
t257 = m(5) * (t173 * t67 - t176 * t68);
t155 = pkin(8) * t233;
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t226 = t176 * t174;
t125 = -t173 * t171 + t175 * t226;
t227 = t176 * t171;
t126 = t173 * t174 + t175 * t227;
t201 = -t126 * rSges(6,1) - t125 * rSges(6,2);
t236 = qJ(4) * t175;
t49 = t155 + (rSges(6,3) * t172 - t236) * t176 + (-pkin(4) + t250) * t173 + t201 + t222;
t230 = t173 * t175;
t130 = -qJ(4) * t230 + t153;
t133 = t176 * pkin(4) + pkin(8) * t234;
t224 = -t130 - t133;
t127 = -t174 * t230 - t227;
t128 = -t171 * t230 + t226;
t79 = t128 * rSges(6,1) + t127 * rSges(6,2) + rSges(6,3) * t234;
t254 = -t79 + t224;
t50 = t214 - t254;
t256 = m(6) * (t173 * t49 - t176 * t50);
t157 = t174 * pkin(5) + pkin(4);
t177 = -pkin(9) - pkin(8);
t170 = qJ(5) + qJ(6);
t158 = sin(t170);
t159 = cos(t170);
t228 = t176 * t159;
t114 = -t173 * t158 + t175 * t228;
t229 = t176 * t158;
t115 = t173 * t159 + t175 * t229;
t200 = -t115 * rSges(7,1) - t114 * rSges(7,2);
t244 = pkin(5) * t171;
t203 = (-qJ(4) - t244) * t175;
t40 = (-t157 + t250) * t173 + (t203 + (rSges(7,3) - t177) * t172) * t176 + t200 + t222;
t146 = t176 * t157;
t232 = t172 * t177;
t116 = -t159 * t230 - t229;
t117 = -t158 * t230 + t228;
t66 = t117 * rSges(7,1) + t116 * rSges(7,2) + rSges(7,3) * t234;
t41 = t146 + (t203 - t232) * t173 + t206 + t66;
t255 = m(7) * (t173 * t40 - t176 * t41);
t253 = t260 * t173 + t176 * t263;
t252 = (rSges(4,1) * t172 + rSges(4,2) * t175) * t176;
t251 = t173 * t263 - t260 * t176;
t167 = t173 ^ 2;
t169 = t176 ^ 2;
t59 = Icges(7,5) * t115 + Icges(7,6) * t114 - Icges(7,3) * t233;
t61 = Icges(7,4) * t115 + Icges(7,2) * t114 - Icges(7,6) * t233;
t63 = Icges(7,1) * t115 + Icges(7,4) * t114 - Icges(7,5) * t233;
t24 = t175 * t59 + (t158 * t63 + t159 * t61) * t172;
t60 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t234;
t62 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t234;
t64 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t234;
t25 = t175 * t60 + (t158 * t64 + t159 * t62) * t172;
t91 = Icges(7,3) * t175 + (Icges(7,5) * t158 + Icges(7,6) * t159) * t172;
t92 = Icges(7,6) * t175 + (Icges(7,4) * t158 + Icges(7,2) * t159) * t172;
t93 = Icges(7,5) * t175 + (Icges(7,1) * t158 + Icges(7,4) * t159) * t172;
t219 = t175 * t91 + (t158 * t93 + t159 * t92) * t172;
t45 = t219 * t175;
t20 = t116 * t61 + t117 * t63 + t234 * t59;
t21 = t116 * t62 + t117 * t64 + t234 * t60;
t38 = t116 * t92 + t117 * t93 + t234 * t91;
t5 = t38 * t175 + (t173 * t21 - t176 * t20) * t172;
t249 = t5 * t234 + t175 * (t45 + (t173 * t25 - t176 * t24) * t172);
t144 = t175 * rSges(4,1) - t172 * rSges(4,2);
t245 = m(4) * t144;
t65 = -rSges(7,3) * t233 - t200;
t39 = t66 * t233 + t65 * t234;
t179 = t175 * t244 + t232;
t80 = t155 + t179 * t176 + (-pkin(4) + t157) * t173;
t243 = t65 + t80;
t235 = t171 * t172;
t94 = t175 * rSges(7,3) + (rSges(7,1) * t158 + rSges(7,2) * t159) * t172;
t241 = pkin(5) * t235 + (-pkin(8) - t177) * t175 + t94;
t113 = t175 * rSges(6,3) + (rSges(6,1) * t171 + rSges(6,2) * t174) * t172;
t231 = t173 * t113;
t120 = t176 * (t176 * t236 - t156);
t225 = t120 + t176 * (t173 * pkin(4) - t155);
t142 = t175 * pkin(3) + t172 * qJ(4);
t131 = t173 * t142;
t223 = pkin(8) * t230 + t131;
t220 = t167 + t169;
t101 = Icges(6,6) * t175 + (Icges(6,4) * t171 + Icges(6,2) * t174) * t172;
t104 = Icges(6,5) * t175 + (Icges(6,1) * t171 + Icges(6,4) * t174) * t172;
t98 = Icges(6,3) * t175 + (Icges(6,5) * t171 + Icges(6,6) * t174) * t172;
t218 = t172 * t174 * t101 + t104 * t235 + t175 * t98;
t72 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t233;
t74 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t233;
t76 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t233;
t33 = t175 * t72 + (t171 * t76 + t174 * t74) * t172;
t42 = t125 * t101 + t126 * t104 - t233 * t98;
t217 = -t42 / 0.2e1 - t33 / 0.2e1;
t73 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t234;
t75 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t234;
t77 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t234;
t34 = t175 * t73 + (t171 * t77 + t174 * t75) * t172;
t43 = t127 * t101 + t128 * t104 + t234 * t98;
t216 = t43 / 0.2e1 + t34 / 0.2e1;
t215 = rSges(4,1) * t234 + rSges(4,2) * t230 + t176 * rSges(4,3);
t213 = t234 / 0.2e1;
t212 = -t233 / 0.2e1;
t211 = Icges(4,5) * t247 - Icges(5,4) * t175 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(5,5) / 0.2e1) * t172;
t210 = -t175 * pkin(8) - t142;
t209 = t172 * t241;
t208 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t18 = t114 * t61 + t115 * t63 - t233 * t59;
t19 = t114 * t62 + t115 * t64 - t233 * t60;
t11 = t18 * t173 + t19 * t176;
t12 = t20 * t173 + t21 * t176;
t37 = t114 * t92 + t115 * t93 - t233 * t91;
t4 = t37 * t175 + (t173 * t19 - t176 * t18) * t172;
t207 = t11 * t212 + t12 * t213 + t5 * t246 + t4 * t248 + (t24 * t173 + t25 * t176) * t247;
t205 = t45 + (t25 + t38) * t213 + (t24 + t37) * t212;
t204 = -t233 * t4 + t249;
t199 = rSges(5,3) * t175 + t242;
t58 = t175 * t66;
t81 = -t179 * t173 - t133 + t146;
t31 = -t173 * t209 + t175 * t81 + t58;
t32 = -t175 * t243 - t176 * t209;
t198 = t32 * t173 - t31 * t176;
t46 = -t234 * t94 + t58;
t47 = -t175 * t65 - t233 * t94;
t196 = t47 * t173 - t46 * t176;
t52 = -t172 * t231 + t175 * t79;
t78 = -rSges(6,3) * t233 - t201;
t53 = -t113 * t233 - t175 * t78;
t194 = t53 * t173 - t52 * t176;
t54 = t173 * t241 + t223;
t55 = (t210 - t241) * t176;
t193 = t54 * t173 - t55 * t176;
t70 = t223 + t231;
t71 = (-t113 + t210) * t176;
t191 = t70 * t173 - t71 * t176;
t143 = -t175 * rSges(5,2) + t172 * rSges(5,3);
t87 = t173 * t143 + t131;
t88 = (-t142 - t143) * t176;
t190 = t87 * t173 - t88 * t176;
t145 = t176 * rSges(2,1) - t173 * rSges(2,2);
t141 = -t173 * rSges(2,1) - t176 * rSges(2,2);
t119 = -t176 * rSges(3,2) + t173 * rSges(3,3) + t221;
t118 = t176 * rSges(3,3) + t161 + (rSges(3,2) - pkin(1)) * t173;
t83 = t214 + t215;
t82 = t161 + t252 + (-rSges(4,3) + t250) * t173;
t69 = -t173 * t215 + (t173 * rSges(4,3) - t252) * t176;
t51 = t218 * t175;
t48 = t120 + t199 * t169 + (t173 * t199 - t130) * t173;
t44 = (t173 * t78 + t176 * t79) * t172;
t30 = t254 * t173 + t176 * t78 + t225;
t29 = t127 * t75 + t128 * t77 + t234 * t73;
t28 = t127 * t74 + t128 * t76 + t234 * t72;
t27 = t125 * t75 + t126 * t77 - t233 * t73;
t26 = t125 * t74 + t126 * t76 - t233 * t72;
t17 = (t173 * t80 + t176 * t81) * t172 + t39;
t16 = t243 * t176 + (-t66 - t81 + t224) * t173 + t225;
t15 = t28 * t173 + t29 * t176;
t14 = t26 * t173 + t27 * t176;
t8 = t43 * t175 + (t173 * t29 - t176 * t28) * t172;
t7 = t42 * t175 + (t173 * t27 - t176 * t26) * t172;
t1 = [Icges(3,1) + Icges(2,3) + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2) + m(4) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t67 ^ 2 + t68 ^ 2) + m(3) * (t118 ^ 2 + t119 ^ 2) + m(2) * (t141 ^ 2 + t145 ^ 2) + t218 + t219 + (t267 * t175 + t264) * t175 + (t266 * t172 + t265) * t172; t255 + t256 + m(4) * (t173 * t82 - t176 * t83) + t257 + m(3) * (t173 * t118 - t176 * t119); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t208) * t220; m(7) * (t54 * t40 + t55 * t41) + m(6) * (t70 * t49 + t71 * t50) + m(5) * (t87 * t67 + t88 * t68) + (t38 / 0.2e1 + t25 / 0.2e1 - t83 * t245 + t211 * t176 + (Icges(5,4) * t258 + Icges(4,5) * t246 + t248 * t261) * t175 + (Icges(5,5) * t246 + Icges(4,6) * t258 + t259 * t262) * t172 + t216) * t176 + (t37 / 0.2e1 + t24 / 0.2e1 + t82 * t245 + t211 * t173 + (Icges(5,4) * t259 + Icges(4,5) * t248 + t258 * t261) * t175 + (Icges(5,5) * t248 + Icges(4,6) * t259 + t246 * t262) * t172 - t217) * t173; m(5) * t190 + m(6) * t191 + m(7) * t193 + t220 * t245; m(7) * (t16 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t30 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t48 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t144 ^ 2 * t220 + t69 ^ 2) + (t253 * t169 + t12 + t15) * t176 + (t11 + t14 + t251 * t167 + (t253 * t173 + t251 * t176) * t176) * t173; 0.2e1 * (-t255 / 0.2e1 - t256 / 0.2e1 - t257 / 0.2e1) * t175; -0.2e1 * t208 * t220 * t175; m(7) * (t172 * t16 - t175 * t193) + m(6) * (t172 * t30 - t175 * t191) + m(5) * (t172 * t48 - t175 * t190); 0.2e1 * t208 * (t175 ^ 2 * t220 + t172 ^ 2); t51 + m(7) * (t31 * t41 + t32 * t40) + m(6) * (t53 * t49 + t52 * t50) + (t173 * t216 + t176 * t217) * t172 + t205; m(6) * t194 + m(7) * t198; (t33 * t173 + t34 * t176) * t247 + t8 * t246 + t7 * t248 + (t14 * t258 + t15 * t248) * t172 + m(7) * (t17 * t16 + t31 * t55 + t32 * t54) + m(6) * (t44 * t30 + t52 * t71 + t53 * t70) + t207; m(6) * (t44 * t172 - t175 * t194) + m(7) * (t17 * t172 - t175 * t198); t175 * t51 + m(7) * (t17 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t44 ^ 2 + t52 ^ 2 + t53 ^ 2) + ((t175 * t34 + t8) * t173 + (-t175 * t33 - t4 - t7) * t176) * t172 + t249; m(7) * (t47 * t40 + t46 * t41) + t205; m(7) * t196; m(7) * (t39 * t16 + t46 * t55 + t47 * t54) + t207; m(7) * (t39 * t172 - t175 * t196); m(7) * (t39 * t17 + t46 * t31 + t32 * t47) + t204; m(7) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + t204;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
