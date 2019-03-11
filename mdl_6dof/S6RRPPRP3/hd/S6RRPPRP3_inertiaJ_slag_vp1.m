% Calculate joint inertia matrix for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:36
% EndTime: 2019-03-09 08:33:41
% DurationCPUTime: 3.18s
% Computational Cost: add. (2694->377), mult. (6060->534), div. (0->0), fcn. (6255->6), ass. (0->174)
t273 = Icges(3,1) + Icges(4,1);
t272 = Icges(5,1) + Icges(4,3);
t268 = Icges(5,5) + Icges(3,6);
t264 = Icges(5,6) + Icges(4,4) + Icges(3,5);
t165 = cos(qJ(2));
t271 = (Icges(5,4) - Icges(4,5)) * t165;
t162 = sin(qJ(2));
t270 = (Icges(3,4) - Icges(4,5)) * t162;
t160 = -qJ(6) - pkin(8);
t269 = rSges(7,3) - t160;
t267 = t272 * t162 - t271;
t266 = t273 * t165 - t270;
t265 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t263 = -t264 * t165 + (-Icges(4,6) + t268) * t162;
t163 = sin(qJ(1));
t262 = -t163 / 0.2e1;
t261 = t163 / 0.2e1;
t166 = cos(qJ(1));
t260 = -t166 / 0.2e1;
t259 = t166 / 0.2e1;
t161 = sin(qJ(5));
t164 = cos(qJ(5));
t79 = Icges(7,3) * t162 + (-Icges(7,5) * t164 + Icges(7,6) * t161) * t165;
t80 = Icges(6,3) * t162 + (-Icges(6,5) * t164 + Icges(6,6) * t161) * t165;
t87 = Icges(7,6) * t162 + (-Icges(7,4) * t164 + Icges(7,2) * t161) * t165;
t88 = Icges(6,6) * t162 + (-Icges(6,4) * t164 + Icges(6,2) * t161) * t165;
t258 = (t87 + t88) * t161 * t165 + (t79 + t80) * t162;
t95 = Icges(7,5) * t162 + (-Icges(7,1) * t164 + Icges(7,4) * t161) * t165;
t96 = Icges(6,5) * t162 + (-Icges(6,1) * t164 + Icges(6,4) * t161) * t165;
t257 = (-t95 - t96) * t164;
t213 = t164 * t166;
t217 = t163 * t161;
t116 = -t162 * t217 + t213;
t216 = t163 * t164;
t220 = t161 * t166;
t117 = t162 * t216 + t220;
t215 = t163 * t165;
t48 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t215;
t52 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t215;
t56 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t215;
t12 = t116 * t52 + t117 * t56 + t215 * t48;
t219 = t162 * t166;
t118 = -t161 * t219 - t216;
t119 = t162 * t213 - t217;
t212 = t165 * t166;
t49 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t212;
t53 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t212;
t57 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t212;
t13 = t116 * t53 + t117 * t57 + t215 * t49;
t50 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t215;
t54 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t215;
t58 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t215;
t14 = t116 * t54 + t117 * t58 + t215 * t50;
t51 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t212;
t55 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t212;
t59 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t212;
t15 = t116 * t55 + t117 * t59 + t215 * t51;
t29 = t116 * t87 + t117 * t95 + t215 * t79;
t30 = t116 * t88 + t117 * t96 + t215 * t80;
t256 = ((t13 + t15) * t166 + (t12 + t14) * t163) * t165 + (t29 + t30) * t162;
t16 = t118 * t52 + t119 * t56 + t212 * t48;
t17 = t118 * t53 + t119 * t57 + t212 * t49;
t18 = t118 * t54 + t119 * t58 + t212 * t50;
t19 = t118 * t55 + t119 * t59 + t212 * t51;
t31 = t118 * t87 + t119 * t95 + t212 * t79;
t32 = t118 * t88 + t119 * t96 + t212 * t80;
t255 = ((t17 + t19) * t166 + (t16 + t18) * t163) * t165 + (t31 + t32) * t162;
t254 = t162 / 0.2e1;
t22 = t162 * t48 + (t161 * t52 - t164 * t56) * t165;
t24 = t162 * t50 + (t161 * t54 - t164 * t58) * t165;
t252 = t22 + t24;
t23 = t162 * t49 + (t161 * t53 - t164 * t57) * t165;
t25 = t162 * t51 + (t161 * t55 - t164 * t59) * t165;
t251 = t23 + t25;
t150 = pkin(5) * t164 + pkin(4);
t250 = t119 * rSges(7,1) + t118 * rSges(7,2) + t150 * t219 + t269 * t212;
t249 = t263 * t163 + t265 * t166;
t248 = t265 * t163 - t263 * t166;
t157 = t163 ^ 2;
t159 = t166 ^ 2;
t247 = m(4) / 0.2e1;
t246 = m(5) / 0.2e1;
t245 = m(6) / 0.2e1;
t244 = -pkin(2) - pkin(3);
t133 = rSges(3,1) * t162 + rSges(3,2) * t165;
t240 = m(3) * t133;
t239 = pkin(4) - t150;
t238 = -pkin(8) - t160;
t237 = (t165 * t257 + t258) * t162;
t192 = -t117 * rSges(7,1) - t116 * rSges(7,2);
t236 = pkin(5) * t220 + (-t162 * t239 + t165 * t238) * t163 + rSges(7,3) * t215 - t192;
t209 = pkin(4) * t219 + pkin(8) * t212;
t235 = -pkin(5) * t217 - t209 + t250;
t232 = t166 * rSges(4,2);
t231 = t166 * rSges(3,3);
t230 = -rSges(5,3) - qJ(4);
t229 = (-rSges(7,1) * t164 + rSges(7,2) * t161 + t239) * t165 + (t238 + rSges(7,3)) * t162;
t227 = Icges(3,4) * t165;
t226 = Icges(5,4) * t162;
t222 = qJ(4) * t166;
t218 = t163 * qJ(4);
t208 = pkin(2) * t212 + qJ(3) * t219;
t211 = t157 * (pkin(2) * t165 + qJ(3) * t162) + t166 * t208;
t131 = pkin(2) * t162 - qJ(3) * t165;
t210 = -rSges(4,1) * t162 + rSges(4,3) * t165 - t131;
t207 = t166 * pkin(1) + t163 * pkin(7);
t206 = t157 + t159;
t65 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t212;
t205 = rSges(4,1) * t212 + t163 * rSges(4,2) + rSges(4,3) * t219;
t204 = -pkin(3) * t162 - t131;
t203 = -pkin(5) * t161 - qJ(4);
t202 = t246 + t245 + m(7) / 0.2e1;
t148 = pkin(3) * t212;
t201 = t163 * (pkin(3) * t215 + t222) + t166 * (t148 - t218) + t211;
t200 = t207 + t208;
t199 = rSges(5,1) * t165 + rSges(5,2) * t162 + t204;
t198 = pkin(4) * t165 - pkin(8) * t162 + t204;
t197 = rSges(5,1) * t219 - rSges(5,2) * t212;
t195 = t264 * t254 + (-Icges(4,6) / 0.2e1 + t268 / 0.2e1) * t165;
t194 = rSges(3,1) * t165 - rSges(3,2) * t162;
t193 = -t117 * rSges(6,1) - t116 * rSges(6,2);
t20 = -t162 * t236 + t215 * t229;
t21 = t162 * t235 - t212 * t229;
t188 = t163 * t21 + t166 * t20;
t167 = t198 - t229;
t36 = t167 * t163;
t37 = t167 * t166;
t187 = t163 * t36 + t166 * t37;
t186 = t148 + t200;
t179 = -Icges(3,2) * t162 + t227;
t177 = -Icges(5,2) * t165 + t226;
t104 = rSges(6,3) * t162 + (-rSges(6,1) * t164 + rSges(6,2) * t161) * t165;
t173 = -t104 + t198;
t172 = rSges(3,1) * t212 - rSges(3,2) * t219 + t163 * rSges(3,3);
t171 = t24 / 0.2e1 + t22 / 0.2e1 + t30 / 0.2e1 + t29 / 0.2e1;
t170 = t25 / 0.2e1 + t23 / 0.2e1 + t32 / 0.2e1 + t31 / 0.2e1;
t169 = t157 * (pkin(4) * t162 + pkin(8) * t165) + t166 * t209 + t201;
t154 = t166 * pkin(7);
t27 = t154 + t203 * t166 + (-pkin(1) + (-qJ(3) - t150) * t162 + (t244 - t269) * t165) * t163 + t192;
t28 = t163 * t203 + t186 + t250;
t168 = m(7) * (t163 * t28 + t166 * t27);
t158 = t165 ^ 2;
t156 = t162 ^ 2;
t136 = rSges(2,1) * t166 - t163 * rSges(2,2);
t134 = -t163 * rSges(2,1) - rSges(2,2) * t166;
t73 = t210 * t166;
t72 = t210 * t163;
t71 = t172 + t207;
t70 = t231 + t154 + (-pkin(1) - t194) * t163;
t69 = t199 * t166;
t68 = t199 * t163;
t67 = t200 + t205;
t66 = t232 + t154 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t165 + (-rSges(4,3) - qJ(3)) * t162) * t163;
t63 = rSges(6,3) * t215 - t193;
t47 = t166 * t172 + (t163 * t194 - t231) * t163;
t46 = t163 * t230 + t186 + t197;
t45 = t154 + t230 * t166 + (-pkin(1) + (-rSges(5,1) - qJ(3)) * t162 + (rSges(5,2) + t244) * t165) * t163;
t44 = t173 * t166;
t43 = t173 * t163;
t42 = -t104 * t212 + t162 * t65;
t41 = t104 * t215 - t162 * t63;
t38 = t166 * t205 + (-t232 + (rSges(4,1) * t165 + rSges(4,3) * t162) * t163) * t163 + t211;
t35 = t186 + t209 + t65 - t218;
t34 = -t222 + t154 + (-pkin(1) + (-pkin(4) - qJ(3)) * t162 + (-rSges(6,3) - pkin(8) + t244) * t165) * t163 + t193;
t33 = (-t163 * t65 + t166 * t63) * t165;
t26 = t166 * t197 + (rSges(5,1) * t162 - rSges(5,2) * t165) * t157 + t201;
t11 = (-t163 * t235 + t166 * t236) * t165;
t10 = t163 * t63 + t166 * t65 + t169;
t9 = t163 * t236 + t166 * t235 + t169;
t8 = t19 * t163 - t166 * t18;
t7 = -t16 * t166 + t17 * t163;
t6 = -t14 * t166 + t15 * t163;
t5 = -t12 * t166 + t13 * t163;
t1 = [Icges(2,3) + m(7) * (t27 ^ 2 + t28 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2) + m(5) * (t45 ^ 2 + t46 ^ 2) + m(4) * (t66 ^ 2 + t67 ^ 2) + m(3) * (t70 ^ 2 + t71 ^ 2) + m(2) * (t134 ^ 2 + t136 ^ 2) + (t226 + t257 + (Icges(3,2) + t272) * t165 + t270) * t165 + (t227 + (Icges(5,2) + t273) * t162 + t271) * t162 + t258; m(5) * (t45 * t69 + t46 * t68) + m(4) * (t66 * t73 + t67 * t72) + m(6) * (t34 * t44 + t35 * t43) + m(7) * (t27 * t37 + t28 * t36) + (-t70 * t240 + t195 * t166 + (Icges(4,6) * t260 + t179 * t262 + t268 * t259 + t267 * t261) * t165 - t171) * t166 + (-t71 * t240 + t195 * t163 + (Icges(4,6) * t262 + t179 * t259 + t267 * t260 + t268 * t261) * t165 + t170) * t163 + (t177 * t260 * t163 + t266 * t166 * t262 + (t264 * t163 + t166 * t177) * t261 + (t266 * t163 + t264 * t166) * t259) * t162; m(7) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t26 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t38 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(3) * (t133 ^ 2 * t206 + t47 ^ 2) + (t249 * t159 - t5 - t6) * t166 + (t7 + t8 + t248 * t157 + (t249 * t163 + t248 * t166) * t166) * t163; 0.2e1 * (t168 / 0.2e1 + (t163 * t35 + t166 * t34) * t245 + (t163 * t46 + t166 * t45) * t246 + (t163 * t67 + t166 * t66) * t247) * t162; m(7) * (t162 * t187 - t165 * t9) + m(6) * (-t165 * t10 + (t163 * t43 + t166 * t44) * t162) + m(5) * (-t165 * t26 + (t163 * t68 + t166 * t69) * t162) + m(4) * (-t165 * t38 + (t163 * t72 + t166 * t73) * t162); 0.2e1 * (t247 + t202) * (t156 * t206 + t158); m(7) * (-t163 * t27 + t166 * t28) + m(6) * (-t163 * t34 + t166 * t35) + m(5) * (-t163 * t45 + t166 * t46); m(7) * (-t163 * t37 + t166 * t36) + m(6) * (-t163 * t44 + t166 * t43) + m(5) * (-t163 * t69 + t166 * t68); 0; 0.2e1 * t202 * t206; m(7) * (t20 * t27 + t21 * t28) + m(6) * (t34 * t41 + t35 * t42) + (t163 * t171 + t166 * t170) * t165 + t237; m(7) * (t11 * t9 + t20 * t37 + t21 * t36) + m(6) * (t10 * t33 + t41 * t44 + t42 * t43) + ((t8 / 0.2e1 + t7 / 0.2e1) * t166 + (t6 / 0.2e1 + t5 / 0.2e1) * t163) * t165 + (t251 * t163 - t252 * t166) * t254 + t255 * t261 + t256 * t260; m(6) * (-t33 * t165 + (t163 * t42 + t166 * t41) * t162) + m(7) * (-t11 * t165 + t162 * t188); m(6) * (-t41 * t163 + t166 * t42) + m(7) * (-t20 * t163 + t166 * t21); t237 * t162 + m(7) * (t11 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t33 ^ 2 + t41 ^ 2 + t42 ^ 2) + (t255 * t166 + t256 * t163 + (t252 * t163 + t251 * t166) * t162) * t165; t165 * t168; m(7) * (t162 * t9 + t165 * t187); m(7) * (-0.1e1 + t206) * t165 * t162; 0; m(7) * (t162 * t11 + t165 * t188); m(7) * (t158 * t206 + t156);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
