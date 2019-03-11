% Calculate joint inertia matrix for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:46
% EndTime: 2019-03-09 03:57:52
% DurationCPUTime: 2.48s
% Computational Cost: add. (6515->365), mult. (6805->523), div. (0->0), fcn. (7216->10), ass. (0->181)
t258 = Icges(4,3) + Icges(5,3);
t169 = qJ(3) + pkin(10);
t160 = sin(t169);
t161 = cos(t169);
t175 = sin(qJ(3));
t178 = cos(qJ(3));
t257 = Icges(4,5) * t175 + Icges(5,5) * t160 + Icges(4,6) * t178 + Icges(5,6) * t161;
t256 = Icges(4,5) * t178 + Icges(5,5) * t161 - Icges(4,6) * t175 - Icges(5,6) * t160;
t176 = sin(qJ(1));
t179 = cos(qJ(1));
t255 = t257 * t176 + t258 * t179;
t254 = (rSges(4,1) * t175 + rSges(4,2) * t178) * t179;
t253 = t258 * t176 - t257 * t179;
t170 = t176 ^ 2;
t171 = t179 ^ 2;
t227 = t161 * t179;
t172 = qJ(5) + qJ(6);
t162 = sin(t172);
t163 = cos(t172);
t219 = t179 * t163;
t225 = t176 * t162;
t108 = -t160 * t225 + t219;
t220 = t179 * t162;
t224 = t176 * t163;
t109 = t160 * t224 + t220;
t228 = t161 * t176;
t57 = Icges(7,5) * t109 + Icges(7,6) * t108 - Icges(7,3) * t228;
t59 = Icges(7,4) * t109 + Icges(7,2) * t108 - Icges(7,6) * t228;
t61 = Icges(7,1) * t109 + Icges(7,4) * t108 - Icges(7,5) * t228;
t26 = t160 * t57 + (-t162 * t59 + t163 * t61) * t161;
t110 = t160 * t220 + t224;
t111 = -t160 * t219 + t225;
t58 = Icges(7,5) * t111 + Icges(7,6) * t110 + Icges(7,3) * t227;
t60 = Icges(7,4) * t111 + Icges(7,2) * t110 + Icges(7,6) * t227;
t62 = Icges(7,1) * t111 + Icges(7,4) * t110 + Icges(7,5) * t227;
t27 = t160 * t58 + (-t162 * t60 + t163 * t62) * t161;
t90 = Icges(7,6) * t160 + (Icges(7,4) * t163 - Icges(7,2) * t162) * t161;
t236 = t162 * t90;
t89 = Icges(7,3) * t160 + (Icges(7,5) * t163 - Icges(7,6) * t162) * t161;
t91 = Icges(7,5) * t160 + (Icges(7,1) * t163 - Icges(7,4) * t162) * t161;
t239 = t161 * t163 * t91 + t160 * t89;
t42 = (-t161 * t236 + t239) * t160;
t20 = t110 * t59 + t111 * t61 + t227 * t57;
t21 = t110 * t60 + t111 * t62 + t227 * t58;
t38 = t110 * t90 + t111 * t91 + t227 * t89;
t5 = t38 * t160 + (-t176 * t20 + t179 * t21) * t161;
t252 = t5 * t227 + t160 * (t42 + (-t176 * t26 + t179 * t27) * t161);
t250 = t160 / 0.2e1;
t247 = t176 / 0.2e1;
t245 = t179 / 0.2e1;
t145 = t178 * rSges(4,1) - t175 * rSges(4,2);
t244 = m(4) * t145;
t243 = pkin(3) * t178;
t180 = -pkin(9) - pkin(8);
t242 = -pkin(8) - t180;
t63 = t109 * rSges(7,1) + t108 * rSges(7,2) - rSges(7,3) * t228;
t229 = t160 * t176;
t149 = pkin(4) * t229;
t121 = -pkin(8) * t228 + t149;
t177 = cos(qJ(5));
t159 = t177 * pkin(5) + pkin(4);
t174 = sin(qJ(5));
t218 = t179 * t174;
t208 = pkin(5) * t218 + t159 * t229 + t180 * t228;
t67 = -t121 + t208;
t241 = -t63 - t67;
t150 = t179 * t160 * pkin(4);
t223 = t176 * t174;
t230 = t159 * t160;
t192 = -t111 * rSges(7,1) - t110 * rSges(7,2);
t64 = rSges(7,3) * t227 - t192;
t240 = t64 + pkin(5) * t223 + t150 + (t161 * t242 - t230) * t179;
t92 = t160 * rSges(7,3) + (rSges(7,1) * t163 - rSges(7,2) * t162) * t161;
t47 = t160 * t63 + t92 * t228;
t93 = Icges(6,3) * t160 + (Icges(6,5) * t177 - Icges(6,6) * t174) * t161;
t95 = Icges(6,5) * t160 + (Icges(6,1) * t177 - Icges(6,4) * t174) * t161;
t238 = t161 * t177 * t95 + t160 * t93;
t85 = (-pkin(4) + t159) * t161 + t242 * t160;
t237 = t85 + t92;
t94 = Icges(6,6) * t160 + (Icges(6,4) * t177 - Icges(6,2) * t174) * t161;
t235 = t174 * t94;
t166 = t179 * rSges(5,3);
t226 = t175 * t176;
t222 = t176 * t177;
t221 = t176 * t178;
t217 = t179 * t177;
t173 = -qJ(4) - pkin(7);
t212 = t179 * t175 * pkin(3) + t176 * t173;
t112 = t179 * (-t176 * pkin(7) - t212);
t216 = t179 * (pkin(8) * t227 - t150) + t112;
t125 = -t160 * t223 + t217;
t126 = t160 * t222 + t218;
t215 = t126 * rSges(6,1) + t125 * rSges(6,2);
t154 = pkin(3) * t226;
t124 = t154 + (-pkin(7) - t173) * t179;
t214 = -t121 - t124;
t135 = t161 * pkin(4) + t160 * pkin(8);
t155 = pkin(3) * t221;
t213 = t176 * t135 + t155;
t211 = t179 * pkin(1) + t176 * qJ(2);
t151 = t170 + t171;
t69 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t228;
t71 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t228;
t73 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t228;
t33 = t160 * t69 + (-t174 * t71 + t177 * t73) * t161;
t40 = t125 * t94 + t126 * t95 - t228 * t93;
t210 = -t33 / 0.2e1 - t40 / 0.2e1;
t127 = t160 * t218 + t222;
t128 = -t160 * t217 + t223;
t70 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t227;
t72 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t227;
t74 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t227;
t34 = t160 * t70 + (-t174 * t72 + t177 * t74) * t161;
t41 = t127 * t94 + t128 * t95 + t227 * t93;
t209 = t34 / 0.2e1 + t41 / 0.2e1;
t207 = (m(5) + m(6) + m(7)) * t151;
t206 = -rSges(5,1) * t229 - rSges(5,2) * t228 - t166;
t205 = rSges(4,1) * t226 + rSges(4,2) * t221 + t179 * rSges(4,3);
t165 = t179 * qJ(2);
t204 = t165 + t212;
t203 = t161 * (-rSges(6,3) - pkin(8));
t202 = -t228 / 0.2e1;
t201 = t227 / 0.2e1;
t199 = -t135 - t243;
t18 = t108 * t59 + t109 * t61 - t228 * t57;
t19 = t108 * t60 + t109 * t62 - t228 * t58;
t11 = t19 * t176 + t18 * t179;
t12 = t21 * t176 + t20 * t179;
t37 = t108 * t90 + t109 * t91 - t228 * t89;
t4 = t37 * t160 + (-t176 * t18 + t179 * t19) * t161;
t198 = t11 * t202 + t12 * t201 + t4 * t245 + t5 * t247 + (t27 * t176 + t26 * t179) * t250;
t197 = t42 + (t26 + t37) * t202 + (t27 + t38) * t201;
t196 = -t228 * t4 + t252;
t194 = rSges(5,1) * t160 + rSges(5,2) * t161;
t193 = -t128 * rSges(6,1) - t127 * rSges(6,2);
t181 = -t179 * t173 + t154 + t211;
t146 = t179 * rSges(2,1) - t176 * rSges(2,2);
t144 = -t176 * rSges(2,1) - t179 * rSges(2,2);
t133 = t161 * rSges(5,1) - t160 * rSges(5,2);
t123 = -t179 * rSges(3,2) + t176 * rSges(3,3) + t211;
t122 = t179 * rSges(3,3) + t165 + (rSges(3,2) - pkin(1)) * t176;
t98 = (-t133 - t243) * t179;
t97 = t176 * t133 + t155;
t96 = t160 * rSges(6,3) + (rSges(6,1) * t177 - rSges(6,2) * t174) * t161;
t87 = t179 * pkin(7) + t205 + t211;
t86 = t165 + t254 + (-rSges(4,3) - pkin(1) - pkin(7)) * t176;
t82 = t92 * t227;
t79 = t181 - t206;
t78 = t194 * t179 + (-rSges(5,3) - pkin(1)) * t176 + t204;
t77 = -t176 * t205 + (t176 * rSges(4,3) - t254) * t179;
t76 = rSges(6,3) * t227 - t193;
t75 = -rSges(6,3) * t228 + t215;
t66 = (t199 - t96) * t179;
t65 = t176 * t96 + t213;
t55 = t176 * t203 + t149 + t181 + t215;
t54 = -t176 * pkin(1) + t179 * t203 + t150 + t193 + t204;
t53 = -t160 * t76 + t227 * t96;
t52 = t160 * t75 + t228 * t96;
t51 = (t199 - t237) * t179;
t50 = t176 * t237 + t213;
t49 = t112 - t194 * t171 + (-t124 + t206 + t166) * t176;
t48 = -t160 * t64 + t82;
t46 = (-t161 * t235 + t238) * t160;
t45 = t181 + t63 + t208;
t44 = (-pkin(5) * t174 - pkin(1)) * t176 + (t230 + (-rSges(7,3) + t180) * t161) * t179 + t192 + t204;
t43 = (-t176 * t76 - t179 * t75) * t161;
t39 = (-t176 * t64 - t179 * t63) * t161;
t32 = t179 * t76 + (-t75 + t214) * t176 + t216;
t31 = t127 * t72 + t128 * t74 + t227 * t70;
t30 = t127 * t71 + t128 * t73 + t227 * t69;
t29 = t125 * t72 + t126 * t74 - t228 * t70;
t28 = t125 * t71 + t126 * t73 - t228 * t69;
t25 = -t160 * t240 + t227 * t85 + t82;
t24 = t160 * t67 + t228 * t85 + t47;
t17 = (-t176 * t240 + t179 * t241) * t161;
t16 = t240 * t179 + (t214 + t241) * t176 + t216;
t15 = t31 * t176 + t30 * t179;
t14 = t29 * t176 + t28 * t179;
t8 = t41 * t160 + (-t176 * t30 + t179 * t31) * t161;
t7 = t40 * t160 + (-t176 * t28 + t179 * t29) * t161;
t1 = [Icges(4,1) * t178 ^ 2 + Icges(3,1) + Icges(2,3) + (Icges(5,1) * t161 - t235 - t236) * t161 + m(7) * (t44 ^ 2 + t45 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t86 ^ 2 + t87 ^ 2) + m(3) * (t122 ^ 2 + t123 ^ 2) + m(2) * (t144 ^ 2 + t146 ^ 2) + t238 + t239 + (-0.2e1 * Icges(4,4) * t178 + Icges(4,2) * t175) * t175 + (-0.2e1 * Icges(5,4) * t161 + Icges(5,2) * t160) * t160; m(7) * (t176 * t44 - t179 * t45) + m(6) * (t176 * t54 - t179 * t55) + m(5) * (t176 * t78 - t179 * t79) + m(4) * (t176 * t86 - t179 * t87) + m(3) * (t176 * t122 - t179 * t123); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t151 + t207; m(7) * (t50 * t44 + t51 * t45) + m(6) * (t65 * t54 + t66 * t55) + m(5) * (t97 * t78 + t98 * t79) + (t26 / 0.2e1 + t37 / 0.2e1 - t87 * t244 - t210 + t256 * t179) * t179 + (t27 / 0.2e1 + t38 / 0.2e1 + t86 * t244 + t209 + t256 * t176) * t176; m(5) * (t97 * t176 - t98 * t179) + m(6) * (t65 * t176 - t66 * t179) + m(7) * (t50 * t176 - t51 * t179) + t151 * t244; m(7) * (t16 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t32 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t49 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t145 ^ 2 * t151 + t77 ^ 2) + (t255 * t171 + t11 + t14) * t179 + (t12 + t15 + t253 * t170 + (t255 * t176 + t253 * t179) * t179) * t176; m(7) * (t176 * t45 + t179 * t44) + m(6) * (t176 * t55 + t179 * t54) + m(5) * (t176 * t79 + t179 * t78); 0; m(7) * (t176 * t51 + t179 * t50) + m(6) * (t176 * t66 + t179 * t65) + m(5) * (t176 * t98 + t179 * t97); t207; t46 + m(7) * (t24 * t45 + t25 * t44) + m(6) * (t52 * t55 + t53 * t54) + (t176 * t210 + t179 * t209) * t161 + t197; m(6) * (t53 * t176 - t52 * t179) + m(7) * (t25 * t176 - t24 * t179); (t34 * t176 + t33 * t179) * t250 + t7 * t245 + t8 * t247 + (t15 * t245 - t176 * t14 / 0.2e1) * t161 + m(7) * (t17 * t16 + t24 * t51 + t25 * t50) + m(6) * (t43 * t32 + t52 * t66 + t53 * t65) + t198; m(6) * (t52 * t176 + t53 * t179) + m(7) * (t24 * t176 + t25 * t179); t160 * t46 + m(7) * (t17 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t43 ^ 2 + t52 ^ 2 + t53 ^ 2) + ((t160 * t34 + t8) * t179 + (-t160 * t33 - t4 - t7) * t176) * t161 + t252; m(7) * (t48 * t44 + t47 * t45) + t197; m(7) * (t48 * t176 - t47 * t179); m(7) * (t39 * t16 + t47 * t51 + t48 * t50) + t198; m(7) * (t47 * t176 + t48 * t179); m(7) * (t39 * t17 + t47 * t24 + t48 * t25) + t196; m(7) * (t39 ^ 2 + t47 ^ 2 + t48 ^ 2) + t196;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
