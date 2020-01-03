% Calculate joint inertia matrix for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:27
% EndTime: 2019-12-31 21:32:35
% DurationCPUTime: 2.76s
% Computational Cost: add. (4918->398), mult. (12271->574), div. (0->0), fcn. (14526->8), ass. (0->195)
t175 = sin(qJ(2));
t250 = Icges(3,5) * t175;
t249 = t250 / 0.2e1;
t174 = sin(qJ(3));
t178 = cos(qJ(3));
t179 = cos(qJ(2));
t125 = -Icges(4,3) * t179 + (Icges(4,5) * t178 - Icges(4,6) * t174) * t175;
t128 = -Icges(5,2) * t179 + (Icges(5,4) * t178 + Icges(5,6) * t174) * t175;
t248 = t125 + t128;
t176 = sin(qJ(1));
t180 = cos(qJ(1));
t215 = t176 * t179;
t149 = t174 * t215 + t178 * t180;
t150 = -t174 * t180 + t178 * t215;
t218 = t175 * t176;
t85 = Icges(5,5) * t150 + Icges(5,6) * t218 + Icges(5,3) * t149;
t89 = Icges(5,4) * t150 + Icges(5,2) * t218 + Icges(5,6) * t149;
t93 = Icges(5,1) * t150 + Icges(5,4) * t218 + Icges(5,5) * t149;
t35 = t149 * t85 + t150 * t93 + t89 * t218;
t214 = t179 * t180;
t151 = t174 * t214 - t176 * t178;
t152 = t176 * t174 + t178 * t214;
t216 = t175 * t180;
t86 = Icges(5,5) * t152 + Icges(5,6) * t216 + Icges(5,3) * t151;
t90 = Icges(5,4) * t152 + Icges(5,2) * t216 + Icges(5,6) * t151;
t94 = Icges(5,1) * t152 + Icges(5,4) * t216 + Icges(5,5) * t151;
t36 = t149 * t86 + t150 * t94 + t90 * t218;
t87 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t218;
t91 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t218;
t95 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t218;
t37 = -t149 * t91 + t150 * t95 + t87 * t218;
t88 = Icges(4,5) * t152 - Icges(4,6) * t151 + Icges(4,3) * t216;
t92 = Icges(4,4) * t152 - Icges(4,2) * t151 + Icges(4,6) * t216;
t96 = Icges(4,1) * t152 - Icges(4,4) * t151 + Icges(4,5) * t216;
t38 = -t149 * t92 + t150 * t96 + t88 * t218;
t124 = -Icges(5,6) * t179 + (Icges(5,5) * t178 + Icges(5,3) * t174) * t175;
t132 = -Icges(5,4) * t179 + (Icges(5,1) * t178 + Icges(5,5) * t174) * t175;
t54 = t124 * t149 + t128 * t218 + t132 * t150;
t129 = -Icges(4,6) * t179 + (Icges(4,4) * t178 - Icges(4,2) * t174) * t175;
t133 = -Icges(4,5) * t179 + (Icges(4,1) * t178 - Icges(4,4) * t174) * t175;
t55 = t125 * t218 - t129 * t149 + t133 * t150;
t247 = (-t54 - t55) * t179 + ((t36 + t38) * t180 + (t35 + t37) * t176) * t175;
t39 = t151 * t85 + t152 * t93 + t89 * t216;
t40 = t151 * t86 + t152 * t94 + t90 * t216;
t41 = -t151 * t91 + t152 * t95 + t87 * t216;
t42 = -t151 * t92 + t152 * t96 + t88 * t216;
t56 = t151 * t124 + t128 * t216 + t152 * t132;
t57 = t125 * t216 - t151 * t129 + t152 * t133;
t246 = (-t56 - t57) * t179 + ((t40 + t42) * t180 + (t39 + t41) * t176) * t175;
t43 = -t179 * t89 + (t174 * t85 + t178 * t93) * t175;
t45 = -t179 * t87 + (-t174 * t91 + t178 * t95) * t175;
t245 = -t43 - t45;
t44 = -t179 * t90 + (t174 * t86 + t178 * t94) * t175;
t46 = -t179 * t88 + (-t174 * t92 + t178 * t96) * t175;
t244 = t44 + t46;
t219 = t174 * t175;
t217 = t175 * t178;
t242 = t124 * t219 + (t132 + t133) * t217;
t243 = (t129 * t219 + t248 * t179 - t242) * t179;
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t108 = t149 * t177 - t150 * t173;
t109 = t149 * t173 + t150 * t177;
t61 = Icges(6,5) * t109 + Icges(6,6) * t108 - Icges(6,3) * t218;
t63 = Icges(6,4) * t109 + Icges(6,2) * t108 - Icges(6,6) * t218;
t65 = Icges(6,1) * t109 + Icges(6,4) * t108 - Icges(6,5) * t218;
t12 = t108 * t63 + t109 * t65 - t61 * t218;
t110 = t151 * t177 - t152 * t173;
t111 = t151 * t173 + t152 * t177;
t62 = Icges(6,5) * t111 + Icges(6,6) * t110 - Icges(6,3) * t216;
t64 = Icges(6,4) * t111 + Icges(6,2) * t110 - Icges(6,6) * t216;
t66 = Icges(6,1) * t111 + Icges(6,4) * t110 - Icges(6,5) * t216;
t13 = t108 * t64 + t109 * t66 - t62 * t218;
t138 = (-t173 * t178 + t174 * t177) * t175;
t139 = (t173 * t174 + t177 * t178) * t175;
t81 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t179;
t82 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t179;
t83 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t179;
t27 = t108 * t82 + t109 * t83 - t81 * t218;
t1 = -t27 * t179 + (t12 * t176 + t13 * t180) * t175;
t14 = t110 * t63 + t111 * t65 - t61 * t216;
t15 = t110 * t64 + t111 * t66 - t62 * t216;
t28 = t110 * t82 + t111 * t83 - t81 * t216;
t2 = -t28 * t179 + (t14 * t176 + t15 * t180) * t175;
t21 = t138 * t63 + t139 * t65 + t179 * t61;
t22 = t138 * t64 + t139 * t66 + t179 * t62;
t204 = t138 * t82 + t139 * t83 + t179 * t81;
t31 = t204 * t179;
t3 = -t31 + (t21 * t176 + t22 * t180) * t175;
t241 = -(t176 * t1 + t180 * t2) * t175 + t179 * t3;
t240 = t176 ^ 2;
t239 = t180 ^ 2;
t238 = -t176 / 0.2e1;
t237 = t176 / 0.2e1;
t236 = -t179 / 0.2e1;
t235 = t179 / 0.2e1;
t234 = -t180 / 0.2e1;
t233 = t180 / 0.2e1;
t232 = -rSges(6,3) - pkin(8);
t159 = rSges(3,1) * t175 + rSges(3,2) * t179;
t231 = m(3) * t159;
t230 = pkin(2) * t179;
t227 = rSges(5,3) * t149;
t226 = t180 * rSges(3,3);
t115 = t152 * pkin(3) + t151 * qJ(4);
t99 = t152 * rSges(5,1) + rSges(5,2) * t216 + t151 * rSges(5,3);
t225 = -t115 - t99;
t191 = -rSges(6,1) * t109 - rSges(6,2) * t108;
t67 = -rSges(6,3) * t218 - t191;
t224 = pkin(4) * t150 - pkin(8) * t218 + t67;
t147 = t152 * pkin(4);
t213 = t111 * rSges(6,1) + t110 * rSges(6,2);
t68 = -rSges(6,3) * t216 + t213;
t223 = -pkin(8) * t216 + t147 + t68;
t84 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t179;
t222 = pkin(4) * t217 + pkin(8) * t179 + t84;
t220 = Icges(3,4) * t179;
t142 = t149 * qJ(4);
t114 = pkin(3) * t150 + t142;
t153 = (pkin(3) * t178 + qJ(4) * t174) * t175;
t212 = t179 * t114 + t153 * t218;
t136 = -t179 * rSges(5,2) + (rSges(5,1) * t178 + rSges(5,3) * t174) * t175;
t210 = -t136 - t153;
t137 = -t179 * rSges(4,3) + (rSges(4,1) * t178 - rSges(4,2) * t174) * t175;
t162 = pkin(2) * t175 - pkin(7) * t179;
t209 = -t137 - t162;
t207 = pkin(2) * t214 + pkin(7) * t216;
t208 = t240 * (pkin(7) * t175 + t230) + t180 * t207;
t206 = t180 * pkin(1) + t176 * pkin(6);
t171 = t180 * pkin(6);
t205 = t171 - t142;
t203 = -t21 / 0.2e1 - t27 / 0.2e1;
t202 = -t22 / 0.2e1 - t28 / 0.2e1;
t201 = -t115 - t223;
t200 = -t153 - t222;
t199 = -t162 + t210;
t100 = t152 * rSges(4,1) - t151 * rSges(4,2) + rSges(4,3) * t216;
t198 = -pkin(1) - t230;
t197 = -t162 + t200;
t196 = t176 * t114 + t180 * t115 + t208;
t195 = t206 + t207;
t193 = rSges(3,1) * t179 - rSges(3,2) * t175;
t192 = -rSges(4,1) * t150 + rSges(4,2) * t149;
t189 = -Icges(3,2) * t175 + t220;
t188 = Icges(3,5) * t179 - Icges(3,6) * t175;
t185 = rSges(3,1) * t214 - rSges(3,2) * t216 + t176 * rSges(3,3);
t184 = t115 + t195;
t183 = t1 * t233 + t2 * t238 + (t22 * t176 - t21 * t180) * t235;
t182 = t45 / 0.2e1 + t43 / 0.2e1 + t54 / 0.2e1 + t55 / 0.2e1 - t203;
t181 = t57 / 0.2e1 + t56 / 0.2e1 + t46 / 0.2e1 + t44 / 0.2e1 - t202;
t161 = rSges(2,1) * t180 - t176 * rSges(2,2);
t160 = -t176 * rSges(2,1) - rSges(2,2) * t180;
t156 = Icges(3,6) * t179 + t250;
t127 = Icges(3,3) * t176 + t188 * t180;
t126 = -Icges(3,3) * t180 + t188 * t176;
t117 = t185 + t206;
t116 = t226 + t171 + (-pkin(1) - t193) * t176;
t113 = t209 * t180;
t112 = t209 * t176;
t101 = t114 * t216;
t98 = rSges(4,3) * t218 - t192;
t97 = rSges(5,1) * t150 + rSges(5,2) * t218 + t227;
t80 = t180 * t185 + (t193 * t176 - t226) * t176;
t78 = t199 * t180;
t77 = t199 * t176;
t74 = t195 + t100;
t73 = t171 + ((-rSges(4,3) - pkin(7)) * t175 + t198) * t176 + t192;
t72 = -t179 * t100 - t137 * t216;
t71 = t137 * t218 + t179 * t98;
t60 = t197 * t180;
t59 = t197 * t176;
t58 = (-t100 * t176 + t180 * t98) * t175;
t53 = t184 + t99;
t52 = -t227 + (-rSges(5,1) - pkin(3)) * t150 + ((-rSges(5,2) - pkin(7)) * t175 + t198) * t176 + t205;
t51 = t100 * t180 + t176 * t98 + t208;
t50 = t225 * t179 + t210 * t216;
t49 = t136 * t218 + t179 * t97 + t212;
t48 = t179 * t68 + t84 * t216;
t47 = -t179 * t67 - t84 * t218;
t34 = t232 * t216 + t147 + t184 + t213;
t33 = (-pkin(3) - pkin(4)) * t150 + ((-pkin(7) - t232) * t175 + t198) * t176 + t191 + t205;
t32 = t101 + (t225 * t176 + t180 * t97) * t175;
t30 = (t176 * t68 - t180 * t67) * t175;
t29 = t176 * t97 + t180 * t99 + t196;
t24 = t201 * t179 + t200 * t216;
t23 = t224 * t179 + t222 * t218 + t212;
t20 = t101 + (t201 * t176 + t224 * t180) * t175;
t19 = t42 * t176 - t180 * t41;
t18 = t40 * t176 - t180 * t39;
t17 = t38 * t176 - t180 * t37;
t16 = t36 * t176 - t180 * t35;
t11 = t224 * t176 + t223 * t180 + t196;
t5 = -t14 * t180 + t15 * t176;
t4 = -t12 * t180 + t13 * t176;
t6 = [Icges(2,3) + (Icges(3,1) * t175 - t129 * t174 + t220) * t175 + (Icges(3,4) * t175 + Icges(3,2) * t179 - t248) * t179 + m(6) * (t33 ^ 2 + t34 ^ 2) + m(5) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2) + m(2) * (t160 ^ 2 + t161 ^ 2) + t204 + t242; m(6) * (t33 * t60 + t34 * t59) + m(5) * (t52 * t78 + t77 * t53) + m(4) * (t112 * t74 + t113 * t73) + ((-Icges(3,6) * t180 + t189 * t176) * t236 + t180 * t249 - t116 * t231 + t156 * t233 - t182) * t180 + ((Icges(3,6) * t176 + t189 * t180) * t235 + t176 * t249 - t117 * t231 + t156 * t237 + t181) * t176; m(6) * (t11 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t29 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t112 ^ 2 + t113 ^ 2 + t51 ^ 2) + m(3) * (t80 ^ 2 + (t239 + t240) * t159 ^ 2) + (-t239 * t126 - t16 - t17 - t4) * t180 + (t240 * t127 + t18 + t19 + t5 + (-t176 * t126 + t180 * t127) * t180) * t176; -t31 + t243 + m(6) * (t23 * t33 + t24 * t34) + m(5) * (t49 * t52 + t50 * t53) + m(4) * (t71 * t73 + t72 * t74) + (t182 * t176 + t181 * t180) * t175; m(6) * (t11 * t20 + t23 * t60 + t24 * t59) + m(5) * (t32 * t29 + t49 * t78 + t50 * t77) + m(4) * (t112 * t72 + t113 * t71 + t58 * t51) + ((t5 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1) * t180 + (t4 / 0.2e1 + t16 / 0.2e1 + t17 / 0.2e1) * t176) * t175 - t183 + t246 * t237 + (t244 * t176 + t245 * t180) * t236 + t247 * t234; m(6) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(4) * (t58 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t32 ^ 2 + t49 ^ 2 + t50 ^ 2) + (-t3 - t243) * t179 + ((-t244 * t179 + t2 + t246) * t180 + (t245 * t179 + t1 + t247) * t176) * t175; m(6) * (t149 * t34 + t151 * t33) + m(5) * (t149 * t53 + t151 * t52); m(6) * (t11 * t219 + t149 * t59 + t151 * t60) + m(5) * (t149 * t77 + t151 * t78 + t29 * t219); m(6) * (t149 * t24 + t151 * t23 + t20 * t219) + m(5) * (t149 * t50 + t151 * t49 + t32 * t219); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t174 ^ 2 * t175 ^ 2 + t149 ^ 2 + t151 ^ 2); m(6) * (t33 * t47 + t34 * t48) + t31 + (t203 * t176 + t202 * t180) * t175; m(6) * (t11 * t30 + t47 * t60 + t48 * t59) + (t5 * t234 + t4 * t238) * t175 + t183; m(6) * (t20 * t30 + t23 * t47 + t24 * t48) + t241; m(6) * (t149 * t48 + t151 * t47 + t30 * t219); m(6) * (t30 ^ 2 + t47 ^ 2 + t48 ^ 2) - t241;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
