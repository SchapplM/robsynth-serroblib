% Calculate joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:44:54
% DurationCPUTime: 2.79s
% Computational Cost: add. (3961->249), mult. (3520->357), div. (0->0), fcn. (3246->8), ass. (0->138)
t231 = Icges(5,4) + Icges(6,4);
t230 = Icges(5,1) + Icges(6,1);
t229 = Icges(5,2) + Icges(6,2);
t131 = qJ(2) + qJ(3);
t121 = qJ(4) + t131;
t116 = cos(t121);
t228 = t231 * t116;
t115 = sin(t121);
t227 = t231 * t115;
t226 = Icges(5,5) + Icges(6,5);
t225 = Icges(5,6) + Icges(6,6);
t224 = -t229 * t115 + t228;
t223 = t230 * t116 - t227;
t133 = sin(qJ(1));
t135 = cos(qJ(1));
t222 = t224 * t133 - t225 * t135;
t221 = t225 * t133 + t224 * t135;
t220 = t223 * t133 - t226 * t135;
t219 = t226 * t133 + t223 * t135;
t218 = Icges(5,3) + Icges(6,3);
t217 = -t225 * t115 + t226 * t116;
t216 = t229 * t116 + t227;
t215 = t230 * t115 + t228;
t214 = -t217 * t133 + t218 * t135;
t213 = t218 * t133 + t217 * t135;
t212 = t222 * t115 - t220 * t116;
t211 = t221 * t115 - t219 * t116;
t129 = t133 ^ 2;
t210 = t133 * pkin(6);
t209 = t226 * t115 + t225 * t116;
t130 = t135 ^ 2;
t118 = sin(t131);
t119 = cos(t131);
t188 = Icges(4,4) * t119;
t152 = -Icges(4,2) * t118 + t188;
t75 = Icges(4,6) * t133 + t152 * t135;
t189 = Icges(4,4) * t118;
t156 = Icges(4,1) * t119 - t189;
t77 = Icges(4,5) * t133 + t156 * t135;
t162 = -t118 * t75 + t119 * t77;
t74 = -Icges(4,6) * t135 + t152 * t133;
t76 = -Icges(4,5) * t135 + t156 * t133;
t163 = t118 * t74 - t119 * t76;
t203 = t214 * t130 + (t211 * t133 + (-t212 + t213) * t135) * t133;
t148 = Icges(4,5) * t119 - Icges(4,6) * t118;
t72 = -Icges(4,3) * t135 + t148 * t133;
t73 = Icges(4,3) * t133 + t148 * t135;
t208 = -t130 * t72 - (t162 * t133 + (t163 - t73) * t135) * t133 + t203;
t182 = t116 * t135;
t183 = t115 * t135;
t134 = cos(qJ(2));
t117 = t134 * pkin(2) + pkin(1);
t101 = pkin(3) * t119 + t117;
t85 = pkin(4) * t116 + t101;
t207 = rSges(6,1) * t182 - rSges(6,2) * t183 + t133 * rSges(6,3) + t135 * t85;
t206 = -rSges(6,1) * t116 + rSges(6,2) * t115 - t85;
t205 = -t216 * t115 + t215 * t116;
t172 = rSges(4,1) * t119 - rSges(4,2) * t118;
t136 = -pkin(7) - pkin(6);
t204 = (t213 * t129 + ((-t211 + t214) * t133 + t212 * t135) * t135) * t133;
t202 = t133 / 0.2e1;
t201 = -t135 / 0.2e1;
t132 = sin(qJ(2));
t200 = pkin(2) * t132;
t199 = pkin(3) * t118;
t111 = t135 * t117;
t95 = t135 * t101;
t198 = t135 * (-t111 + t95) + (t101 - t117) * t129;
t144 = rSges(5,1) * t182 - rSges(5,2) * t183 + t133 * rSges(5,3);
t171 = rSges(5,1) * t116 - rSges(5,2) * t115;
t25 = t133 * (-t135 * rSges(5,3) + t171 * t133) + t135 * t144;
t127 = t135 * pkin(6);
t197 = t133 * (t127 + (-pkin(1) + t117) * t133) + t135 * (-t135 * pkin(1) + t111 - t210);
t145 = t133 * rSges(4,3) + t172 * t135;
t30 = t133 * (-t135 * rSges(4,3) + t172 * t133) + t135 * t145;
t196 = rSges(3,1) * t134;
t194 = rSges(3,2) * t132;
t192 = t135 * rSges(3,3);
t191 = Icges(3,4) * t132;
t190 = Icges(3,4) * t134;
t181 = t133 * rSges(3,3) + t135 * t196;
t178 = t129 + t130;
t128 = -pkin(8) + t136;
t177 = t133 * (t129 * t73 + (t163 * t135 + (t162 - t72) * t133) * t135) + t204;
t94 = t115 * rSges(5,1) + t116 * rSges(5,2);
t176 = -t94 - t199;
t175 = -t116 * rSges(6,2) + (-rSges(6,1) - pkin(4)) * t115;
t100 = t118 * rSges(4,1) + t119 * rSges(4,2);
t174 = -t100 - t200;
t10 = (-t95 + t207) * t135 + (-t135 * rSges(6,3) + (-t101 - t206) * t133) * t133;
t11 = t25 + t198;
t173 = -t194 + t196;
t98 = Icges(4,2) * t119 + t189;
t99 = Icges(4,1) * t118 + t188;
t161 = -t118 * t98 + t119 * t99;
t158 = t203 * t135 + t204;
t2 = t10 + t198;
t157 = Icges(3,1) * t134 - t191;
t153 = -Icges(3,2) * t132 + t190;
t149 = Icges(3,5) * t134 - Icges(3,6) * t132;
t142 = t176 - t200;
t141 = t175 - t199;
t140 = t208 * t135 + t177;
t139 = (t219 * t115 + t221 * t116 + t209 * t133 + t205 * t135) * t202 + (t220 * t115 + t222 * t116 + t205 * t133 - t209 * t135) * t201;
t138 = t141 - t200;
t97 = Icges(4,5) * t118 + Icges(4,6) * t119;
t137 = t139 + (t118 * t77 + t119 * t75 + t133 * t97 + t161 * t135) * t202 + (t118 * t76 + t119 * t74 + t161 * t133 - t135 * t97) * t201;
t120 = -qJ(5) + t128;
t110 = t135 * rSges(2,1) - t133 * rSges(2,2);
t109 = -t133 * rSges(2,1) - t135 * rSges(2,2);
t108 = t132 * rSges(3,1) + t134 * rSges(3,2);
t80 = Icges(3,3) * t133 + t149 * t135;
t79 = -Icges(3,3) * t135 + t149 * t133;
t71 = t174 * t135;
t70 = t174 * t133;
t53 = t210 + (pkin(1) - t194) * t135 + t181;
t52 = t192 + t127 + (-pkin(1) - t173) * t133;
t51 = t176 * t135;
t50 = t176 * t133;
t45 = t175 * t135;
t44 = t175 * t133;
t43 = t142 * t135;
t42 = t142 * t133;
t41 = -t133 * t136 + t111 + t145;
t40 = (rSges(4,3) - t136) * t135 + (-t117 - t172) * t133;
t39 = t141 * t135;
t38 = t141 * t133;
t35 = t138 * t135;
t34 = t138 * t133;
t33 = t135 * (-t135 * t194 + t181) + (t173 * t133 - t192) * t133;
t32 = -t133 * t128 + t144 + t95;
t31 = (rSges(5,3) - t128) * t135 + (-t101 - t171) * t133;
t29 = -t133 * t120 + t207;
t28 = (rSges(6,3) - t120) * t135 + t206 * t133;
t12 = t30 + t197;
t7 = t11 + t197;
t1 = t2 + t197;
t3 = [t134 * (Icges(3,2) * t134 + t191) + t132 * (Icges(3,1) * t132 + t190) + t118 * t99 + t119 * t98 + Icges(2,3) + t216 * t116 + t215 * t115 + m(6) * (t28 ^ 2 + t29 ^ 2) + m(5) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t52 ^ 2 + t53 ^ 2) + m(2) * (t109 ^ 2 + t110 ^ 2); (t130 / 0.2e1 + t129 / 0.2e1) * (Icges(3,5) * t132 + Icges(3,6) * t134) + m(6) * (t35 * t28 + t34 * t29) + m(5) * (t43 * t31 + t42 * t32) + m(4) * (t71 * t40 + t70 * t41) + m(3) * (-t133 * t53 - t135 * t52) * t108 + (t132 * (-Icges(3,5) * t135 + t157 * t133) + t134 * (-Icges(3,6) * t135 + t153 * t133)) * t201 + (t132 * (Icges(3,5) * t133 + t157 * t135) + t134 * (Icges(3,6) * t133 + t153 * t135)) * t202 + t137; m(6) * (t1 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2 + t7 ^ 2) + m(4) * (t12 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(3) * (t178 * t108 ^ 2 + t33 ^ 2) + t133 * t129 * t80 + t177 + (-t130 * t79 + (-t133 * t79 + t135 * t80) * t133 + t208) * t135; m(4) * (-t133 * t41 - t135 * t40) * t100 + m(6) * (t39 * t28 + t38 * t29) + m(5) * (t51 * t31 + t50 * t32) + t137; m(6) * (t2 * t1 + t38 * t34 + t39 * t35) + m(5) * (t11 * t7 + t50 * t42 + t51 * t43) + m(4) * (t30 * t12 + (-t133 * t70 - t135 * t71) * t100) + t140; m(6) * (t2 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t11 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t178 * t100 ^ 2 + t30 ^ 2) + t140; m(6) * (t45 * t28 + t44 * t29) + m(5) * (-t133 * t32 - t135 * t31) * t94 + t139; m(6) * (t10 * t1 + t44 * t34 + t45 * t35) + m(5) * (t25 * t7 + (-t133 * t42 - t135 * t43) * t94) + t158; m(6) * (t10 * t2 + t44 * t38 + t45 * t39) + m(5) * (t25 * t11 + (-t133 * t50 - t135 * t51) * t94) + t158; m(5) * (t178 * t94 ^ 2 + t25 ^ 2) + m(6) * (t10 ^ 2 + t44 ^ 2 + t45 ^ 2) + t158; m(6) * (t133 * t28 - t135 * t29); m(6) * (t133 * t35 - t135 * t34); m(6) * (t133 * t39 - t135 * t38); m(6) * (t133 * t45 - t135 * t44); m(6) * t178;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
