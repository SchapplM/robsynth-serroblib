% Calculate joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP10_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:26
% EndTime: 2019-12-31 20:09:32
% DurationCPUTime: 2.12s
% Computational Cost: add. (2187->322), mult. (5054->463), div. (0->0), fcn. (5317->6), ass. (0->154)
t138 = -qJ(5) - pkin(7);
t222 = rSges(6,3) - t138;
t221 = Icges(4,1) + Icges(3,3);
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t220 = (-Icges(4,4) + Icges(3,5)) * t143 + (Icges(4,5) - Icges(3,6)) * t140;
t141 = sin(qJ(1));
t219 = -t141 / 0.2e1;
t218 = t141 / 0.2e1;
t144 = cos(qJ(1));
t217 = -t144 / 0.2e1;
t216 = t144 / 0.2e1;
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t71 = Icges(6,3) * t140 + (-Icges(6,5) * t139 - Icges(6,6) * t142) * t143;
t72 = Icges(5,3) * t140 + (-Icges(5,5) * t139 - Icges(5,6) * t142) * t143;
t215 = (t71 + t72) * t140;
t75 = Icges(6,6) * t140 + (-Icges(6,4) * t139 - Icges(6,2) * t142) * t143;
t76 = Icges(5,6) * t140 + (-Icges(5,4) * t139 - Icges(5,2) * t142) * t143;
t79 = Icges(6,5) * t140 + (-Icges(6,1) * t139 - Icges(6,4) * t142) * t143;
t80 = Icges(5,5) * t140 + (-Icges(5,1) * t139 - Icges(5,4) * t142) * t143;
t214 = (-t75 - t76) * t142 + (-t79 - t80) * t139;
t183 = t140 * t144;
t172 = t139 * t183;
t181 = t141 * t142;
t100 = t172 + t181;
t178 = t143 * t144;
t179 = t142 * t144;
t182 = t141 * t139;
t99 = t140 * t179 - t182;
t45 = Icges(6,5) * t100 + Icges(6,6) * t99 + Icges(6,3) * t178;
t49 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t178;
t53 = Icges(6,1) * t100 + Icges(6,4) * t99 + Icges(6,5) * t178;
t11 = t100 * t53 + t45 * t178 + t99 * t49;
t101 = t139 * t144 + t140 * t181;
t102 = t140 * t182 - t179;
t180 = t141 * t143;
t46 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t180;
t50 = Icges(6,4) * t102 + Icges(6,2) * t101 + Icges(6,6) * t180;
t54 = Icges(6,1) * t102 + Icges(6,4) * t101 + Icges(6,5) * t180;
t12 = t100 * t54 + t46 * t178 + t99 * t50;
t47 = Icges(5,5) * t100 + Icges(5,6) * t99 + Icges(5,3) * t178;
t51 = Icges(5,4) * t100 + Icges(5,2) * t99 + Icges(5,6) * t178;
t55 = Icges(5,1) * t100 + Icges(5,4) * t99 + Icges(5,5) * t178;
t13 = t100 * t55 + t47 * t178 + t99 * t51;
t48 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t180;
t52 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t180;
t56 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t180;
t14 = t100 * t56 + t48 * t178 + t99 * t52;
t26 = t100 * t79 + t71 * t178 + t99 * t75;
t27 = t100 * t80 + t72 * t178 + t99 * t76;
t213 = ((t11 + t13) * t144 + (t12 + t14) * t141) * t143 + (t26 + t27) * t140;
t15 = t101 * t49 + t102 * t53 + t45 * t180;
t16 = t101 * t50 + t102 * t54 + t46 * t180;
t17 = t101 * t51 + t102 * t55 + t47 * t180;
t18 = t101 * t52 + t102 * t56 + t48 * t180;
t28 = t101 * t75 + t102 * t79 + t71 * t180;
t29 = t101 * t76 + t102 * t80 + t72 * t180;
t212 = ((t15 + t17) * t144 + (t16 + t18) * t141) * t143 + (t28 + t29) * t140;
t211 = t140 / 0.2e1;
t20 = t140 * t45 + (-t139 * t53 - t142 * t49) * t143;
t22 = t140 * t47 + (-t139 * t55 - t142 * t51) * t143;
t210 = t20 + t22;
t21 = t140 * t46 + (-t139 * t54 - t142 * t50) * t143;
t23 = t140 * t48 + (-t139 * t56 - t142 * t52) * t143;
t209 = t21 + t23;
t208 = -t220 * t141 + t221 * t144;
t207 = t221 * t141 + t220 * t144;
t126 = pkin(4) * t142 + pkin(3);
t206 = t100 * rSges(6,1) + t99 * rSges(6,2) + pkin(4) * t172 + t141 * t126 + t222 * t178;
t205 = -t102 * rSges(6,1) - t101 * rSges(6,2) + t126 * t144;
t135 = t141 ^ 2;
t137 = t144 ^ 2;
t204 = m(4) / 0.2e1;
t203 = m(5) / 0.2e1;
t114 = rSges(3,1) * t140 + rSges(3,2) * t143;
t199 = m(3) * t114;
t198 = pkin(4) * t139;
t197 = -pkin(7) - t138;
t196 = (t214 * t143 + t215) * t140;
t174 = t141 * pkin(3) + pkin(7) * t178;
t195 = -t174 + t206;
t132 = t144 * pkin(3);
t194 = rSges(6,3) * t180 + t132 + (t140 * t198 + t197 * t143) * t141 - t205;
t193 = (-rSges(6,1) * t139 - rSges(6,2) * t142 - t198) * t143 + (rSges(6,3) + t197) * t140;
t176 = pkin(2) * t178 + qJ(3) * t183;
t185 = qJ(3) * t140;
t192 = t135 * (pkin(2) * t143 + t185) + t144 * t176;
t191 = t144 * rSges(4,1);
t190 = t144 * rSges(3,3);
t189 = Icges(3,4) * t140;
t188 = Icges(3,4) * t143;
t187 = Icges(4,6) * t140;
t186 = Icges(4,6) * t143;
t112 = pkin(2) * t140 - qJ(3) * t143;
t177 = rSges(4,2) * t140 + rSges(4,3) * t143 - t112;
t175 = t144 * pkin(1) + t141 * pkin(6);
t173 = t135 + t137;
t58 = t100 * rSges(5,1) + t99 * rSges(5,2) + rSges(5,3) * t178;
t171 = Icges(3,5) * t211 - Icges(4,4) * t140 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,5) / 0.2e1) * t143;
t170 = -pkin(7) * t140 - t112;
t169 = t141 * (pkin(7) * t180 - t132) + t144 * t174 + t192;
t168 = t175 + t176;
t90 = rSges(5,3) * t140 + (-rSges(5,1) * t139 - rSges(5,2) * t142) * t143;
t167 = t170 - t90;
t166 = rSges(3,1) * t143 - rSges(3,2) * t140;
t165 = -rSges(5,1) * t102 - rSges(5,2) * t101;
t24 = -t194 * t140 + t193 * t180;
t25 = t195 * t140 - t193 * t178;
t159 = t141 * t25 + t144 * t24;
t157 = t170 - t193;
t40 = t157 * t141;
t41 = t157 * t144;
t158 = t141 * t40 + t144 * t41;
t156 = Icges(3,1) * t143 - t189;
t155 = -Icges(3,2) * t140 + t188;
t152 = -Icges(4,2) * t143 + t187;
t151 = Icges(4,3) * t140 - t186;
t150 = rSges(3,1) * t178 - rSges(3,2) * t183 + t141 * rSges(3,3);
t149 = t141 * rSges(4,1) - rSges(4,2) * t178 + rSges(4,3) * t183;
t147 = t27 / 0.2e1 + t26 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t146 = t29 / 0.2e1 + t28 / 0.2e1 + t23 / 0.2e1 + t21 / 0.2e1;
t131 = t144 * pkin(6);
t31 = t131 + (-pkin(1) + (-qJ(3) - t198) * t140 + (-pkin(2) - t222) * t143) * t141 + t205;
t32 = t168 + t206;
t145 = m(6) * (t141 * t32 + t144 * t31);
t136 = t143 ^ 2;
t134 = t140 ^ 2;
t116 = rSges(2,1) * t144 - t141 * rSges(2,2);
t115 = -t141 * rSges(2,1) - rSges(2,2) * t144;
t68 = t177 * t144;
t67 = t177 * t141;
t66 = t150 + t175;
t65 = t190 + t131 + (-pkin(1) - t166) * t141;
t62 = t149 + t168;
t61 = t191 + t131 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t143 + (-rSges(4,3) - qJ(3)) * t140) * t141;
t60 = rSges(5,3) * t180 - t165;
t44 = t167 * t144;
t43 = t167 * t141;
t42 = t144 * t150 + (t141 * t166 - t190) * t141;
t39 = t140 * t58 - t90 * t178;
t38 = -t140 * t60 + t90 * t180;
t37 = t168 + t58 + t174;
t36 = t131 + t132 + (-t185 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(7)) * t143) * t141 + t165;
t33 = t144 * t149 + (-t191 + (-rSges(4,2) * t143 + rSges(4,3) * t140) * t141) * t141 + t192;
t30 = (-t141 * t58 + t144 * t60) * t143;
t19 = t141 * t60 + t144 * t58 + t169;
t10 = (-t195 * t141 + t194 * t144) * t143;
t9 = t194 * t141 + t195 * t144 + t169;
t8 = t17 * t141 - t144 * t18;
t7 = t15 * t141 - t144 * t16;
t6 = t13 * t141 - t14 * t144;
t5 = t11 * t141 - t12 * t144;
t1 = [Icges(2,3) + m(6) * (t31 ^ 2 + t32 ^ 2) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t65 ^ 2 + t66 ^ 2) + m(4) * (t61 ^ 2 + t62 ^ 2) + m(2) * (t115 ^ 2 + t116 ^ 2) + (t187 + t189 + (Icges(4,3) + Icges(3,2)) * t143 + t214) * t143 + (t186 + t188 + (Icges(3,1) + Icges(4,2)) * t140) * t140 + t215; m(6) * (t31 * t41 + t32 * t40) + m(5) * (t36 * t44 + t37 * t43) + m(4) * (t61 * t68 + t62 * t67) + (-t65 * t199 + t171 * t144 + (Icges(4,5) * t217 + Icges(3,6) * t216 + t151 * t218 + t155 * t219) * t143 + (Icges(4,4) * t217 + Icges(3,5) * t216 + t152 * t218 + t156 * t219) * t140 - t146) * t144 + (-t66 * t199 + t171 * t141 + (Icges(4,5) * t219 + Icges(3,6) * t218 + t151 * t217 + t155 * t216) * t143 + (Icges(4,4) * t219 + Icges(3,5) * t218 + t152 * t217 + t156 * t216) * t140 + t147) * t141; m(6) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(5) * (t19 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(4) * (t33 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t114 ^ 2 * t173 + t42 ^ 2) + (t137 * t208 - t7 - t8) * t144 + (t5 + t6 + t207 * t135 + (t141 * t208 + t144 * t207) * t144) * t141; 0.2e1 * (t145 / 0.2e1 + (t141 * t37 + t144 * t36) * t203 + (t141 * t62 + t144 * t61) * t204) * t140; m(6) * (t140 * t158 - t143 * t9) + m(5) * (-t143 * t19 + (t141 * t43 + t144 * t44) * t140) + m(4) * (-t143 * t33 + (t141 * t67 + t144 * t68) * t140); 0.2e1 * (t204 + t203 + m(6) / 0.2e1) * (t134 * t173 + t136); m(6) * (t24 * t31 + t25 * t32) + m(5) * (t36 * t38 + t37 * t39) + (t141 * t146 + t144 * t147) * t143 + t196; m(6) * (t10 * t9 + t24 * t41 + t25 * t40) + m(5) * (t19 * t30 + t38 * t44 + t39 * t43) + ((t6 / 0.2e1 + t5 / 0.2e1) * t144 + (t7 / 0.2e1 + t8 / 0.2e1) * t141) * t143 + (t141 * t210 - t144 * t209) * t211 + t213 * t218 + t212 * t217; m(5) * (-t30 * t143 + (t141 * t39 + t144 * t38) * t140) + m(6) * (-t10 * t143 + t140 * t159); t196 * t140 + m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t30 ^ 2 + t38 ^ 2 + t39 ^ 2) + (t213 * t144 + t212 * t141 + (t141 * t209 + t144 * t210) * t140) * t143; t143 * t145; m(6) * (t140 * t9 + t143 * t158); m(6) * (-0.1e1 + t173) * t143 * t140; m(6) * (t140 * t10 + t143 * t159); m(6) * (t136 * t173 + t134);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
