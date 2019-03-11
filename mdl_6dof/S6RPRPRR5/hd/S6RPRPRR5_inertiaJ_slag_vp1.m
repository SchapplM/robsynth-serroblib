% Calculate joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:00
% EndTime: 2019-03-09 03:48:06
% DurationCPUTime: 2.90s
% Computational Cost: add. (7099->335), mult. (9283->511), div. (0->0), fcn. (10871->10), ass. (0->170)
t242 = Icges(4,1) + Icges(5,1);
t240 = Icges(4,5) + Icges(5,4);
t135 = pkin(10) + qJ(3);
t131 = sin(t135);
t241 = (Icges(4,4) - Icges(5,5)) * t131;
t239 = Icges(5,2) + Icges(4,3);
t132 = cos(t135);
t238 = t240 * t132 + (-Icges(4,6) + Icges(5,6)) * t131;
t237 = t242 * t132 - t241;
t143 = sin(qJ(1));
t236 = -t143 / 0.2e1;
t230 = t143 / 0.2e1;
t145 = cos(qJ(1));
t235 = -t145 / 0.2e1;
t234 = t145 / 0.2e1;
t142 = sin(qJ(5));
t221 = cos(qJ(5));
t105 = t131 * t142 + t132 * t221;
t100 = t105 * t143;
t194 = Icges(6,3) * t145;
t195 = Icges(6,5) * t145;
t200 = Icges(6,1) * t100;
t182 = t131 * t221;
t191 = t132 * t143;
t99 = t142 * t191 - t143 * t182;
t205 = Icges(6,6) * t99;
t210 = Icges(6,4) * t99;
t226 = t99 ^ 2;
t190 = t132 * t145;
t101 = t142 * t190 - t145 * t182;
t102 = t105 * t145;
t141 = sin(qJ(6));
t144 = cos(qJ(6));
t76 = -t102 * t141 - t143 * t144;
t77 = t102 * t144 - t141 * t143;
t39 = Icges(7,5) * t77 + Icges(7,6) * t76 + Icges(7,3) * t101;
t40 = Icges(7,4) * t77 + Icges(7,2) * t76 + Icges(7,6) * t101;
t41 = Icges(7,1) * t77 + Icges(7,4) * t76 + Icges(7,5) * t101;
t74 = -t100 * t141 + t144 * t145;
t75 = t100 * t144 + t141 * t145;
t11 = t39 * t99 + t40 * t74 + t41 * t75;
t204 = Icges(7,6) * t99;
t206 = Icges(7,2) * t74;
t209 = Icges(7,4) * t75;
t148 = (0.2e1 * t204 + t206 + 0.2e1 * t209) * t74;
t207 = Icges(7,5) * t99;
t211 = Icges(7,1) * t75;
t3 = t11 * t143 - (Icges(7,3) * t226 + (0.2e1 * t207 + t211) * t75 + t148) * t145;
t55 = Icges(6,5) * t102 - Icges(6,6) * t101 - Icges(6,3) * t143;
t56 = Icges(6,4) * t102 - Icges(6,2) * t101 - Icges(6,6) * t143;
t57 = Icges(6,1) * t102 - Icges(6,4) * t101 - Icges(6,5) * t143;
t223 = -(t100 * t57 + t145 * t55 - t99 * t56) * t143 + (t226 * Icges(6,2) + (t194 - 0.2e1 * t205) * t145 + (0.2e1 * t195 + t200 - 0.2e1 * t210) * t100) * t145 - t3;
t150 = Icges(6,4) * t100 - Icges(6,2) * t99 + Icges(6,6) * t145;
t151 = t195 + t200 - t210;
t12 = t101 * t39 + t40 * t76 + t41 * t77;
t203 = Icges(7,3) * t99;
t208 = Icges(7,5) * t75;
t152 = Icges(7,6) * t74 + t203 + t208;
t153 = t204 + t206 + t209;
t154 = Icges(7,4) * t74 + t207 + t211;
t146 = t101 * t152 + t76 * t153 + t77 * t154;
t4 = t12 * t143 - t146 * t145;
t222 = (-t101 * t56 + t102 * t57 - t143 * t55) * t143 - (t102 * t151 - t101 * t150 - t143 * (Icges(6,5) * t100 + t194 - t205)) * t145 + t4;
t233 = t101 / 0.2e1;
t232 = t105 / 0.2e1;
t229 = t239 * t143 + t238 * t145;
t228 = -t238 * t143 + t239 * t145;
t43 = t77 * rSges(7,1) + t76 * rSges(7,2) + t101 * rSges(7,3);
t216 = t102 * pkin(5) + pkin(9) * t101 + t43;
t219 = t100 * pkin(5);
t174 = -t75 * rSges(7,1) - t74 * rSges(7,2);
t42 = t99 * rSges(7,3) - t174;
t15 = -(t99 * pkin(9) + t219 + t42) * t143 - t216 * t145;
t227 = -t222 * t143 - t223 * t145;
t136 = t143 ^ 2;
t137 = t145 ^ 2;
t225 = m(5) / 0.2e1;
t224 = m(7) / 0.2e1;
t116 = rSges(4,1) * t131 + rSges(4,2) * t132;
t220 = m(4) * t116;
t140 = -pkin(7) - qJ(2);
t218 = -pkin(8) - t140;
t106 = -t132 * t142 + t182;
t48 = Icges(7,3) * t105 + (Icges(7,5) * t144 - Icges(7,6) * t141) * t106;
t50 = Icges(7,5) * t105 + (Icges(7,1) * t144 - Icges(7,4) * t141) * t106;
t215 = t106 * t144 * t50 + t105 * t48;
t51 = rSges(7,3) * t105 + (rSges(7,1) * t144 - rSges(7,2) * t141) * t106;
t214 = pkin(5) * t106 + pkin(9) * t105 + t51;
t213 = t102 * rSges(6,1) - t101 * rSges(6,2);
t192 = t131 * t145;
t188 = pkin(3) * t190 + qJ(4) * t192;
t193 = qJ(4) * t131;
t212 = t136 * (pkin(3) * t132 + t193) + t145 * t188;
t49 = Icges(7,6) * t105 + (Icges(7,4) * t144 - Icges(7,2) * t141) * t106;
t202 = t141 * t49;
t201 = rSges(3,3) + qJ(2);
t198 = Icges(4,4) * t132;
t196 = Icges(5,5) * t132;
t114 = pkin(3) * t131 - qJ(4) * t132;
t189 = -rSges(5,1) * t131 + rSges(5,3) * t132 - t114;
t187 = t136 + t137;
t186 = -rSges(6,3) + t218;
t13 = t105 * t152 + (-t141 * t153 + t144 * t154) * t106;
t16 = t48 * t99 + t49 * t74 + t50 * t75;
t185 = t13 / 0.2e1 + t16 / 0.2e1;
t14 = t105 * t39 + (-t141 * t40 + t144 * t41) * t106;
t17 = t101 * t48 + t49 * t76 + t50 * t77;
t184 = t17 / 0.2e1 + t14 / 0.2e1;
t183 = rSges(5,1) * t190 + t143 * rSges(5,2) + rSges(5,3) * t192;
t181 = t240 * t131 / 0.2e1 + (Icges(4,6) / 0.2e1 - Icges(5,6) / 0.2e1) * t132;
t180 = -pkin(4) * t131 - t114;
t139 = cos(pkin(10));
t129 = pkin(2) * t139 + pkin(1);
t122 = t145 * t129;
t179 = -t143 * t140 + t122;
t126 = pkin(4) * t190;
t178 = t143 * pkin(4) * t191 + t145 * t126 + t212;
t177 = t225 + m(6) / 0.2e1 + t224;
t176 = t122 + t126 + t188;
t68 = rSges(6,1) * t106 - rSges(6,2) * t105;
t175 = t180 - t68;
t173 = -t100 * rSges(6,1) + t99 * rSges(6,2);
t172 = rSges(4,1) * t132 - rSges(4,2) * t131;
t46 = t175 * t143;
t47 = t175 * t145;
t167 = t143 * t46 + t145 * t47;
t33 = t143 * t173 - t145 * t213;
t166 = t180 - t214;
t163 = -Icges(4,2) * t131 + t198;
t160 = Icges(5,3) * t131 + t196;
t159 = rSges(4,1) * t190 - rSges(4,2) * t192 + t143 * rSges(4,3);
t138 = sin(pkin(10));
t158 = rSges(3,1) * t139 - rSges(3,2) * t138 + pkin(1);
t65 = Icges(6,5) * t106 - Icges(6,6) * t105;
t66 = Icges(6,4) * t106 - Icges(6,2) * t105;
t67 = Icges(6,1) * t106 - Icges(6,4) * t105;
t157 = -t105 * t150 / 0.2e1 + t106 * t151 / 0.2e1 + t100 * t67 / 0.2e1 + t65 * t234 - t99 * t66 / 0.2e1 + t185;
t156 = t66 * t233 - t102 * t67 / 0.2e1 + t65 * t230 + t56 * t232 - t106 * t57 / 0.2e1 - t184;
t149 = (-t193 - t129 + (-pkin(3) - pkin(4)) * t132) * t143;
t36 = t186 * t145 + t149 + t173;
t37 = t186 * t143 + t176 + t213;
t155 = m(6) * (t143 * t37 + t145 * t36);
t1 = t11 * t101 + t16 * t105 + (Icges(7,1) * t75 ^ 2 + (t203 + 0.2e1 * t208) * t99 + t148) * t99;
t2 = t12 * t101 + t17 * t105 + t146 * t99;
t147 = t4 * t233 + (-t13 * t145 + t14 * t143) * t232 + t2 * t230 + t1 * t235 + t99 * t3 / 0.2e1;
t120 = rSges(2,1) * t145 - t143 * rSges(2,2);
t119 = -t143 * rSges(2,1) - rSges(2,2) * t145;
t79 = t201 * t143 + t158 * t145;
t78 = -t158 * t143 + t201 * t145;
t71 = t189 * t145;
t70 = t189 * t143;
t63 = t159 + t179;
t62 = (rSges(4,3) - t140) * t145 + (-t129 - t172) * t143;
t54 = t145 * t159 + (-t145 * rSges(4,3) + t172 * t143) * t143;
t53 = t179 + t183 + t188;
t52 = (rSges(5,2) - t140) * t145 + (-t129 + (-rSges(5,1) - pkin(3)) * t132 + (-rSges(5,3) - qJ(4)) * t131) * t143;
t38 = t145 * t183 + (-t145 * rSges(5,2) + (rSges(5,1) * t132 + rSges(5,3) * t131) * t143) * t143 + t212;
t35 = t214 * t145;
t34 = t214 * t143;
t32 = t166 * t145;
t31 = t166 * t143;
t26 = t218 * t143 + t176 + t216;
t25 = -t219 + (-rSges(7,3) - pkin(9)) * t99 + t218 * t145 + t149 + t174;
t24 = -t33 + t178;
t23 = -t101 * t51 + t105 * t43;
t22 = -t105 * t42 + t51 * t99;
t19 = t101 * t42 - t43 * t99;
t18 = (-t106 * t202 + t215) * t105;
t8 = -t15 + t178;
t5 = [Icges(3,2) * t139 ^ 2 - t105 * t66 + Icges(2,3) + (Icges(3,1) * t138 + 0.2e1 * Icges(3,4) * t139) * t138 + (t67 - t202) * t106 + m(7) * (t25 ^ 2 + t26 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2) + m(5) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(2) * (t119 ^ 2 + t120 ^ 2) + t215 + ((Icges(4,2) + Icges(5,3)) * t132 + t241) * t132 + (t242 * t131 - t196 + t198) * t131; m(7) * (t143 * t25 - t145 * t26) + m(6) * (t143 * t36 - t145 * t37) + m(5) * (t143 * t52 - t145 * t53) + m(4) * (t143 * t62 - t145 * t63) + m(3) * (t143 * t78 - t145 * t79); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t177) * t187; m(7) * (t25 * t32 + t26 * t31) + m(6) * (t36 * t47 + t37 * t46) + m(5) * (t52 * t71 + t53 * t70) + (-t62 * t220 + t181 * t145 + (Icges(4,6) * t234 + Icges(5,6) * t235 + t160 * t230 + t163 * t236) * t132 + (t240 * t234 + t237 * t236) * t131 - t157) * t145 + (-t63 * t220 + t181 * t143 + (Icges(4,6) * t230 + Icges(5,6) * t236 + t160 * t235 + t163 * t234) * t132 + (t240 * t230 + t237 * t234) * t131 - t156) * t143; m(5) * (t71 * t143 - t145 * t70) + m(6) * (t47 * t143 - t145 * t46) + m(7) * (t32 * t143 - t145 * t31); m(7) * (t31 ^ 2 + t32 ^ 2 + t8 ^ 2) + m(6) * (t24 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t38 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t187 * t116 ^ 2 + t54 ^ 2) + (t137 * t228 + t223) * t145 + (t229 * t136 + (t143 * t228 + t145 * t229) * t145 + t222) * t143; 0.2e1 * ((t143 * t26 + t145 * t25) * t224 + t155 / 0.2e1 + (t143 * t53 + t145 * t52) * t225) * t131; 0; m(7) * (-t132 * t8 + (t143 * t31 + t145 * t32) * t131) + m(6) * (t167 * t131 - t132 * t24) + m(5) * (-t132 * t38 + (t143 * t70 + t145 * t71) * t131); 0.2e1 * t177 * (t187 * t131 ^ 2 + t132 ^ 2); t157 * t145 + t156 * t143 + m(7) * (t25 * t35 + t26 * t34) + t68 * t155; m(7) * (t35 * t143 - t145 * t34); m(7) * (t15 * t8 + t31 * t34 + t32 * t35) + m(6) * (t167 * t68 + t33 * t24) + t227; m(6) * (t187 * t68 * t131 - t132 * t33) + m(7) * (-t15 * t132 + (t143 * t34 + t145 * t35) * t131); m(6) * (t187 * t68 ^ 2 + t33 ^ 2) + m(7) * (t15 ^ 2 + t34 ^ 2 + t35 ^ 2) - t227; m(7) * (t22 * t25 + t23 * t26) + t18 + t185 * t99 + t184 * t101; m(7) * (t22 * t143 - t145 * t23); m(7) * (t19 * t8 + t22 * t32 + t23 * t31) + t147; m(7) * (-t19 * t132 + (t143 * t23 + t145 * t22) * t131); m(7) * (t15 * t19 + t22 * t35 + t23 * t34) - t147; t101 * t2 + t99 * t1 + t105 * (t14 * t101 + t13 * t99 + t18) + m(7) * (t19 ^ 2 + t22 ^ 2 + t23 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
