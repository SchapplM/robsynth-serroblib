% Calculate joint inertia matrix for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:50
% EndTime: 2019-03-09 03:07:54
% DurationCPUTime: 2.38s
% Computational Cost: add. (6404->367), mult. (5810->543), div. (0->0), fcn. (6210->10), ass. (0->173)
t226 = rSges(7,3) + qJ(6);
t227 = rSges(7,1) + pkin(5);
t151 = pkin(10) + qJ(5);
t146 = sin(t151);
t148 = cos(t151);
t152 = qJ(1) + pkin(9);
t149 = cos(t152);
t147 = sin(t152);
t159 = cos(qJ(3));
t195 = t147 * t159;
t94 = t146 * t195 + t149 * t148;
t95 = -t149 * t146 + t148 * t195;
t229 = -t226 * t94 - t227 * t95;
t215 = t147 / 0.2e1;
t228 = t149 / 0.2e1;
t157 = sin(qJ(3));
t100 = -Icges(7,2) * t159 + (Icges(7,4) * t148 + Icges(7,6) * t146) * t157;
t99 = -Icges(6,3) * t159 + (Icges(6,5) * t148 - Icges(6,6) * t146) * t157;
t225 = -t100 - t99;
t196 = t147 * t157;
t47 = Icges(7,5) * t95 + Icges(7,6) * t196 + Icges(7,3) * t94;
t51 = Icges(7,4) * t95 + Icges(7,2) * t196 + Icges(7,6) * t94;
t55 = Icges(7,1) * t95 + Icges(7,4) * t196 + Icges(7,5) * t94;
t12 = t196 * t51 + t94 * t47 + t95 * t55;
t192 = t149 * t157;
t191 = t149 * t159;
t96 = t146 * t191 - t147 * t148;
t97 = t147 * t146 + t148 * t191;
t48 = Icges(7,5) * t97 + Icges(7,6) * t192 + Icges(7,3) * t96;
t52 = Icges(7,4) * t97 + Icges(7,2) * t192 + Icges(7,6) * t96;
t56 = Icges(7,1) * t97 + Icges(7,4) * t192 + Icges(7,5) * t96;
t13 = t196 * t52 + t94 * t48 + t95 * t56;
t49 = Icges(6,5) * t95 - Icges(6,6) * t94 + Icges(6,3) * t196;
t53 = Icges(6,4) * t95 - Icges(6,2) * t94 + Icges(6,6) * t196;
t57 = Icges(6,1) * t95 - Icges(6,4) * t94 + Icges(6,5) * t196;
t14 = t196 * t49 - t94 * t53 + t95 * t57;
t50 = Icges(6,5) * t97 - Icges(6,6) * t96 + Icges(6,3) * t192;
t54 = Icges(6,4) * t97 - Icges(6,2) * t96 + Icges(6,6) * t192;
t58 = Icges(6,1) * t97 - Icges(6,4) * t96 + Icges(6,5) * t192;
t15 = t196 * t50 - t94 * t54 + t95 * t58;
t102 = -Icges(7,4) * t159 + (Icges(7,1) * t148 + Icges(7,5) * t146) * t157;
t98 = -Icges(7,6) * t159 + (Icges(7,5) * t148 + Icges(7,3) * t146) * t157;
t30 = t100 * t196 + t95 * t102 + t94 * t98;
t101 = -Icges(6,6) * t159 + (Icges(6,4) * t148 - Icges(6,2) * t146) * t157;
t103 = -Icges(6,5) * t159 + (Icges(6,1) * t148 - Icges(6,4) * t146) * t157;
t31 = -t94 * t101 + t95 * t103 + t196 * t99;
t224 = (-t30 - t31) * t159 + ((t13 + t15) * t149 + (t12 + t14) * t147) * t157;
t16 = t192 * t51 + t96 * t47 + t97 * t55;
t17 = t192 * t52 + t96 * t48 + t97 * t56;
t18 = t192 * t49 - t96 * t53 + t97 * t57;
t19 = t192 * t50 - t96 * t54 + t97 * t58;
t32 = t100 * t192 + t97 * t102 + t96 * t98;
t33 = -t96 * t101 + t97 * t103 + t192 * t99;
t223 = (-t32 - t33) * t159 + ((t17 + t19) * t149 + (t16 + t18) * t147) * t157;
t20 = -t159 * t51 + (t146 * t47 + t148 * t55) * t157;
t22 = -t159 * t49 + (-t146 * t53 + t148 * t57) * t157;
t222 = -t20 - t22;
t21 = -t159 * t52 + (t146 * t48 + t148 * t56) * t157;
t23 = -t159 * t50 + (-t146 * t54 + t148 * t58) * t157;
t221 = t21 + t23;
t198 = t146 * t157;
t220 = t98 * t198 + (t102 + t103) * t148 * t157;
t144 = t147 ^ 2;
t145 = t149 ^ 2;
t219 = t159 ^ 2;
t218 = m(5) / 0.2e1;
t217 = m(6) / 0.2e1;
t216 = m(7) / 0.2e1;
t214 = -t149 / 0.2e1;
t128 = t157 * rSges(4,1) + t159 * rSges(4,2);
t212 = m(4) * t128;
t211 = pkin(3) * t159;
t158 = sin(qJ(1));
t210 = t158 * pkin(1);
t155 = cos(pkin(10));
t143 = t155 * pkin(4) + pkin(3);
t209 = -pkin(3) + t143;
t208 = -t101 * t198 + t225 * t159 + t220;
t207 = rSges(7,2) * t196 - t229;
t206 = rSges(7,2) * t192 + t226 * t96 + t227 * t97;
t184 = pkin(3) * t191 + qJ(4) * t192;
t199 = qJ(4) * t157;
t204 = t144 * (t199 + t211) + t149 * t184;
t203 = t149 * rSges(4,3);
t127 = t157 * pkin(3) - t159 * qJ(4);
t156 = -pkin(8) - qJ(4);
t202 = -t127 - (qJ(4) + t156) * t159 - t209 * t157;
t201 = Icges(4,4) * t157;
t200 = Icges(4,4) * t159;
t154 = sin(pkin(10));
t197 = t147 * t154;
t193 = t149 * t154;
t190 = t154 * t159;
t189 = t155 * t159;
t188 = t156 * t157;
t187 = -t159 * rSges(7,2) + (t226 * t146 + t227 * t148) * t157;
t186 = t159 * rSges(5,3) - (rSges(5,1) * t155 - rSges(5,2) * t154) * t157 - t127;
t185 = -pkin(4) * t193 - t147 * t188;
t183 = t144 + t145;
t182 = -m(5) - m(6) - m(7);
t62 = t97 * rSges(6,1) - t96 * rSges(6,2) + rSges(6,3) * t192;
t105 = -t159 * rSges(6,3) + (rSges(6,1) * t148 - rSges(6,2) * t146) * t157;
t181 = -t105 + t202;
t111 = t147 * t155 - t149 * t190;
t112 = t149 * t189 + t197;
t180 = t112 * rSges(5,1) + t111 * rSges(5,2) + rSges(5,3) * t192;
t160 = cos(qJ(1));
t150 = t160 * pkin(1);
t179 = t149 * pkin(2) + t147 * pkin(7) + t150;
t178 = t149 * pkin(7) - t210;
t177 = -t143 * t159 - pkin(2);
t164 = pkin(4) * t197 + t143 * t191 - t149 * t188;
t176 = t147 * ((t159 * t209 - t199) * t147 + t185) + t149 * (t164 - t184) + t204;
t175 = -t187 + t202;
t174 = -t95 * rSges(6,1) + t94 * rSges(6,2);
t173 = rSges(4,1) * t159 - rSges(4,2) * t157;
t109 = -t147 * t190 - t149 * t155;
t110 = t147 * t189 - t193;
t172 = -t110 * rSges(5,1) - t109 * rSges(5,2);
t169 = Icges(4,1) * t159 - t201;
t168 = -Icges(4,2) * t157 + t200;
t167 = Icges(4,5) * t159 - Icges(4,6) * t157;
t166 = t178 - t185;
t165 = rSges(4,1) * t191 - rSges(4,2) * t192 + t147 * rSges(4,3);
t163 = t23 / 0.2e1 + t21 / 0.2e1 + t33 / 0.2e1 + t32 / 0.2e1;
t162 = t31 / 0.2e1 + t30 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t161 = t164 + t179;
t153 = t157 ^ 2;
t130 = t160 * rSges(2,1) - t158 * rSges(2,2);
t129 = -t158 * rSges(2,1) - t160 * rSges(2,2);
t124 = Icges(4,5) * t157 + Icges(4,6) * t159;
t119 = t149 * rSges(3,1) - t147 * rSges(3,2) + t150;
t118 = -t147 * rSges(3,1) - t149 * rSges(3,2) - t210;
t115 = -Icges(5,5) * t159 + (Icges(5,1) * t155 - Icges(5,4) * t154) * t157;
t114 = -Icges(5,6) * t159 + (Icges(5,4) * t155 - Icges(5,2) * t154) * t157;
t81 = Icges(4,3) * t147 + t149 * t167;
t80 = -Icges(4,3) * t149 + t147 * t167;
t76 = t186 * t149;
t75 = t186 * t147;
t74 = t165 + t179;
t73 = t203 + (-pkin(2) - t173) * t147 + t178;
t70 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t192;
t69 = Icges(5,1) * t110 + Icges(5,4) * t109 + Icges(5,5) * t196;
t68 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t192;
t67 = Icges(5,4) * t110 + Icges(5,2) * t109 + Icges(5,6) * t196;
t66 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t192;
t65 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t196;
t60 = rSges(6,3) * t196 - t174;
t46 = t149 * t165 + (t147 * t173 - t203) * t147;
t45 = t181 * t149;
t44 = t181 * t147;
t43 = t179 + t180 + t184;
t42 = (-t211 - pkin(2) + (-rSges(5,3) - qJ(4)) * t157) * t147 + t172 + t178;
t39 = -t105 * t192 - t159 * t62;
t38 = t105 * t196 + t159 * t60;
t37 = t175 * t149;
t36 = t175 * t147;
t35 = t161 + t62;
t34 = (-rSges(6,3) * t157 + t177) * t147 + t166 + t174;
t29 = (-t147 * t62 + t149 * t60) * t157;
t28 = t161 + t206;
t27 = (-rSges(7,2) * t157 + t177) * t147 + t166 + t229;
t26 = t147 * (rSges(5,3) * t196 - t172) + t149 * t180 + t204;
t25 = -t159 * t206 - t187 * t192;
t24 = t159 * t207 + t187 * t196;
t11 = (-t147 * t206 + t149 * t207) * t157;
t10 = t147 * t60 + t149 * t62 + t176;
t9 = t147 * t207 + t149 * t206 + t176;
t8 = t19 * t147 - t18 * t149;
t7 = t17 * t147 - t16 * t149;
t6 = -t14 * t149 + t15 * t147;
t5 = -t12 * t149 + t13 * t147;
t1 = [Icges(2,3) + Icges(3,3) + (-(Icges(5,5) * t155 - Icges(5,6) * t154) * t157 + t201 + (Icges(5,3) + Icges(4,2)) * t159 + t225) * t159 + (Icges(4,1) * t157 - t146 * t101 - t154 * t114 + t155 * t115 + t200) * t157 + m(6) * (t34 ^ 2 + t35 ^ 2) + m(7) * (t27 ^ 2 + t28 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t118 ^ 2 + t119 ^ 2) + m(2) * (t129 ^ 2 + t130 ^ 2) + t220; 0; m(3) + m(4) - t182; m(6) * (t45 * t34 + t44 * t35) + m(7) * (t27 * t37 + t36 * t28) + m(5) * (t76 * t42 + t75 * t43) + (-t73 * t212 - t109 * t114 / 0.2e1 - t110 * t115 / 0.2e1 + t124 * t228 + (t65 / 0.2e1 + Icges(4,6) * t228 - t147 * t168 / 0.2e1) * t159 - t162) * t149 + (-t74 * t212 + t111 * t114 / 0.2e1 + t112 * t115 / 0.2e1 + t124 * t215 + (Icges(4,6) * t215 + t168 * t228 - t66 / 0.2e1) * t159 + t163) * t147 + ((Icges(4,5) * t147 + t149 * t169 - t154 * t68 + t155 * t70) * t215 + (-Icges(4,5) * t149 + t147 * t169 - t154 * t67 + t155 * t69) * t214) * t157; m(4) * t46 + m(5) * t26 + m(6) * t10 + m(7) * t9; m(7) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t26 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t128 ^ 2 * t183 + t46 ^ 2) + (-t145 * t80 - t5 - t6 + (t109 * t67 + t110 * t69 + t65 * t196) * t149) * t149 + (t7 + t8 + t144 * t81 + (t111 * t68 + t112 * t70 + t66 * t192) * t147 + (-t109 * t68 - t110 * t70 - t111 * t67 - t112 * t69 - t147 * t80 + t149 * t81 - t192 * t65 - t196 * t66) * t149) * t147; 0.2e1 * ((t147 * t35 + t149 * t34) * t217 + (t147 * t28 + t149 * t27) * t216 + (t147 * t43 + t149 * t42) * t218) * t157; t182 * t159; m(7) * (-t159 * t9 + (t147 * t36 + t149 * t37) * t157) + m(6) * (-t159 * t10 + (t147 * t44 + t149 * t45) * t157) + m(5) * (-t159 * t26 + (t147 * t75 + t149 * t76) * t157); 0.2e1 * (t218 + t217 + t216) * (t153 * t183 + t219); -t208 * t159 + m(6) * (t38 * t34 + t39 * t35) + m(7) * (t24 * t27 + t25 * t28) + (t147 * t162 + t149 * t163) * t157; m(6) * t29 + m(7) * t11; m(7) * (t11 * t9 + t24 * t37 + t25 * t36) + m(6) * (t29 * t10 + t38 * t45 + t39 * t44) + ((t7 / 0.2e1 + t8 / 0.2e1) * t149 + (t5 / 0.2e1 + t6 / 0.2e1) * t147) * t157 + t223 * t215 + t224 * t214 - (t221 * t147 + t222 * t149) * t159 / 0.2e1; m(6) * (-t29 * t159 + (t147 * t39 + t149 * t38) * t157) + m(7) * (-t11 * t159 + (t147 * t25 + t149 * t24) * t157); m(7) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t29 ^ 2 + t38 ^ 2 + t39 ^ 2) + t208 * t219 + ((-t221 * t159 + t223) * t149 + (t222 * t159 + t224) * t147) * t157; m(7) * (t96 * t27 + t94 * t28); m(7) * t198; m(7) * (t198 * t9 + t94 * t36 + t96 * t37); m(7) * (-t146 * t159 + t147 * t94 + t149 * t96) * t157; m(7) * (t11 * t198 + t96 * t24 + t94 * t25); m(7) * (t153 * t146 ^ 2 + t94 ^ 2 + t96 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
