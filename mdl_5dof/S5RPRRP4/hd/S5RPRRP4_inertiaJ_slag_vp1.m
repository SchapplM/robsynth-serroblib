% Calculate joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:48
% EndTime: 2022-01-23 09:31:51
% DurationCPUTime: 2.74s
% Computational Cost: add. (4547->283), mult. (6061->394), div. (0->0), fcn. (6444->8), ass. (0->137)
t224 = Icges(5,1) + Icges(6,1);
t223 = Icges(5,4) + Icges(6,4);
t222 = Icges(5,5) + Icges(6,5);
t221 = Icges(5,2) + Icges(6,2);
t220 = Icges(5,6) + Icges(6,6);
t219 = Icges(6,3) + Icges(5,3);
t149 = qJ(3) + qJ(4);
t139 = sin(t149);
t140 = cos(t149);
t150 = sin(pkin(8));
t151 = cos(pkin(8));
t213 = t220 * t151 + (t221 * t139 - t223 * t140) * t150;
t218 = t213 * t139;
t217 = -t222 * t151 + (-t223 * t139 + t224 * t140) * t150;
t155 = cos(qJ(1));
t153 = sin(qJ(1));
t182 = t151 * t153;
t115 = -t139 * t182 - t140 * t155;
t116 = -t139 * t155 + t140 * t182;
t177 = t153 * t150;
t203 = t220 * t115 + t222 * t116 + t219 * t177;
t181 = t151 * t155;
t117 = -t139 * t181 + t140 * t153;
t118 = t139 * t153 + t140 * t181;
t174 = t155 * t150;
t202 = t220 * t117 + t222 * t118 + t219 * t174;
t201 = t221 * t115 + t223 * t116 + t220 * t177;
t200 = t221 * t117 + t223 * t118 + t220 * t174;
t199 = t223 * t115 + t224 * t116 + t222 * t177;
t198 = t223 * t117 + t224 * t118 + t222 * t174;
t214 = t219 * t151 + (t220 * t139 - t222 * t140) * t150;
t215 = t217 * t140 * t150;
t209 = t150 * t218 + t214 * t151 + t215;
t216 = t209 * t151;
t156 = pkin(7) + pkin(6);
t145 = -qJ(5) - t156;
t208 = rSges(6,3) - t145;
t211 = t203 * t151 + (t201 * t139 - t199 * t140) * t150;
t210 = t202 * t151 + (t200 * t139 - t198 * t140) * t150;
t142 = t155 * qJ(2);
t207 = -rSges(6,1) * t116 - rSges(6,2) * t115 + t142;
t206 = t213 * t115 - t116 * t217 + t214 * t177;
t205 = t213 * t117 - t118 * t217 + t214 * t174;
t193 = pkin(2) * t151 + pkin(1);
t154 = cos(qJ(3));
t194 = pkin(3) * t154;
t114 = t150 * t156 + t151 * t194 + t193;
t152 = sin(qJ(3));
t195 = pkin(3) * t152;
t129 = pkin(4) * t139 + t195;
t138 = qJ(2) + t195;
t168 = -pkin(2) - t194;
t128 = pkin(4) * t140 - t168;
t165 = -t128 * t151 - pkin(1);
t204 = (-t129 + t138) * t155 + (-t145 * t150 - t114 - t165) * t153 + rSges(6,3) * t177 - t207;
t197 = (t156 - t208) * t151 + (rSges(6,1) * t140 - rSges(6,2) * t139 + t128 + t168) * t150;
t196 = t216 + (t211 * t153 + t210 * t155) * t150;
t141 = t153 * qJ(2);
t172 = t155 * pkin(1) + t141;
t36 = t118 * rSges(6,1) + t117 * rSges(6,2) + t128 * t181 + t153 * t129 + t208 * t174 + t172;
t191 = t204 * t174;
t173 = t114 * t155 + t138 * t153;
t190 = t173 - t36;
t183 = t150 * t154;
t119 = pkin(3) * t183 + (pkin(6) - t156) * t151;
t130 = pkin(6) * t150 + t193;
t186 = t138 * t155;
t61 = -t186 + t142 + (t114 - t130) * t153;
t189 = t119 * t177 + t151 * t61;
t164 = -t130 * t155 - t141;
t62 = t164 + t173;
t78 = t118 * rSges(5,1) + t117 * rSges(5,2) + rSges(5,3) * t174;
t188 = -t62 - t78;
t104 = -rSges(5,3) * t151 + (rSges(5,1) * t140 - rSges(5,2) * t139) * t150;
t161 = -rSges(5,1) * t116 - rSges(5,2) * t115;
t76 = rSges(5,3) * t177 - t161;
t44 = t104 * t177 + t151 * t76;
t111 = -Icges(4,6) * t151 + (Icges(4,4) * t154 - Icges(4,2) * t152) * t150;
t180 = t152 * t111;
t179 = t152 * t153;
t178 = t152 * t155;
t176 = t153 * t154;
t175 = t154 * t155;
t171 = t153 ^ 2 + t155 ^ 2;
t170 = -t62 + t190;
t169 = ((t201 * t115 + t199 * t116 + t203 * t177) * t177 + t206 * t151) * t177 + ((t200 * t117 + t198 * t118 + t202 * t174) * t174 + t205 * t151 + (t200 * t115 + t198 * t116 + t201 * t117 + t199 * t118 + t203 * t174 + t202 * t177) * t177) * t174;
t125 = -t151 * t178 + t176;
t126 = t151 * t175 + t179;
t87 = t126 * rSges(4,1) + t125 * rSges(4,2) + rSges(4,3) * t174;
t12 = t204 * t151 + t197 * t177;
t163 = rSges(3,1) * t151 - rSges(3,2) * t150;
t123 = -t151 * t179 - t175;
t124 = t151 * t176 - t178;
t162 = -rSges(4,1) * t124 - rSges(4,2) * t123;
t159 = t196 * t151 + t169;
t158 = (-t206 - t211) * t177 / 0.2e1 + (-t205 - t210) * t174 / 0.2e1;
t132 = rSges(2,1) * t155 - rSges(2,2) * t153;
t131 = -rSges(2,1) * t153 - rSges(2,2) * t155;
t113 = -rSges(4,3) * t151 + (rSges(4,1) * t154 - rSges(4,2) * t152) * t150;
t112 = -Icges(4,5) * t151 + (Icges(4,1) * t154 - Icges(4,4) * t152) * t150;
t110 = -Icges(4,3) * t151 + (Icges(4,5) * t154 - Icges(4,6) * t152) * t150;
t95 = t112 * t183;
t94 = rSges(3,3) * t153 + t163 * t155 + t172;
t93 = rSges(3,3) * t155 + t142 + (-pkin(1) - t163) * t153;
t86 = rSges(4,3) * t177 - t162;
t85 = Icges(4,1) * t126 + Icges(4,4) * t125 + Icges(4,5) * t174;
t84 = Icges(4,1) * t124 + Icges(4,4) * t123 + Icges(4,5) * t177;
t83 = Icges(4,4) * t126 + Icges(4,2) * t125 + Icges(4,6) * t174;
t82 = Icges(4,4) * t124 + Icges(4,2) * t123 + Icges(4,6) * t177;
t81 = Icges(4,5) * t126 + Icges(4,6) * t125 + Icges(4,3) * t174;
t80 = Icges(4,5) * t124 + Icges(4,6) * t123 + Icges(4,3) * t177;
t57 = t76 * t174;
t55 = t61 * t174;
t54 = t142 + (-rSges(4,3) * t150 - t130) * t153 + t162;
t53 = -t164 + t87;
t52 = -t113 * t174 - t151 * t87;
t51 = t113 * t177 + t151 * t86;
t50 = -t151 * t110 - t150 * t180 + t95;
t49 = t78 + t173;
t48 = t186 + (-rSges(5,3) * t150 - t114) * t153 + t161;
t45 = -t104 * t174 - t151 * t78;
t39 = (-t153 * t87 + t155 * t86) * t150;
t38 = t110 * t174 + t111 * t125 + t112 * t126;
t37 = t110 * t177 + t111 * t123 + t112 * t124;
t35 = t129 * t155 + (-t208 * t150 + t165) * t153 + t207;
t34 = -t78 * t177 + t57;
t25 = -t151 * t81 + (-t152 * t83 + t154 * t85) * t150;
t24 = -t151 * t80 + (-t152 * t82 + t154 * t84) * t150;
t23 = t188 * t151 + (-t104 - t119) * t174;
t22 = t44 + t189;
t13 = t190 * t151 - t174 * t197;
t11 = t188 * t177 + t55 + t57;
t10 = t190 * t177 + t191;
t9 = t170 * t151 + (-t119 - t197) * t174;
t8 = t12 + t189;
t7 = t170 * t177 + t191 + t55;
t1 = [Icges(2,3) + t95 + (Icges(3,2) * t151 - t110 + t214) * t151 + (Icges(3,1) * t150 + 0.2e1 * Icges(3,4) * t151 - t180 + t218) * t150 + m(6) * (t35 ^ 2 + t36 ^ 2) + m(5) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2) + m(3) * (t93 ^ 2 + t94 ^ 2) + m(2) * (t131 ^ 2 + t132 ^ 2) + t215; m(6) * (t153 * t35 - t155 * t36) + m(5) * (t153 * t48 - t155 * t49) + m(4) * (t153 * t54 - t155 * t53) + m(3) * (t153 * t93 - t155 * t94); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t171; (-t50 - t209) * t151 + m(6) * (t35 * t8 + t36 * t9) + m(5) * (t22 * t48 + t23 * t49) + m(4) * (t51 * t54 + t52 * t53) + ((t38 / 0.2e1 + t25 / 0.2e1) * t155 + (t24 / 0.2e1 + t37 / 0.2e1) * t153) * t150 + t158; m(4) * (t153 * t51 - t155 * t52) + m(5) * (t153 * t22 - t155 * t23) + m(6) * (t153 * t8 - t155 * t9); ((t123 * t83 + t124 * t85 + t81 * t177) * t174 + (t123 * t82 + t124 * t84 + t80 * t177) * t177) * t177 + ((t125 * t83 + t126 * t85 + t81 * t174) * t174 + (t125 * t82 + t126 * t84 + t80 * t174) * t177) * t174 + m(6) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(4) * (t39 ^ 2 + t51 ^ 2 + t52 ^ 2) + t169 + (-(t153 * t24 + t155 * t25) * t150 - t37 * t177 - t38 * t174 + t50 * t151 + t196) * t151; -t216 + m(6) * (t12 * t35 + t13 * t36) + m(5) * (t44 * t48 + t45 * t49) + t158; m(5) * (t153 * t44 - t155 * t45) + m(6) * (t12 * t153 - t13 * t155); m(6) * (t10 * t7 + t12 * t8 + t13 * t9) + m(5) * (t11 * t34 + t22 * t44 + t23 * t45) + t159; m(5) * (t34 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t10 ^ 2 + t12 ^ 2 + t13 ^ 2) + t159; m(6) * (t153 * t36 + t155 * t35) * t150; 0; m(6) * (-t151 * t7 + (t153 * t9 + t155 * t8) * t150); m(6) * (-t10 * t151 + (t12 * t155 + t13 * t153) * t150); m(6) * (t171 * t150 ^ 2 + t151 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
