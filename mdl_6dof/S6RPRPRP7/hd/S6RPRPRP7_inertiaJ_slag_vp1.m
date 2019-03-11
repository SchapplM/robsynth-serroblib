% Calculate joint inertia matrix for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:13
% EndTime: 2019-03-09 03:21:18
% DurationCPUTime: 2.09s
% Computational Cost: add. (4099->334), mult. (5262->470), div. (0->0), fcn. (5467->8), ass. (0->155)
t150 = -qJ(6) - pkin(8);
t239 = -rSges(7,3) + t150;
t238 = Icges(4,3) + Icges(5,3);
t147 = qJ(3) + pkin(9);
t140 = sin(t147);
t141 = cos(t147);
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t237 = Icges(4,5) * t153 + Icges(5,5) * t140 + Icges(4,6) * t156 + Icges(5,6) * t141;
t236 = Icges(4,5) * t156 + Icges(5,5) * t141 - Icges(4,6) * t153 - Icges(5,6) * t140;
t155 = cos(qJ(5));
t152 = sin(qJ(5));
t72 = Icges(7,3) * t140 + (Icges(7,5) * t155 - Icges(7,6) * t152) * t141;
t73 = Icges(6,3) * t140 + (Icges(6,5) * t155 - Icges(6,6) * t152) * t141;
t76 = Icges(7,5) * t140 + (Icges(7,1) * t155 - Icges(7,4) * t152) * t141;
t77 = Icges(6,5) * t140 + (Icges(6,1) * t155 - Icges(6,4) * t152) * t141;
t235 = (t76 + t77) * t141 * t155 + (t72 + t73) * t140;
t74 = Icges(7,6) * t140 + (Icges(7,4) * t155 - Icges(7,2) * t152) * t141;
t75 = Icges(6,6) * t140 + (Icges(6,4) * t155 - Icges(6,2) * t152) * t141;
t234 = (-t74 - t75) * t152;
t157 = cos(qJ(1));
t191 = t155 * t157;
t154 = sin(qJ(1));
t194 = t154 * t152;
t104 = -t140 * t194 + t191;
t193 = t154 * t155;
t196 = t152 * t157;
t105 = t140 * t193 + t196;
t199 = t141 * t154;
t46 = Icges(7,5) * t105 + Icges(7,6) * t104 - Icges(7,3) * t199;
t50 = Icges(7,4) * t105 + Icges(7,2) * t104 - Icges(7,6) * t199;
t54 = Icges(7,1) * t105 + Icges(7,4) * t104 - Icges(7,5) * t199;
t11 = t104 * t50 + t105 * t54 - t46 * t199;
t106 = t140 * t196 + t193;
t107 = -t140 * t191 + t194;
t197 = t141 * t157;
t47 = Icges(7,5) * t107 + Icges(7,6) * t106 + Icges(7,3) * t197;
t51 = Icges(7,4) * t107 + Icges(7,2) * t106 + Icges(7,6) * t197;
t55 = Icges(7,1) * t107 + Icges(7,4) * t106 + Icges(7,5) * t197;
t12 = t104 * t51 + t105 * t55 - t47 * t199;
t48 = Icges(6,5) * t105 + Icges(6,6) * t104 - Icges(6,3) * t199;
t52 = Icges(6,4) * t105 + Icges(6,2) * t104 - Icges(6,6) * t199;
t56 = Icges(6,1) * t105 + Icges(6,4) * t104 - Icges(6,5) * t199;
t13 = t104 * t52 + t105 * t56 - t48 * t199;
t49 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t197;
t53 = Icges(6,4) * t107 + Icges(6,2) * t106 + Icges(6,6) * t197;
t57 = Icges(6,1) * t107 + Icges(6,4) * t106 + Icges(6,5) * t197;
t14 = t104 * t53 + t105 * t57 - t49 * t199;
t26 = t104 * t74 + t105 * t76 - t72 * t199;
t27 = t104 * t75 + t105 * t77 - t73 * t199;
t233 = ((t12 + t14) * t157 + (-t11 - t13) * t154) * t141 + (t26 + t27) * t140;
t15 = t106 * t50 + t107 * t54 + t46 * t197;
t16 = t106 * t51 + t107 * t55 + t47 * t197;
t17 = t106 * t52 + t107 * t56 + t48 * t197;
t18 = t106 * t53 + t107 * t57 + t49 * t197;
t28 = t106 * t74 + t107 * t76 + t72 * t197;
t29 = t106 * t75 + t107 * t77 + t73 * t197;
t232 = ((t16 + t18) * t157 + (-t15 - t17) * t154) * t141 + (t28 + t29) * t140;
t22 = t140 * t46 + (-t152 * t50 + t155 * t54) * t141;
t24 = t140 * t48 + (-t152 * t52 + t155 * t56) * t141;
t231 = t22 + t24;
t23 = t140 * t47 + (-t152 * t51 + t155 * t55) * t141;
t25 = t140 * t49 + (-t152 * t53 + t155 * t57) * t141;
t230 = t23 + t25;
t229 = t238 * t154 - t237 * t157;
t228 = t237 * t154 + t238 * t157;
t138 = pkin(5) * t155 + pkin(4);
t201 = t140 * t154;
t227 = t105 * rSges(7,1) + t104 * rSges(7,2) + pkin(5) * t196 + t138 * t201 + t239 * t199;
t226 = (rSges(4,1) * t153 + rSges(4,2) * t156) * t157;
t148 = t154 ^ 2;
t149 = t157 ^ 2;
t124 = rSges(4,1) * t156 - rSges(4,2) * t153;
t218 = m(4) * t124;
t217 = m(7) * t141;
t216 = pkin(3) * t156;
t215 = -pkin(8) - t150;
t214 = (t141 * t234 + t235) * t140;
t128 = pkin(4) * t201;
t100 = -pkin(8) * t199 + t128;
t213 = -t100 + t227;
t129 = t157 * t140 * pkin(4);
t174 = -t107 * rSges(7,1) - t106 * rSges(7,2);
t202 = t138 * t140;
t212 = pkin(5) * t194 + t129 + (t215 * t141 - t202) * t157 + rSges(7,3) * t197 - t174;
t209 = (rSges(7,1) * t155 - rSges(7,2) * t152 - pkin(4) + t138) * t141 + (t215 + rSges(7,3)) * t140;
t151 = -qJ(4) - pkin(7);
t188 = t157 * t153 * pkin(3) + t154 * t151;
t89 = t157 * (-t154 * pkin(7) - t188);
t208 = t157 * (pkin(8) * t197 - t129) + t89;
t207 = t105 * rSges(6,1) + t104 * rSges(6,2);
t144 = t157 * rSges(5,3);
t195 = t153 * t154;
t192 = t154 * t156;
t133 = pkin(3) * t195;
t103 = t133 + (-pkin(7) - t151) * t157;
t190 = -t100 - t103;
t114 = pkin(4) * t141 + pkin(8) * t140;
t134 = pkin(3) * t192;
t189 = t154 * t114 + t134;
t187 = t157 * pkin(1) + t154 * qJ(2);
t130 = t148 + t149;
t185 = (m(5) + m(6) + m(7)) * t130;
t184 = -rSges(5,1) * t201 - rSges(5,2) * t199 - t144;
t183 = rSges(4,1) * t195 + rSges(4,2) * t192 + t157 * rSges(4,3);
t143 = t157 * qJ(2);
t182 = t143 + t188;
t181 = t141 * (-rSges(6,3) - pkin(8));
t179 = t209 * t154;
t178 = -t114 - t216;
t176 = rSges(5,1) * t140 + rSges(5,2) * t141;
t175 = -t107 * rSges(6,1) - t106 * rSges(6,2);
t20 = t213 * t140 + t141 * t179;
t21 = -t212 * t140 + t209 * t197;
t169 = t21 * t154 - t157 * t20;
t33 = (-pkin(5) * t152 - pkin(1)) * t154 + (t239 * t141 + t202) * t157 + t174 + t182;
t160 = -t151 * t157 + t133 + t187;
t34 = t160 + t227;
t168 = t154 * t33 - t157 * t34;
t36 = t179 + t189;
t37 = (t178 - t209) * t157;
t167 = t36 * t154 - t157 * t37;
t159 = -t27 / 0.2e1 - t26 / 0.2e1 - t24 / 0.2e1 - t22 / 0.2e1;
t158 = t29 / 0.2e1 + t28 / 0.2e1 + t25 / 0.2e1 + t23 / 0.2e1;
t125 = rSges(2,1) * t157 - t154 * rSges(2,2);
t123 = -t154 * rSges(2,1) - rSges(2,2) * t157;
t112 = rSges(5,1) * t141 - rSges(5,2) * t140;
t102 = -rSges(3,2) * t157 + t154 * rSges(3,3) + t187;
t101 = rSges(3,3) * t157 + t143 + (rSges(3,2) - pkin(1)) * t154;
t81 = (-t112 - t216) * t157;
t80 = t112 * t154 + t134;
t79 = rSges(6,3) * t140 + (rSges(6,1) * t155 - rSges(6,2) * t152) * t141;
t69 = pkin(7) * t157 + t183 + t187;
t68 = t143 + t226 + (-rSges(4,3) - pkin(1) - pkin(7)) * t154;
t64 = t160 - t184;
t63 = t176 * t157 + (-rSges(5,3) - pkin(1)) * t154 + t182;
t62 = -t154 * t183 + (t154 * rSges(4,3) - t226) * t157;
t61 = rSges(6,3) * t197 - t175;
t59 = -rSges(6,3) * t199 + t207;
t43 = (t178 - t79) * t157;
t42 = t154 * t79 + t189;
t41 = t154 * t181 + t128 + t160 + t207;
t40 = -t154 * pkin(1) + t157 * t181 + t129 + t175 + t182;
t39 = -t140 * t61 + t79 * t197;
t38 = t140 * t59 + t79 * t199;
t35 = t89 - t176 * t149 + (-t103 + t184 + t144) * t154;
t30 = (-t154 * t61 - t157 * t59) * t141;
t19 = t157 * t61 + (-t59 + t190) * t154 + t208;
t10 = (-t212 * t154 - t213 * t157) * t141;
t9 = t212 * t157 + (t190 - t213) * t154 + t208;
t8 = t18 * t154 + t157 * t17;
t7 = t15 * t157 + t16 * t154;
t6 = t13 * t157 + t14 * t154;
t5 = t11 * t157 + t12 * t154;
t1 = [Icges(4,1) * t156 ^ 2 + Icges(3,1) + Icges(2,3) + (Icges(5,1) * t141 + t234) * t141 + m(7) * (t33 ^ 2 + t34 ^ 2) + m(6) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t63 ^ 2 + t64 ^ 2) + m(4) * (t68 ^ 2 + t69 ^ 2) + m(3) * (t101 ^ 2 + t102 ^ 2) + m(2) * (t123 ^ 2 + t125 ^ 2) + t235 + (-0.2e1 * Icges(4,4) * t156 + Icges(4,2) * t153) * t153 + (-0.2e1 * Icges(5,4) * t141 + Icges(5,2) * t140) * t140; m(7) * t168 + m(6) * (t154 * t40 - t157 * t41) + m(5) * (t154 * t63 - t157 * t64) + m(4) * (t154 * t68 - t157 * t69) + m(3) * (t154 * t101 - t102 * t157); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t130 + t185; m(7) * (t33 * t36 + t34 * t37) + m(6) * (t40 * t42 + t41 * t43) + m(5) * (t63 * t80 + t64 * t81) + (t236 * t157 - t69 * t218 - t159) * t157 + (t236 * t154 + t68 * t218 + t158) * t154; m(5) * (t80 * t154 - t157 * t81) + m(6) * (t42 * t154 - t157 * t43) + m(7) * t167 + t130 * t218; m(7) * (t36 ^ 2 + t37 ^ 2 + t9 ^ 2) + m(6) * (t19 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t35 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t130 * t124 ^ 2 + t62 ^ 2) + (t228 * t149 + t5 + t6) * t157 + (t7 + t8 + t229 * t148 + (t228 * t154 + t229 * t157) * t157) * t154; m(7) * (t154 * t34 + t157 * t33) + m(6) * (t154 * t41 + t157 * t40) + m(5) * (t154 * t64 + t157 * t63); 0; m(7) * (t154 * t37 + t157 * t36) + m(6) * (t154 * t43 + t157 * t42) + m(5) * (t154 * t81 + t157 * t80); t185; m(7) * (t20 * t34 + t21 * t33) + m(6) * (t38 * t41 + t39 * t40) + (t159 * t154 + t158 * t157) * t141 + t214; m(6) * (t39 * t154 - t157 * t38) + m(7) * t169; m(7) * (t10 * t9 + t20 * t37 + t21 * t36) + m(6) * (t19 * t30 + t38 * t43 + t39 * t42) + ((t7 / 0.2e1 + t8 / 0.2e1) * t157 + (-t6 / 0.2e1 - t5 / 0.2e1) * t154) * t141 + (t230 * t154 + t231 * t157) * t140 / 0.2e1 + t232 * t154 / 0.2e1 + t233 * t157 / 0.2e1; m(6) * (t38 * t154 + t157 * t39) + m(7) * (t20 * t154 + t157 * t21); t214 * t140 + m(7) * (t10 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t30 ^ 2 + t38 ^ 2 + t39 ^ 2) + (t232 * t157 - t233 * t154 + (-t231 * t154 + t230 * t157) * t140) * t141; -t168 * t217; -t130 * t217; m(7) * (t140 * t9 - t167 * t141); 0; m(7) * (t140 * t10 - t169 * t141); m(7) * (t130 * t141 ^ 2 + t140 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
