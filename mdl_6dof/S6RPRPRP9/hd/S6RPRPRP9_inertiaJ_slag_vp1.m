% Calculate joint inertia matrix for
% S6RPRPRP9
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
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:27:35
% EndTime: 2019-03-09 03:27:41
% DurationCPUTime: 2.40s
% Computational Cost: add. (4171->372), mult. (6009->551), div. (0->0), fcn. (6369->8), ass. (0->175)
t153 = pkin(9) + qJ(5);
t146 = sin(t153);
t147 = cos(t153);
t161 = sin(qJ(1));
t202 = t161 * t147;
t160 = sin(qJ(3));
t163 = cos(qJ(1));
t204 = t160 * t163;
t103 = t146 * t204 + t202;
t197 = t163 * t147;
t203 = t161 * t146;
t105 = -t160 * t197 + t203;
t162 = cos(qJ(3));
t198 = t162 * t163;
t237 = rSges(7,3) + qJ(6);
t239 = rSges(7,1) + pkin(5);
t217 = rSges(7,2) * t198 - t103 * t237 + t105 * t239;
t240 = -t161 / 0.2e1;
t221 = t163 / 0.2e1;
t238 = rSges(5,3) + qJ(4);
t208 = t146 * t162;
t84 = Icges(7,6) * t160 + (Icges(7,5) * t147 + Icges(7,3) * t146) * t162;
t85 = Icges(6,3) * t160 + (Icges(6,5) * t147 - Icges(6,6) * t146) * t162;
t86 = Icges(7,2) * t160 + (Icges(7,4) * t147 + Icges(7,6) * t146) * t162;
t88 = Icges(7,4) * t160 + (Icges(7,1) * t147 + Icges(7,5) * t146) * t162;
t89 = Icges(6,5) * t160 + (Icges(6,1) * t147 - Icges(6,4) * t146) * t162;
t236 = t84 * t208 + (t88 + t89) * t147 * t162 + (t85 + t86) * t160;
t101 = t160 * t203 - t197;
t102 = t146 * t163 + t160 * t202;
t199 = t161 * t162;
t46 = Icges(7,5) * t102 - Icges(7,6) * t199 + Icges(7,3) * t101;
t50 = Icges(7,4) * t102 - Icges(7,2) * t199 + Icges(7,6) * t101;
t54 = Icges(7,1) * t102 - Icges(7,4) * t199 + Icges(7,5) * t101;
t12 = t101 * t46 + t102 * t54 - t50 * t199;
t47 = Icges(7,5) * t105 + Icges(7,6) * t198 - Icges(7,3) * t103;
t51 = Icges(7,4) * t105 + Icges(7,2) * t198 - Icges(7,6) * t103;
t55 = Icges(7,1) * t105 + Icges(7,4) * t198 - Icges(7,5) * t103;
t13 = t101 * t47 + t102 * t55 - t51 * t199;
t48 = Icges(6,5) * t102 - Icges(6,6) * t101 - Icges(6,3) * t199;
t52 = Icges(6,4) * t102 - Icges(6,2) * t101 - Icges(6,6) * t199;
t56 = Icges(6,1) * t102 - Icges(6,4) * t101 - Icges(6,5) * t199;
t14 = -t101 * t52 + t102 * t56 - t48 * t199;
t49 = Icges(6,5) * t105 + Icges(6,6) * t103 + Icges(6,3) * t198;
t53 = Icges(6,4) * t105 + Icges(6,2) * t103 + Icges(6,6) * t198;
t57 = Icges(6,1) * t105 + Icges(6,4) * t103 + Icges(6,5) * t198;
t15 = -t101 * t53 + t102 * t57 - t49 * t199;
t29 = t101 * t84 + t102 * t88 - t86 * t199;
t87 = Icges(6,6) * t160 + (Icges(6,4) * t147 - Icges(6,2) * t146) * t162;
t30 = -t101 * t87 + t102 * t89 - t85 * t199;
t235 = ((t13 + t15) * t163 + (-t12 - t14) * t161) * t162 + (t29 + t30) * t160;
t16 = -t103 * t46 + t105 * t54 + t50 * t198;
t17 = -t103 * t47 + t105 * t55 + t51 * t198;
t18 = t103 * t52 + t105 * t56 + t48 * t198;
t19 = t103 * t53 + t105 * t57 + t49 * t198;
t31 = -t103 * t84 + t105 * t88 + t86 * t198;
t32 = t103 * t87 + t105 * t89 + t85 * t198;
t234 = ((t17 + t19) * t163 + (-t16 - t18) * t161) * t162 + (t31 + t32) * t160;
t144 = pkin(3) * t204;
t149 = t163 * qJ(2);
t157 = sin(pkin(9));
t158 = cos(pkin(9));
t200 = t161 * t158;
t121 = t157 * t204 + t200;
t201 = t161 * t157;
t122 = -t158 * t204 + t201;
t182 = -t122 * rSges(5,1) - t121 * rSges(5,2);
t187 = t238 * t162;
t224 = -pkin(1) - pkin(7);
t42 = t224 * t161 - t163 * t187 + t144 + t149 + t182;
t194 = t163 * pkin(1) + t161 * qJ(2);
t188 = t163 * pkin(7) + t194;
t119 = t158 * t163 - t160 * t201;
t206 = t157 * t163;
t120 = t160 * t200 + t206;
t205 = t160 * t161;
t227 = -t120 * rSges(5,1) - t119 * rSges(5,2) - pkin(3) * t205;
t43 = -t161 * t187 + t188 - t227;
t233 = m(5) * (t161 * t42 - t163 * t43);
t145 = pkin(4) * t158 + pkin(3);
t159 = -pkin(8) - qJ(4);
t195 = t145 * t204 + t159 * t198;
t164 = t149 + (-pkin(4) * t157 + t224) * t161 + t195;
t61 = t105 * rSges(6,1) + t103 * rSges(6,2) + rSges(6,3) * t198;
t36 = t164 - t61;
t190 = pkin(4) * t206 + t145 * t205 + t159 * t199;
t165 = t188 + t190;
t59 = t102 * rSges(6,1) - t101 * rSges(6,2) - rSges(6,3) * t199;
t37 = t165 + t59;
t232 = m(6) * (t161 * t36 - t163 * t37);
t26 = t164 - t217;
t218 = -rSges(7,2) * t199 + t101 * t237 + t102 * t239;
t27 = t165 + t218;
t231 = m(7) * (t161 * t26 - t163 * t27);
t20 = t160 * t50 + (t146 * t46 + t147 * t54) * t162;
t22 = t160 * t48 + (-t146 * t52 + t147 * t56) * t162;
t230 = t20 + t22;
t21 = t160 * t51 + (t146 * t47 + t147 * t55) * t162;
t23 = t160 * t49 + (-t146 * t53 + t147 * t57) * t162;
t229 = t21 + t23;
t228 = (rSges(4,1) * t160 + rSges(4,2) * t162) * t163;
t154 = t161 ^ 2;
t156 = t163 ^ 2;
t98 = Icges(5,6) * t160 + (Icges(5,4) * t158 - Icges(5,2) * t157) * t162;
t226 = t98 / 0.2e1;
t99 = Icges(5,5) * t160 + (Icges(5,1) * t158 - Icges(5,4) * t157) * t162;
t225 = t99 / 0.2e1;
t222 = t161 / 0.2e1;
t134 = rSges(4,1) * t162 - rSges(4,2) * t160;
t220 = m(4) * t134;
t219 = (-t87 * t208 + t236) * t160;
t215 = rSges(7,2) * t160 + (t146 * t237 + t147 * t239) * t162;
t184 = qJ(4) * t198 - t144;
t115 = t163 * t184;
t214 = t115 + t163 * (pkin(4) * t201 - t184 - t195);
t133 = pkin(3) * t162 + qJ(4) * t160;
t124 = t161 * t133;
t83 = (-pkin(3) + t145) * t162 + (-qJ(4) - t159) * t160;
t212 = t161 * t83 + t124;
t211 = -t133 - t83;
t210 = Icges(4,4) * t160;
t209 = Icges(4,4) * t162;
t193 = t154 + t156;
t189 = rSges(4,1) * t205 + rSges(4,2) * t199 + t163 * rSges(4,3);
t186 = t215 * t162;
t185 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t24 = t218 * t160 + t161 * t186;
t25 = -t217 * t160 + t163 * t186;
t181 = t25 * t161 - t163 * t24;
t38 = t215 * t161 + t212;
t39 = (t211 - t215) * t163;
t178 = t38 * t161 - t163 * t39;
t91 = rSges(6,3) * t160 + (rSges(6,1) * t147 - rSges(6,2) * t146) * t162;
t40 = t160 * t59 + t91 * t199;
t41 = -t160 * t61 + t91 * t198;
t177 = t41 * t161 - t163 * t40;
t44 = t161 * t91 + t212;
t45 = (-t91 + t211) * t163;
t175 = t44 * t161 - t163 * t45;
t100 = rSges(5,3) * t160 + (rSges(5,1) * t158 - rSges(5,2) * t157) * t162;
t73 = t100 * t161 + t124;
t74 = (-t100 - t133) * t163;
t174 = t73 * t161 - t163 * t74;
t173 = Icges(4,1) * t160 + t209;
t172 = Icges(4,2) * t162 + t210;
t171 = Icges(4,5) * t160 + Icges(4,6) * t162;
t170 = t101 * t161 + t103 * t163;
t167 = -t30 / 0.2e1 - t29 / 0.2e1 - t22 / 0.2e1 - t20 / 0.2e1;
t166 = t32 / 0.2e1 + t31 / 0.2e1 + t23 / 0.2e1 + t21 / 0.2e1;
t155 = t162 ^ 2;
t135 = rSges(2,1) * t163 - t161 * rSges(2,2);
t132 = -t161 * rSges(2,1) - rSges(2,2) * t163;
t129 = Icges(4,5) * t162 - Icges(4,6) * t160;
t114 = -rSges(3,2) * t163 + t161 * rSges(3,3) + t194;
t113 = rSges(3,3) * t163 + t149 + (rSges(3,2) - pkin(1)) * t161;
t107 = Icges(4,3) * t161 - t171 * t163;
t106 = Icges(4,3) * t163 + t171 * t161;
t79 = t188 + t189;
t78 = t149 + t228 + (-rSges(4,3) + t224) * t161;
t71 = Icges(5,1) * t122 + Icges(5,4) * t121 + Icges(5,5) * t198;
t70 = Icges(5,1) * t120 + Icges(5,4) * t119 - Icges(5,5) * t199;
t69 = Icges(5,4) * t122 + Icges(5,2) * t121 + Icges(5,6) * t198;
t68 = Icges(5,4) * t120 + Icges(5,2) * t119 - Icges(5,6) * t199;
t67 = Icges(5,5) * t122 + Icges(5,6) * t121 + Icges(5,3) * t198;
t66 = Icges(5,5) * t120 + Icges(5,6) * t119 - Icges(5,3) * t199;
t63 = -t161 * t189 + (t161 * rSges(4,3) - t228) * t163;
t33 = (-t161 * t61 - t163 * t59) * t162;
t28 = t115 + t163 * (rSges(5,3) * t198 - t182) + (t199 * t238 + t227) * t161;
t11 = (-t217 * t161 - t218 * t163) * t162;
t10 = t163 * t61 + (-t59 - t190) * t161 + t214;
t9 = t217 * t163 + (-t190 - t218) * t161 + t214;
t8 = t19 * t161 + t163 * t18;
t7 = t16 * t163 + t17 * t161;
t6 = t14 * t163 + t15 * t161;
t5 = t12 * t163 + t13 * t161;
t1 = [Icges(3,1) + Icges(2,3) + ((Icges(5,5) * t158 - Icges(5,6) * t157) * t162 - t209 + (Icges(5,3) + Icges(4,2)) * t160) * t160 + (Icges(4,1) * t162 - t146 * t87 - t157 * t98 + t158 * t99 - t210) * t162 + m(7) * (t26 ^ 2 + t27 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t132 ^ 2 + t135 ^ 2) + t236; t231 + t232 + t233 + m(4) * (t161 * t78 - t163 * t79) + m(3) * (t161 * t113 - t114 * t163); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t185) * t193; m(7) * (t26 * t38 + t27 * t39) + m(6) * (t36 * t44 + t37 * t45) + m(5) * (t42 * t73 + t43 * t74) + (t119 * t226 + t120 * t225 - t79 * t220 + t129 * t221 + (t66 / 0.2e1 - Icges(4,6) * t163 / 0.2e1 + t172 * t240) * t160 - t167) * t163 + (t121 * t226 + t122 * t225 + t78 * t220 + t129 * t222 + (t67 / 0.2e1 + Icges(4,6) * t240 + t172 * t221) * t160 + t166) * t161 + ((Icges(4,5) * t161 - t157 * t69 + t158 * t71 - t173 * t163) * t222 + (Icges(4,5) * t163 - t157 * t68 + t158 * t70 + t173 * t161) * t221) * t162; m(5) * t174 + m(6) * t175 + m(7) * t178 + t193 * t220; m(7) * (t38 ^ 2 + t39 ^ 2 + t9 ^ 2) + m(6) * (t10 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t28 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t193 * t134 ^ 2 + t63 ^ 2) + (t156 * t106 + t5 + t6 + (t119 * t68 + t120 * t70 - t66 * t199) * t163) * t163 + (t7 + t8 + t154 * t107 + (t121 * t69 + t122 * t71 + t67 * t198) * t161 + (t161 * t106 + t163 * t107 + t119 * t69 + t120 * t71 + t121 * t68 + t122 * t70 + t66 * t198 - t67 * t199) * t163) * t161; 0.2e1 * (-t231 / 0.2e1 - t232 / 0.2e1 - t233 / 0.2e1) * t162; -0.2e1 * t185 * t193 * t162; m(7) * (t160 * t9 - t178 * t162) + m(6) * (t160 * t10 - t175 * t162) + m(5) * (t160 * t28 - t174 * t162); 0.2e1 * t185 * (t193 * t155 + t160 ^ 2); m(7) * (t24 * t27 + t25 * t26) + m(6) * (t36 * t41 + t37 * t40) + (t167 * t161 + t166 * t163) * t162 + t219; m(6) * t177 + m(7) * t181; m(7) * (t11 * t9 + t24 * t39 + t25 * t38) + m(6) * (t10 * t33 + t40 * t45 + t41 * t44) + ((t8 / 0.2e1 + t7 / 0.2e1) * t163 + (-t6 / 0.2e1 - t5 / 0.2e1) * t161) * t162 + (t161 * t229 + t163 * t230) * t160 / 0.2e1 + t234 * t222 + t235 * t221; m(6) * (t33 * t160 - t177 * t162) + m(7) * (t11 * t160 - t181 * t162); t219 * t160 + m(7) * (t11 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t33 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t234 * t163 - t235 * t161 + (-t161 * t230 + t163 * t229) * t160) * t162; m(7) * (t101 * t26 - t103 * t27); m(7) * t170; m(7) * (t101 * t38 - t103 * t39 + t9 * t208); m(7) * (t146 * t160 - t170) * t162; m(7) * (t101 * t25 - t103 * t24 + t11 * t208); m(7) * (t146 ^ 2 * t155 + t101 ^ 2 + t103 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
