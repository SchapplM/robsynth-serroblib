% Calculate joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:09
% DurationCPUTime: 2.43s
% Computational Cost: add. (4412->355), mult. (5512->516), div. (0->0), fcn. (5725->8), ass. (0->162)
t145 = -qJ(6) - pkin(8);
t228 = rSges(7,3) - t145;
t227 = Icges(5,1) + Icges(4,3);
t140 = pkin(9) + qJ(3);
t134 = sin(t140);
t135 = cos(t140);
t226 = (-Icges(5,4) + Icges(4,5)) * t135 + (Icges(5,5) - Icges(4,6)) * t134;
t148 = sin(qJ(1));
t225 = -t148 / 0.2e1;
t224 = t148 / 0.2e1;
t150 = cos(qJ(1));
t223 = -t150 / 0.2e1;
t222 = t150 / 0.2e1;
t147 = sin(qJ(5));
t149 = cos(qJ(5));
t73 = Icges(7,3) * t134 + (-Icges(7,5) * t147 - Icges(7,6) * t149) * t135;
t74 = Icges(6,3) * t134 + (-Icges(6,5) * t147 - Icges(6,6) * t149) * t135;
t221 = (t73 + t74) * t134;
t75 = Icges(7,6) * t134 + (-Icges(7,4) * t147 - Icges(7,2) * t149) * t135;
t76 = Icges(6,6) * t134 + (-Icges(6,4) * t147 - Icges(6,2) * t149) * t135;
t77 = Icges(7,5) * t134 + (-Icges(7,1) * t147 - Icges(7,4) * t149) * t135;
t78 = Icges(6,5) * t134 + (-Icges(6,1) * t147 - Icges(6,4) * t149) * t135;
t220 = (-t75 - t76) * t149 + (-t77 - t78) * t147;
t186 = t149 * t150;
t188 = t148 * t147;
t101 = t134 * t186 - t188;
t189 = t147 * t150;
t181 = t134 * t189;
t187 = t148 * t149;
t102 = t181 + t187;
t190 = t135 * t150;
t47 = Icges(7,5) * t102 + Icges(7,6) * t101 + Icges(7,3) * t190;
t51 = Icges(7,4) * t102 + Icges(7,2) * t101 + Icges(7,6) * t190;
t55 = Icges(7,1) * t102 + Icges(7,4) * t101 + Icges(7,5) * t190;
t11 = t101 * t51 + t102 * t55 + t47 * t190;
t103 = t134 * t187 + t189;
t104 = t134 * t188 - t186;
t191 = t135 * t148;
t48 = Icges(7,5) * t104 + Icges(7,6) * t103 + Icges(7,3) * t191;
t52 = Icges(7,4) * t104 + Icges(7,2) * t103 + Icges(7,6) * t191;
t56 = Icges(7,1) * t104 + Icges(7,4) * t103 + Icges(7,5) * t191;
t12 = t101 * t52 + t102 * t56 + t48 * t190;
t49 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t190;
t53 = Icges(6,4) * t102 + Icges(6,2) * t101 + Icges(6,6) * t190;
t57 = Icges(6,1) * t102 + Icges(6,4) * t101 + Icges(6,5) * t190;
t13 = t101 * t53 + t102 * t57 + t49 * t190;
t50 = Icges(6,5) * t104 + Icges(6,6) * t103 + Icges(6,3) * t191;
t54 = Icges(6,4) * t104 + Icges(6,2) * t103 + Icges(6,6) * t191;
t58 = Icges(6,1) * t104 + Icges(6,4) * t103 + Icges(6,5) * t191;
t14 = t101 * t54 + t102 * t58 + t50 * t190;
t26 = t101 * t75 + t102 * t77 + t73 * t190;
t27 = t101 * t76 + t102 * t78 + t74 * t190;
t219 = ((t11 + t13) * t150 + (t12 + t14) * t148) * t135 + (t26 + t27) * t134;
t15 = t103 * t51 + t104 * t55 + t47 * t191;
t16 = t103 * t52 + t104 * t56 + t48 * t191;
t17 = t103 * t53 + t104 * t57 + t49 * t191;
t18 = t103 * t54 + t104 * t58 + t50 * t191;
t28 = t103 * t75 + t104 * t77 + t73 * t191;
t29 = t103 * t76 + t104 * t78 + t74 * t191;
t218 = ((t15 + t17) * t150 + (t16 + t18) * t148) * t135 + (t28 + t29) * t134;
t217 = t134 / 0.2e1;
t22 = t134 * t47 + (-t147 * t55 - t149 * t51) * t135;
t24 = t134 * t49 + (-t147 * t57 - t149 * t53) * t135;
t216 = t22 + t24;
t23 = t134 * t48 + (-t147 * t56 - t149 * t52) * t135;
t25 = t134 * t50 + (-t147 * t58 - t149 * t54) * t135;
t215 = t23 + t25;
t214 = -t226 * t148 + t227 * t150;
t213 = t227 * t148 + t226 * t150;
t131 = pkin(5) * t149 + pkin(4);
t212 = t102 * rSges(7,1) + t101 * rSges(7,2) + pkin(5) * t181 + t148 * t131 + t228 * t190;
t141 = t148 ^ 2;
t142 = t150 ^ 2;
t211 = m(5) / 0.2e1;
t210 = m(6) / 0.2e1;
t116 = rSges(4,1) * t134 + rSges(4,2) * t135;
t206 = m(4) * t116;
t205 = pkin(5) * t147;
t204 = -pkin(8) - t145;
t203 = (t220 * t135 + t221) * t134;
t183 = t148 * pkin(4) + pkin(8) * t190;
t202 = -t183 + t212;
t139 = t150 * pkin(4);
t172 = -t104 * rSges(7,1) - t103 * rSges(7,2);
t201 = rSges(7,3) * t191 - t172 - t131 * t150 + t139 + (t134 * t205 + t204 * t135) * t148;
t200 = (-rSges(7,1) * t147 - rSges(7,2) * t149 - t205) * t135 + (rSges(7,3) + t204) * t134;
t192 = t134 * t150;
t184 = pkin(3) * t190 + qJ(4) * t192;
t193 = qJ(4) * t134;
t199 = t141 * (pkin(3) * t135 + t193) + t150 * t184;
t198 = rSges(3,3) + qJ(2);
t197 = Icges(4,4) * t134;
t196 = Icges(4,4) * t135;
t195 = Icges(5,6) * t134;
t194 = Icges(5,6) * t135;
t114 = pkin(3) * t134 - qJ(4) * t135;
t185 = rSges(5,2) * t134 + rSges(5,3) * t135 - t114;
t182 = t141 + t142;
t60 = t102 * rSges(6,1) + t101 * rSges(6,2) + rSges(6,3) * t190;
t180 = Icges(4,5) * t217 - Icges(5,4) * t134 / 0.2e1 + (Icges(4,6) / 0.2e1 - Icges(5,5) / 0.2e1) * t135;
t179 = -pkin(8) * t134 - t114;
t144 = cos(pkin(9));
t130 = pkin(2) * t144 + pkin(1);
t146 = -pkin(7) - qJ(2);
t178 = t150 * t130 - t148 * t146;
t177 = t148 * (pkin(8) * t191 - t139) + t150 * t183 + t199;
t176 = t211 + t210 + m(7) / 0.2e1;
t80 = rSges(6,3) * t134 + (-rSges(6,1) * t147 - rSges(6,2) * t149) * t135;
t175 = t179 - t80;
t174 = rSges(4,1) * t135 - rSges(4,2) * t134;
t173 = -t104 * rSges(6,1) - t103 * rSges(6,2);
t20 = -t201 * t134 + t200 * t191;
t21 = t202 * t134 - t200 * t190;
t167 = t148 * t21 + t150 * t20;
t165 = t179 - t200;
t40 = t165 * t148;
t41 = t165 * t150;
t166 = t148 * t40 + t150 * t41;
t164 = Icges(4,1) * t135 - t197;
t163 = -Icges(4,2) * t134 + t196;
t160 = -Icges(5,2) * t135 + t195;
t159 = Icges(5,3) * t134 - t194;
t158 = rSges(4,1) * t190 - rSges(4,2) * t192 + t148 * rSges(4,3);
t157 = t148 * rSges(5,1) - rSges(5,2) * t190 + rSges(5,3) * t192;
t143 = sin(pkin(9));
t156 = rSges(3,1) * t144 - rSges(3,2) * t143 + pkin(1);
t154 = t178 + t184;
t153 = t22 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1 + t24 / 0.2e1;
t152 = t23 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1 + t25 / 0.2e1;
t30 = (t131 - t146) * t150 + (-t130 + (-qJ(4) - t205) * t134 + (-pkin(3) - t228) * t135) * t148 + t172;
t31 = t154 + t212;
t151 = m(7) * (t148 * t31 + t150 * t30);
t133 = t135 ^ 2;
t132 = t134 ^ 2;
t119 = rSges(2,1) * t150 - t148 * rSges(2,2);
t118 = -t148 * rSges(2,1) - rSges(2,2) * t150;
t72 = t198 * t148 + t156 * t150;
t71 = -t156 * t148 + t198 * t150;
t68 = t185 * t150;
t67 = t185 * t148;
t66 = t158 + t178;
t65 = (rSges(4,3) - t146) * t150 + (-t130 - t174) * t148;
t62 = rSges(6,3) * t191 - t173;
t46 = t150 * t158 + (-t150 * rSges(4,3) + t174 * t148) * t148;
t45 = t175 * t150;
t44 = t175 * t148;
t43 = t154 + t157;
t42 = (rSges(5,1) - t146) * t150 + (-t130 + (rSges(5,2) - pkin(3)) * t135 + (-rSges(5,3) - qJ(4)) * t134) * t148;
t39 = t134 * t60 - t80 * t190;
t38 = -t134 * t62 + t80 * t191;
t37 = t154 + t60 + t183;
t36 = -t146 * t150 + t139 + (-t193 - t130 + (-rSges(6,3) - pkin(3) - pkin(8)) * t135) * t148 + t173;
t35 = t150 * t157 + (-rSges(5,1) * t150 + (-rSges(5,2) * t135 + rSges(5,3) * t134) * t148) * t148 + t199;
t32 = (-t148 * t60 + t150 * t62) * t135;
t19 = t148 * t62 + t150 * t60 + t177;
t10 = (-t202 * t148 + t201 * t150) * t135;
t9 = t201 * t148 + t202 * t150 + t177;
t8 = t17 * t148 - t150 * t18;
t7 = t15 * t148 - t150 * t16;
t6 = t13 * t148 - t14 * t150;
t5 = t11 * t148 - t12 * t150;
t1 = [Icges(3,2) * t144 ^ 2 + Icges(2,3) + (Icges(3,1) * t143 + 0.2e1 * Icges(3,4) * t144) * t143 + m(7) * (t30 ^ 2 + t31 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t65 ^ 2 + t66 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + m(3) * (t71 ^ 2 + t72 ^ 2) + m(2) * (t118 ^ 2 + t119 ^ 2) + (t195 + t197 + (Icges(5,3) + Icges(4,2)) * t135 + t220) * t135 + (t194 + t196 + (Icges(4,1) + Icges(5,2)) * t134) * t134 + t221; m(7) * (t148 * t30 - t150 * t31) + m(6) * (t148 * t36 - t150 * t37) + m(4) * (t148 * t65 - t150 * t66) + m(5) * (t148 * t42 - t150 * t43) + m(3) * (t148 * t71 - t150 * t72); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t176) * t182; m(7) * (t30 * t41 + t31 * t40) + m(6) * (t36 * t45 + t37 * t44) + m(5) * (t42 * t68 + t43 * t67) + (-t65 * t206 + t180 * t150 + (Icges(5,5) * t223 + Icges(4,6) * t222 + t159 * t224 + t163 * t225) * t135 + (Icges(5,4) * t223 + Icges(4,5) * t222 + t160 * t224 + t164 * t225) * t134 - t152) * t150 + (-t66 * t206 + t180 * t148 + (Icges(5,5) * t225 + Icges(4,6) * t224 + t159 * t223 + t163 * t222) * t135 + (Icges(5,4) * t225 + Icges(4,5) * t224 + t160 * t223 + t164 * t222) * t134 + t153) * t148; m(5) * (t68 * t148 - t150 * t67) + m(6) * (t45 * t148 - t150 * t44) + m(7) * (t41 * t148 - t150 * t40); m(7) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(6) * (t19 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(5) * (t35 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t182 * t116 ^ 2 + t46 ^ 2) + (t214 * t142 - t7 - t8) * t150 + (t5 + t6 + t213 * t141 + (t214 * t148 + t213 * t150) * t150) * t148; 0.2e1 * (t151 / 0.2e1 + (t148 * t37 + t150 * t36) * t210 + (t148 * t43 + t150 * t42) * t211) * t134; 0; m(7) * (t166 * t134 - t135 * t9) + m(6) * (-t135 * t19 + (t148 * t44 + t150 * t45) * t134) + m(5) * (-t135 * t35 + (t148 * t67 + t150 * t68) * t134); 0.2e1 * t176 * (t182 * t132 + t133); m(7) * (t20 * t30 + t21 * t31) + m(6) * (t36 * t38 + t37 * t39) + (t152 * t148 + t153 * t150) * t135 + t203; m(6) * (t38 * t148 - t150 * t39) + m(7) * (t20 * t148 - t150 * t21); m(7) * (t10 * t9 + t20 * t41 + t21 * t40) + m(6) * (t19 * t32 + t38 * t45 + t39 * t44) + ((t6 / 0.2e1 + t5 / 0.2e1) * t150 + (t8 / 0.2e1 + t7 / 0.2e1) * t148) * t135 + (t216 * t148 - t215 * t150) * t217 + t219 * t224 + t218 * t223; m(6) * (-t32 * t135 + (t148 * t39 + t150 * t38) * t134) + m(7) * (-t10 * t135 + t167 * t134); t203 * t134 + m(7) * (t10 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(6) * (t32 ^ 2 + t38 ^ 2 + t39 ^ 2) + (t219 * t150 + t218 * t148 + (t215 * t148 + t216 * t150) * t134) * t135; t135 * t151; 0; m(7) * (t134 * t9 + t166 * t135); m(7) * (-0.1e1 + t182) * t135 * t134; m(7) * (t134 * t10 + t167 * t135); m(7) * (t182 * t133 + t132);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
