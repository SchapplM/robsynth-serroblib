% Calculate joint inertia matrix for
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP11_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP11_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:23
% EndTime: 2019-12-31 20:12:29
% DurationCPUTime: 2.14s
% Computational Cost: add. (2089->315), mult. (5065->464), div. (0->0), fcn. (5412->6), ass. (0->151)
t221 = Icges(4,1) + Icges(3,3);
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t220 = (-Icges(4,4) + Icges(3,5)) * t143 + (Icges(4,5) - Icges(3,6)) * t140;
t139 = sin(qJ(4));
t144 = cos(qJ(1));
t141 = sin(qJ(1));
t142 = cos(qJ(4));
t179 = t141 * t142;
t102 = t139 * t144 + t140 * t179;
t176 = t142 * t144;
t181 = t139 * t141;
t104 = t140 * t181 - t176;
t213 = rSges(6,3) + qJ(5);
t214 = rSges(6,1) + pkin(4);
t219 = t213 * t102 - t214 * t104;
t218 = -t141 / 0.2e1;
t217 = t141 / 0.2e1;
t216 = -t144 / 0.2e1;
t215 = t144 / 0.2e1;
t177 = t142 * t143;
t72 = Icges(6,6) * t140 + (-Icges(6,5) * t139 + Icges(6,3) * t142) * t143;
t73 = Icges(5,3) * t140 + (-Icges(5,5) * t139 - Icges(5,6) * t142) * t143;
t76 = Icges(6,2) * t140 + (-Icges(6,4) * t139 + Icges(6,6) * t142) * t143;
t212 = t72 * t177 + (t73 + t76) * t140;
t77 = Icges(5,6) * t140 + (-Icges(5,4) * t139 - Icges(5,2) * t142) * t143;
t80 = Icges(6,4) * t140 + (-Icges(6,1) * t139 + Icges(6,5) * t142) * t143;
t81 = Icges(5,5) * t140 + (-Icges(5,1) * t139 - Icges(5,4) * t142) * t143;
t211 = -t142 * t77 + (-t80 - t81) * t139;
t100 = -t140 * t176 + t181;
t180 = t140 * t144;
t101 = t139 * t180 + t179;
t175 = t143 * t144;
t45 = Icges(6,5) * t101 + Icges(6,6) * t175 + Icges(6,3) * t100;
t49 = Icges(6,4) * t101 + Icges(6,2) * t175 + Icges(6,6) * t100;
t53 = Icges(6,1) * t101 + Icges(6,4) * t175 + Icges(6,5) * t100;
t11 = t100 * t45 + t101 * t53 + t49 * t175;
t178 = t141 * t143;
t46 = Icges(6,5) * t104 + Icges(6,6) * t178 - Icges(6,3) * t102;
t50 = Icges(6,4) * t104 + Icges(6,2) * t178 - Icges(6,6) * t102;
t54 = Icges(6,1) * t104 + Icges(6,4) * t178 - Icges(6,5) * t102;
t12 = t100 * t46 + t101 * t54 + t50 * t175;
t47 = Icges(5,5) * t101 - Icges(5,6) * t100 + Icges(5,3) * t175;
t51 = Icges(5,4) * t101 - Icges(5,2) * t100 + Icges(5,6) * t175;
t55 = Icges(5,1) * t101 - Icges(5,4) * t100 + Icges(5,5) * t175;
t13 = -t100 * t51 + t101 * t55 + t47 * t175;
t48 = Icges(5,5) * t104 + Icges(5,6) * t102 + Icges(5,3) * t178;
t52 = Icges(5,4) * t104 + Icges(5,2) * t102 + Icges(5,6) * t178;
t56 = Icges(5,1) * t104 + Icges(5,4) * t102 + Icges(5,5) * t178;
t14 = -t100 * t52 + t101 * t56 + t48 * t175;
t28 = t100 * t72 + t101 * t80 + t76 * t175;
t29 = -t100 * t77 + t101 * t81 + t73 * t175;
t210 = ((t11 + t13) * t144 + (t12 + t14) * t141) * t143 + (t28 + t29) * t140;
t15 = -t102 * t45 + t104 * t53 + t49 * t178;
t16 = -t102 * t46 + t104 * t54 + t50 * t178;
t17 = t102 * t51 + t104 * t55 + t47 * t178;
t18 = t102 * t52 + t104 * t56 + t48 * t178;
t30 = -t102 * t72 + t104 * t80 + t76 * t178;
t31 = t102 * t77 + t104 * t81 + t73 * t178;
t209 = ((t15 + t17) * t144 + (t16 + t18) * t141) * t143 + (t30 + t31) * t140;
t208 = t140 / 0.2e1;
t20 = t140 * t49 + (-t139 * t53 + t142 * t45) * t143;
t22 = t140 * t47 + (-t139 * t55 - t142 * t51) * t143;
t207 = t20 + t22;
t21 = t140 * t50 + (-t139 * t54 + t142 * t46) * t143;
t23 = t140 * t48 + (-t139 * t56 - t142 * t52) * t143;
t206 = t21 + t23;
t205 = t221 * t141 + t220 * t144;
t204 = -t220 * t141 + t221 * t144;
t136 = t141 ^ 2;
t138 = t144 ^ 2;
t203 = m(4) / 0.2e1;
t202 = m(5) / 0.2e1;
t201 = m(6) / 0.2e1;
t200 = -pkin(2) - pkin(7);
t117 = rSges(3,1) * t140 + rSges(3,2) * t143;
t196 = m(3) * t117;
t195 = (t211 * t143 + t212) * t140;
t194 = rSges(6,2) * t175 + t213 * t100 + t214 * t101;
t193 = rSges(6,2) * t178 - t219;
t173 = pkin(2) * t175 + qJ(3) * t180;
t182 = qJ(3) * t140;
t191 = t136 * (pkin(2) * t143 + t182) + t144 * t173;
t189 = t144 * rSges(4,1);
t188 = t144 * rSges(3,3);
t187 = rSges(6,2) * t140 + (-t214 * t139 + t213 * t142) * t143;
t186 = Icges(3,4) * t140;
t185 = Icges(3,4) * t143;
t184 = Icges(4,6) * t140;
t183 = Icges(4,6) * t143;
t115 = pkin(2) * t140 - qJ(3) * t143;
t174 = rSges(4,2) * t140 + rSges(4,3) * t143 - t115;
t172 = t144 * pkin(1) + t141 * pkin(6);
t171 = t141 * pkin(3) + pkin(7) * t175;
t132 = t144 * pkin(6);
t133 = t144 * pkin(3);
t170 = t132 + t133;
t169 = t136 + t138;
t58 = t101 * rSges(5,1) - t100 * rSges(5,2) + rSges(5,3) * t175;
t168 = -Icges(4,4) * t140 / 0.2e1 + Icges(3,5) * t208 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t143;
t167 = -pkin(1) - t182;
t166 = -pkin(7) * t140 - t115;
t165 = t141 * (pkin(7) * t178 - t133) + t144 * t171 + t191;
t164 = t172 + t173;
t91 = rSges(5,3) * t140 + (-rSges(5,1) * t139 - rSges(5,2) * t142) * t143;
t163 = t166 - t91;
t162 = rSges(3,1) * t143 - rSges(3,2) * t140;
t161 = -rSges(5,1) * t104 - rSges(5,2) * t102;
t156 = Icges(3,1) * t143 - t186;
t155 = -Icges(3,2) * t140 + t185;
t152 = -Icges(4,2) * t143 + t184;
t151 = Icges(4,3) * t140 - t183;
t150 = t166 - t187;
t149 = rSges(3,1) * t175 - rSges(3,2) * t180 + t141 * rSges(3,3);
t148 = t141 * rSges(4,1) - rSges(4,2) * t175 + rSges(4,3) * t180;
t147 = t23 / 0.2e1 + t21 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1;
t146 = t29 / 0.2e1 + t28 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t145 = t164 + t171;
t137 = t143 ^ 2;
t119 = rSges(2,1) * t144 - t141 * rSges(2,2);
t118 = -t141 * rSges(2,1) - rSges(2,2) * t144;
t68 = t174 * t144;
t67 = t174 * t141;
t66 = t149 + t172;
t65 = t188 + t132 + (-pkin(1) - t162) * t141;
t62 = t148 + t164;
t61 = t189 + t132 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t143 + (-rSges(4,3) - qJ(3)) * t140) * t141;
t60 = rSges(5,3) * t178 - t161;
t44 = t163 * t144;
t43 = t163 * t141;
t42 = t144 * t149 + (t162 * t141 - t188) * t141;
t41 = t150 * t144;
t40 = t150 * t141;
t39 = t140 * t58 - t91 * t175;
t38 = -t140 * t60 + t91 * t178;
t37 = t145 + t58;
t36 = ((-rSges(5,3) + t200) * t143 + t167) * t141 + t161 + t170;
t33 = t144 * t148 + (-t189 + (-rSges(4,2) * t143 + rSges(4,3) * t140) * t141) * t141 + t191;
t32 = (-t141 * t58 + t144 * t60) * t143;
t27 = t145 + t194;
t26 = ((-rSges(6,2) + t200) * t143 + t167) * t141 + t170 + t219;
t25 = t194 * t140 - t187 * t175;
t24 = -t193 * t140 + t187 * t178;
t19 = t141 * t60 + t144 * t58 + t165;
t10 = (-t194 * t141 + t193 * t144) * t143;
t9 = t193 * t141 + t194 * t144 + t165;
t8 = t17 * t141 - t144 * t18;
t7 = t15 * t141 - t144 * t16;
t6 = t13 * t141 - t14 * t144;
t5 = t11 * t141 - t12 * t144;
t1 = [Icges(2,3) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(6) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t61 ^ 2 + t62 ^ 2) + m(3) * (t65 ^ 2 + t66 ^ 2) + m(2) * (t118 ^ 2 + t119 ^ 2) + (t184 + t186 + (Icges(4,3) + Icges(3,2)) * t143 + t211) * t143 + (t183 + t185 + (Icges(3,1) + Icges(4,2)) * t140) * t140 + t212; m(5) * (t36 * t44 + t37 * t43) + m(6) * (t26 * t41 + t27 * t40) + m(4) * (t61 * t68 + t62 * t67) + (-t65 * t196 + t168 * t144 + (Icges(4,5) * t216 + Icges(3,6) * t215 + t151 * t217 + t155 * t218) * t143 + (Icges(4,4) * t216 + Icges(3,5) * t215 + t152 * t217 + t156 * t218) * t140 - t147) * t144 + (-t66 * t196 + t168 * t141 + (Icges(4,5) * t218 + Icges(3,6) * t217 + t151 * t216 + t155 * t215) * t143 + (Icges(4,4) * t218 + Icges(3,5) * t217 + t152 * t216 + t156 * t215) * t140 + t146) * t141; m(6) * (t40 ^ 2 + t41 ^ 2 + t9 ^ 2) + m(5) * (t19 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(4) * (t33 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t169 * t117 ^ 2 + t42 ^ 2) + (t204 * t138 - t7 - t8) * t144 + (t5 + t6 + t205 * t136 + (t204 * t141 + t205 * t144) * t144) * t141; 0.2e1 * ((t141 * t37 + t144 * t36) * t202 + (t141 * t27 + t144 * t26) * t201 + (t141 * t62 + t144 * t61) * t203) * t140; m(6) * (-t143 * t9 + (t141 * t40 + t144 * t41) * t140) + m(5) * (-t143 * t19 + (t141 * t43 + t144 * t44) * t140) + m(4) * (-t143 * t33 + (t141 * t67 + t144 * t68) * t140); 0.2e1 * (t203 + t202 + t201) * (t169 * t140 ^ 2 + t137); m(5) * (t36 * t38 + t37 * t39) + m(6) * (t24 * t26 + t25 * t27) + (t147 * t141 + t146 * t144) * t143 + t195; m(6) * (t10 * t9 + t24 * t41 + t25 * t40) + m(5) * (t19 * t32 + t38 * t44 + t39 * t43) + ((t5 / 0.2e1 + t6 / 0.2e1) * t144 + (t8 / 0.2e1 + t7 / 0.2e1) * t141) * t143 + (t207 * t141 - t206 * t144) * t208 + t210 * t217 + t209 * t216; m(5) * (-t32 * t143 + (t141 * t39 + t144 * t38) * t140) + m(6) * (-t10 * t143 + (t141 * t25 + t144 * t24) * t140); t195 * t140 + m(6) * (t10 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t32 ^ 2 + t38 ^ 2 + t39 ^ 2) + (t210 * t144 + t209 * t141 + (t206 * t141 + t207 * t144) * t140) * t143; m(6) * (t100 * t26 - t102 * t27); m(6) * (t100 * t41 - t102 * t40 + t9 * t177); m(6) * (-t137 * t142 + (t100 * t144 - t102 * t141) * t140); m(6) * (t10 * t177 + t100 * t24 - t102 * t25); m(6) * (t137 * t142 ^ 2 + t100 ^ 2 + t102 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
