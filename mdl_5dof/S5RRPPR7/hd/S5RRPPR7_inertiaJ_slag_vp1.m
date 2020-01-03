% Calculate joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:35:02
% DurationCPUTime: 1.73s
% Computational Cost: add. (2322->268), mult. (3057->401), div. (0->0), fcn. (3101->8), ass. (0->131)
t191 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t111 = qJ(2) + pkin(8);
t104 = sin(t111);
t105 = cos(t111);
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t190 = Icges(3,5) * t119 - Icges(3,6) * t116 + (-Icges(5,4) + Icges(4,5)) * t105 + (Icges(5,5) - Icges(4,6)) * t104;
t117 = sin(qJ(1));
t189 = -t117 / 0.2e1;
t176 = t117 / 0.2e1;
t120 = cos(qJ(1));
t188 = -t120 / 0.2e1;
t187 = t120 / 0.2e1;
t186 = t117 * pkin(6);
t185 = t104 / 0.2e1;
t184 = t116 / 0.2e1;
t183 = t119 / 0.2e1;
t159 = t105 * t120;
t118 = cos(qJ(5));
t155 = t118 * t120;
t115 = sin(qJ(5));
t157 = t117 * t115;
t73 = t104 * t155 - t157;
t156 = t117 * t118;
t158 = t115 * t120;
t74 = t104 * t158 + t156;
t32 = t74 * rSges(6,1) + t73 * rSges(6,2) + rSges(6,3) * t159;
t182 = t117 * pkin(4) + pkin(7) * t159 + t32;
t181 = t191 * t117 + t190 * t120;
t180 = -t190 * t117 + t191 * t120;
t112 = t117 ^ 2;
t113 = t120 ^ 2;
t179 = m(5) / 0.2e1;
t178 = m(6) / 0.2e1;
t90 = rSges(3,1) * t116 + rSges(3,2) * t119;
t177 = m(3) * t90;
t175 = pkin(2) * t116;
t102 = pkin(2) * t119 + pkin(1);
t110 = t120 * pkin(6);
t99 = t120 * t102;
t174 = t117 * (t110 + (-pkin(1) + t102) * t117) + t120 * (-pkin(1) * t120 - t186 + t99);
t161 = t104 * t120;
t173 = pkin(3) * t159 + qJ(4) * t161;
t172 = rSges(3,1) * t119;
t171 = rSges(3,2) * t116;
t170 = t120 * rSges(3,3);
t168 = Icges(3,4) * t116;
t167 = Icges(3,4) * t119;
t166 = Icges(4,4) * t104;
t165 = Icges(4,4) * t105;
t164 = Icges(5,6) * t104;
t163 = Icges(5,6) * t105;
t162 = qJ(4) * t104;
t160 = t105 * t117;
t154 = t117 * rSges(3,3) + t120 * t172;
t153 = t112 + t113;
t152 = t179 + t178;
t26 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t159;
t28 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t159;
t30 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t159;
t10 = t104 * t26 + (-t115 * t30 - t118 * t28) * t105;
t44 = Icges(6,3) * t104 + (-Icges(6,5) * t115 - Icges(6,6) * t118) * t105;
t45 = Icges(6,6) * t104 + (-Icges(6,4) * t115 - Icges(6,2) * t118) * t105;
t46 = Icges(6,5) * t104 + (-Icges(6,1) * t115 - Icges(6,4) * t118) * t105;
t13 = t159 * t44 + t73 * t45 + t74 * t46;
t151 = t10 / 0.2e1 + t13 / 0.2e1;
t75 = t104 * t156 + t158;
t76 = t104 * t157 - t155;
t27 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t160;
t29 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t160;
t31 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t160;
t11 = t104 * t27 + (-t115 * t31 - t118 * t29) * t105;
t14 = t160 * t44 + t45 * t75 + t46 * t76;
t150 = t11 / 0.2e1 + t14 / 0.2e1;
t149 = -pkin(3) * t104 + qJ(4) * t105 - t175;
t148 = -rSges(4,1) * t104 - rSges(4,2) * t105 - t175;
t147 = t112 * (pkin(3) * t105 + t162) + t120 * t173 + t174;
t114 = -qJ(3) - pkin(6);
t146 = -t117 * t114 + t99;
t145 = rSges(5,2) * t104 + rSges(5,3) * t105 + t149;
t144 = Icges(3,5) * t184 + Icges(3,6) * t183 + Icges(4,5) * t185 - Icges(5,4) * t104 / 0.2e1 + (Icges(4,6) / 0.2e1 - Icges(5,5) / 0.2e1) * t105;
t143 = -t76 * rSges(6,1) - t75 * rSges(6,2);
t142 = -t171 + t172;
t141 = rSges(4,1) * t105 - rSges(4,2) * t104;
t136 = -t115 * t46 - t118 * t45;
t133 = Icges(3,1) * t119 - t168;
t132 = Icges(4,1) * t105 - t166;
t131 = -Icges(3,2) * t116 + t167;
t130 = -Icges(4,2) * t104 + t165;
t126 = -Icges(5,2) * t105 + t164;
t125 = Icges(5,3) * t104 - t163;
t124 = t146 + t173;
t123 = rSges(4,1) * t159 - rSges(4,2) * t161 + t117 * rSges(4,3);
t122 = t117 * rSges(5,1) - rSges(5,2) * t159 + rSges(5,3) * t161;
t47 = rSges(6,3) * t104 + (-rSges(6,1) * t115 - rSges(6,2) * t118) * t105;
t121 = -pkin(7) * t104 + t149 - t47;
t92 = rSges(2,1) * t120 - t117 * rSges(2,2);
t91 = -t117 * rSges(2,1) - rSges(2,2) * t120;
t49 = t148 * t120;
t48 = t148 * t117;
t43 = t186 + (pkin(1) - t171) * t120 + t154;
t42 = t170 + t110 + (-pkin(1) - t142) * t117;
t39 = t104 * t44;
t38 = t123 + t146;
t37 = (rSges(4,3) - t114) * t120 + (-t102 - t141) * t117;
t36 = t145 * t120;
t35 = t145 * t117;
t34 = t120 * (-t120 * t171 + t154) + (t117 * t142 - t170) * t117;
t33 = rSges(6,3) * t160 - t143;
t25 = t122 + t124;
t24 = (rSges(5,1) - t114) * t120 + (-t102 + (rSges(5,2) - pkin(3)) * t105 + (-rSges(5,3) - qJ(4)) * t104) * t117;
t23 = t121 * t120;
t22 = t121 * t117;
t21 = t104 * t32 - t159 * t47;
t20 = -t104 * t33 + t160 * t47;
t19 = t124 + t182;
t18 = (pkin(4) - t114) * t120 + (-t162 - t102 + (-rSges(6,3) - pkin(3) - pkin(7)) * t105) * t117 + t143;
t17 = (t105 * t136 + t39) * t104;
t16 = t120 * t123 + (-t120 * rSges(4,3) + t117 * t141) * t117 + t174;
t15 = (-t117 * t32 + t120 * t33) * t105;
t12 = t120 * t122 + (-t120 * rSges(5,1) + (-rSges(5,2) * t105 + rSges(5,3) * t104) * t117) * t117 + t147;
t9 = t160 * t27 + t29 * t75 + t31 * t76;
t8 = t160 * t26 + t28 * t75 + t30 * t76;
t7 = t159 * t27 + t73 * t29 + t74 * t31;
t6 = t159 * t26 + t73 * t28 + t74 * t30;
t5 = t182 * t120 + (-pkin(4) * t120 + pkin(7) * t160 + t33) * t117 + t147;
t4 = t8 * t117 - t120 * t9;
t3 = t6 * t117 - t120 * t7;
t2 = t14 * t104 + (t117 * t9 + t120 * t8) * t105;
t1 = t13 * t104 + (t117 * t7 + t120 * t6) * t105;
t40 = [t116 * (Icges(3,1) * t116 + t167) + t119 * (Icges(3,2) * t119 + t168) + Icges(2,3) + t39 + m(6) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(2) * (t91 ^ 2 + t92 ^ 2) + (t136 + t164 + t166 + (Icges(4,2) + Icges(5,3)) * t105) * t105 + (t163 + t165 + (Icges(4,1) + Icges(5,2)) * t104) * t104; m(6) * (t18 * t23 + t19 * t22) + m(4) * (t37 * t49 + t38 * t48) + m(5) * (t24 * t36 + t25 * t35) + (-t42 * t177 - t116 * (-Icges(3,5) * t120 + t117 * t133) / 0.2e1 - t119 * (-Icges(3,6) * t120 + t117 * t131) / 0.2e1 + t144 * t120 + (Icges(5,5) * t188 + Icges(4,6) * t187 + t125 * t176 + t130 * t189) * t105 - t150) * t120 + ((Icges(3,6) * t117 + t120 * t131) * t183 - t43 * t177 + t144 * t117 + (Icges(5,5) * t189 + Icges(4,6) * t176 + t125 * t188 + t130 * t187) * t105 + (Icges(3,5) * t117 + t120 * t133) * t184 + t151) * t117 + ((Icges(5,4) * t188 + Icges(4,5) * t187 + t126 * t176 + t132 * t189) * t120 + (Icges(5,4) * t189 + Icges(4,5) * t176 + t126 * t188 + t132 * t187) * t117) * t104; m(6) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(5) * (t12 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(4) * (t16 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(3) * (t153 * t90 ^ 2 + t34 ^ 2) + (t180 * t113 - t4) * t120 + (t3 + t181 * t112 + (t180 * t117 + t181 * t120) * t120) * t117; m(6) * (t117 * t18 - t120 * t19) + m(4) * (t117 * t37 - t120 * t38) + m(5) * (t117 * t24 - t120 * t25); m(6) * (t117 * t23 - t120 * t22) + m(5) * (t117 * t36 - t120 * t35) + m(4) * (t117 * t49 - t120 * t48); 0.2e1 * (m(4) / 0.2e1 + t152) * t153; 0.2e1 * ((t117 * t19 + t120 * t18) * t178 + (t117 * t25 + t120 * t24) * t179) * t104; m(6) * (-t105 * t5 + (t117 * t22 + t120 * t23) * t104) + m(5) * (-t105 * t12 + (t117 * t35 + t120 * t36) * t104); 0; 0.2e1 * t152 * (t104 ^ 2 * t153 + t105 ^ 2); m(6) * (t18 * t20 + t19 * t21) + t17 + (t117 * t150 + t120 * t151) * t105; t1 * t176 + t2 * t188 + m(6) * (t15 * t5 + t20 * t23 + t21 * t22) + (t10 * t117 - t11 * t120) * t185 + (t4 * t176 + t3 * t187) * t105; m(6) * (t20 * t117 - t120 * t21); m(6) * (-t15 * t105 + (t117 * t21 + t120 * t20) * t104); m(6) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + t104 * t17 + (t120 * t1 + t117 * t2 + t104 * (t10 * t120 + t11 * t117)) * t105;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t40(1), t40(2), t40(4), t40(7), t40(11); t40(2), t40(3), t40(5), t40(8), t40(12); t40(4), t40(5), t40(6), t40(9), t40(13); t40(7), t40(8), t40(9), t40(10), t40(14); t40(11), t40(12), t40(13), t40(14), t40(15);];
Mq = res;
