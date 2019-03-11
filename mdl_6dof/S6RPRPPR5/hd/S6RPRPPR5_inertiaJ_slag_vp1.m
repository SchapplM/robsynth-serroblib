% Calculate joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:06
% EndTime: 2019-03-09 02:50:10
% DurationCPUTime: 2.26s
% Computational Cost: add. (3868->335), mult. (4023->503), div. (0->0), fcn. (4080->10), ass. (0->156)
t210 = Icges(5,1) + Icges(4,3);
t128 = pkin(9) + qJ(3);
t120 = sin(t128);
t122 = cos(t128);
t209 = (-Icges(5,4) + Icges(4,5)) * t122 + (Icges(5,5) - Icges(4,6)) * t120;
t137 = sin(qJ(1));
t208 = -t137 / 0.2e1;
t196 = t137 / 0.2e1;
t138 = cos(qJ(1));
t207 = -t138 / 0.2e1;
t206 = t138 / 0.2e1;
t205 = t120 / 0.2e1;
t204 = -t209 * t137 + t210 * t138;
t203 = t210 * t137 + t209 * t138;
t133 = cos(pkin(10));
t115 = pkin(5) * t133 + pkin(4);
t135 = -pkin(8) - qJ(5);
t131 = sin(pkin(10));
t182 = t131 * t138;
t168 = t120 * t182;
t183 = t122 * t138;
t127 = pkin(10) + qJ(6);
t119 = sin(t127);
t180 = t137 * t119;
t121 = cos(t127);
t185 = t121 * t138;
t80 = t120 * t185 - t180;
t179 = t137 * t121;
t186 = t120 * t138;
t81 = t119 * t186 + t179;
t34 = t81 * rSges(7,1) + t80 * rSges(7,2) + rSges(7,3) * t183;
t202 = pkin(5) * t168 + t137 * t115 - t135 * t183 + t34;
t201 = -qJ(5) - t135;
t129 = t137 ^ 2;
t130 = t138 ^ 2;
t200 = 0.2e1 * t122;
t199 = m(5) / 0.2e1;
t198 = m(6) / 0.2e1;
t197 = m(7) / 0.2e1;
t195 = pkin(5) * t131;
t176 = pkin(3) * t183 + qJ(4) * t186;
t187 = qJ(4) * t120;
t194 = t129 * (pkin(3) * t122 + t187) + t138 * t176;
t193 = rSges(3,3) + qJ(2);
t99 = pkin(3) * t120 - qJ(4) * t122;
t192 = rSges(5,2) * t120 + rSges(5,3) * t122 - t99;
t191 = Icges(4,4) * t120;
t190 = Icges(4,4) * t122;
t189 = Icges(5,6) * t120;
t188 = Icges(5,6) * t122;
t184 = t122 * t137;
t181 = t133 * t138;
t178 = t137 * t131;
t177 = t137 * t133;
t175 = t137 * pkin(4) + qJ(5) * t183;
t174 = t129 + t130;
t173 = t198 + t197;
t86 = t120 * t181 - t178;
t87 = t168 + t177;
t172 = t87 * rSges(6,1) + t86 * rSges(6,2) + rSges(6,3) * t183;
t26 = Icges(7,5) * t81 + Icges(7,6) * t80 + Icges(7,3) * t183;
t28 = Icges(7,4) * t81 + Icges(7,2) * t80 + Icges(7,6) * t183;
t30 = Icges(7,1) * t81 + Icges(7,4) * t80 + Icges(7,5) * t183;
t10 = t120 * t26 + (-t119 * t30 - t121 * t28) * t122;
t51 = Icges(7,3) * t120 + (-Icges(7,5) * t119 - Icges(7,6) * t121) * t122;
t52 = Icges(7,6) * t120 + (-Icges(7,4) * t119 - Icges(7,2) * t121) * t122;
t53 = Icges(7,5) * t120 + (-Icges(7,1) * t119 - Icges(7,4) * t121) * t122;
t13 = t51 * t183 + t80 * t52 + t81 * t53;
t171 = t10 / 0.2e1 + t13 / 0.2e1;
t82 = t119 * t138 + t120 * t179;
t83 = t120 * t180 - t185;
t27 = Icges(7,5) * t83 + Icges(7,6) * t82 + Icges(7,3) * t184;
t29 = Icges(7,4) * t83 + Icges(7,2) * t82 + Icges(7,6) * t184;
t31 = Icges(7,1) * t83 + Icges(7,4) * t82 + Icges(7,5) * t184;
t11 = t120 * t27 + (-t119 * t31 - t121 * t29) * t122;
t14 = t51 * t184 + t52 * t82 + t53 * t83;
t170 = t11 / 0.2e1 + t14 / 0.2e1;
t169 = -Icges(5,4) * t120 / 0.2e1 + Icges(4,5) * t205 + (-Icges(5,5) / 0.2e1 + Icges(4,6) / 0.2e1) * t122;
t126 = t138 * pkin(4);
t167 = t137 * (qJ(5) * t184 - t126) + t138 * t175 + t194;
t166 = -qJ(5) * t120 - t99;
t134 = cos(pkin(9));
t116 = pkin(2) * t134 + pkin(1);
t136 = -pkin(7) - qJ(2);
t165 = t138 * t116 - t137 * t136;
t164 = t199 + t173;
t163 = t166 - rSges(6,3) * t120 - (-rSges(6,1) * t131 - rSges(6,2) * t133) * t122;
t88 = t120 * t177 + t182;
t89 = t120 * t178 - t181;
t162 = -t89 * rSges(6,1) - t88 * rSges(6,2);
t161 = -t83 * rSges(7,1) - t82 * rSges(7,2);
t160 = rSges(4,1) * t122 - rSges(4,2) * t120;
t159 = -t119 * t53 - t121 * t52;
t35 = rSges(7,3) * t184 - t161;
t54 = rSges(7,3) * t120 + (-rSges(7,1) * t119 - rSges(7,2) * t121) * t122;
t19 = -t120 * t35 + t54 * t184;
t20 = t120 * t34 - t54 * t183;
t154 = t137 * t20 + t138 * t19;
t145 = -t201 * t120 + t122 * t195 + t166 - t54;
t24 = t145 * t137;
t25 = t145 * t138;
t153 = t137 * t24 + t138 * t25;
t36 = t163 * t137;
t37 = t163 * t138;
t152 = t137 * t36 + t138 * t37;
t151 = Icges(4,1) * t122 - t191;
t150 = -Icges(4,2) * t120 + t190;
t147 = -Icges(5,2) * t122 + t189;
t146 = Icges(5,3) * t120 - t188;
t144 = rSges(4,1) * t183 - rSges(4,2) * t186 + t137 * rSges(4,3);
t143 = t137 * rSges(5,1) - rSges(5,2) * t183 + rSges(5,3) * t186;
t132 = sin(pkin(9));
t142 = rSges(3,1) * t134 - rSges(3,2) * t132 + pkin(1);
t140 = t165 + t176;
t17 = (t115 - t136) * t138 + (-t116 + (-qJ(4) - t195) * t120 + (-rSges(7,3) - pkin(3) + t135) * t122) * t137 + t161;
t18 = t140 + t202;
t22 = -t136 * t138 + t126 + (-t187 - t116 + (-rSges(6,3) - pkin(3) - qJ(5)) * t122) * t137 + t162;
t23 = t140 + t172 + t175;
t139 = (t137 * t18 + t138 * t17) * t197 + (t137 * t23 + t138 * t22) * t198;
t118 = t122 ^ 2;
t117 = t120 ^ 2;
t104 = rSges(2,1) * t138 - t137 * rSges(2,2);
t103 = -t137 * rSges(2,1) - rSges(2,2) * t138;
t101 = rSges(4,1) * t120 + rSges(4,2) * t122;
t59 = Icges(6,5) * t120 + (-Icges(6,1) * t131 - Icges(6,4) * t133) * t122;
t58 = Icges(6,6) * t120 + (-Icges(6,4) * t131 - Icges(6,2) * t133) * t122;
t56 = t193 * t137 + t142 * t138;
t55 = -t142 * t137 + t193 * t138;
t50 = t192 * t138;
t49 = t192 * t137;
t48 = t120 * t51;
t46 = t144 + t165;
t45 = (rSges(4,3) - t136) * t138 + (-t116 - t160) * t137;
t44 = Icges(6,1) * t89 + Icges(6,4) * t88 + Icges(6,5) * t184;
t43 = Icges(6,1) * t87 + Icges(6,4) * t86 + Icges(6,5) * t183;
t42 = Icges(6,4) * t89 + Icges(6,2) * t88 + Icges(6,6) * t184;
t41 = Icges(6,4) * t87 + Icges(6,2) * t86 + Icges(6,6) * t183;
t40 = Icges(6,5) * t89 + Icges(6,6) * t88 + Icges(6,3) * t184;
t39 = Icges(6,5) * t87 + Icges(6,6) * t86 + Icges(6,3) * t183;
t38 = t138 * t144 + (-t138 * rSges(4,3) + t160 * t137) * t137;
t33 = t140 + t143;
t32 = (rSges(5,1) - t136) * t138 + (-t116 + (rSges(5,2) - pkin(3)) * t122 + (-rSges(5,3) - qJ(4)) * t120) * t137;
t21 = t138 * t143 + (-t138 * rSges(5,1) + (-rSges(5,2) * t122 + rSges(5,3) * t120) * t137) * t137 + t194;
t16 = (t159 * t122 + t48) * t120;
t15 = (-t137 * t34 + t138 * t35) * t122;
t12 = t137 * (rSges(6,3) * t184 - t162) + t138 * t172 + t167;
t9 = t27 * t184 + t29 * t82 + t31 * t83;
t8 = t26 * t184 + t28 * t82 + t30 * t83;
t7 = t27 * t183 + t80 * t29 + t81 * t31;
t6 = t26 * t183 + t80 * t28 + t81 * t30;
t5 = (-t175 + t202) * t138 + (-t138 * t115 + t126 + t35 + (t120 * t195 + t201 * t122) * t137) * t137 + t167;
t4 = t8 * t137 - t138 * t9;
t3 = t6 * t137 - t138 * t7;
t2 = t14 * t120 + (t137 * t9 + t138 * t8) * t122;
t1 = t13 * t120 + (t137 * t7 + t138 * t6) * t122;
t47 = [Icges(3,2) * t134 ^ 2 + Icges(2,3) + t48 + (Icges(3,1) * t132 + 0.2e1 * Icges(3,4) * t134) * t132 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t45 ^ 2 + t46 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t103 ^ 2 + t104 ^ 2) + (-t131 * t59 - t133 * t58 + t159 + t189 + t191 + (Icges(4,2) + Icges(5,3)) * t122) * t122 + (t188 + (-Icges(6,5) * t131 - Icges(6,6) * t133) * t122 + t190 + (Icges(5,2) + Icges(6,3) + Icges(4,1)) * t120) * t120; m(7) * (t137 * t17 - t138 * t18) + m(5) * (t137 * t32 - t138 * t33) + m(6) * (t137 * t22 - t138 * t23) + m(4) * (t137 * t45 - t138 * t46) + m(3) * (t137 * t55 - t138 * t56); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t164) * t174; (-t58 * t88 / 0.2e1 - t59 * t89 / 0.2e1 + t169 * t138 - t170) * t138 + (t86 * t58 / 0.2e1 + t87 * t59 / 0.2e1 + t169 * t137 + t171) * t137 + m(7) * (t17 * t25 + t18 * t24) + m(5) * (t32 * t50 + t33 * t49) + m(6) * (t22 * t37 + t23 * t36) + m(4) * (-t137 * t46 - t138 * t45) * t101 + ((Icges(5,4) * t207 + t147 * t196 - t40 / 0.2e1 + Icges(4,5) * t206 + t151 * t208) * t138 + (t39 / 0.2e1 + Icges(4,5) * t196 + t151 * t206 + Icges(5,4) * t208 + t147 * t207) * t137) * t120 + ((Icges(5,5) * t207 + t146 * t196 + t131 * t44 / 0.2e1 + t133 * t42 / 0.2e1 + Icges(4,6) * t206 + t150 * t208) * t138 + (-t131 * t43 / 0.2e1 - t133 * t41 / 0.2e1 + Icges(4,6) * t196 + t150 * t206 + Icges(5,5) * t208 + t146 * t207) * t137) * t122; m(5) * (t50 * t137 - t138 * t49) + m(6) * (t37 * t137 - t138 * t36) + m(7) * (t25 * t137 - t138 * t24); m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(5) * (t21 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t12 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t174 * t101 ^ 2 + t38 ^ 2) + (-t4 + (t40 * t184 + t88 * t42 + t89 * t44) * t138 + t204 * t130) * t138 + (t3 + (t39 * t183 + t86 * t41 + t87 * t43) * t137 + t203 * t129 + (t204 * t137 + t203 * t138 - t40 * t183 - t39 * t184 - t41 * t88 - t86 * t42 - t43 * t89 - t87 * t44) * t138) * t137; 0.2e1 * ((t137 * t33 + t138 * t32) * t199 + t139) * t120; 0; m(7) * (t153 * t120 - t122 * t5) + m(5) * (-t122 * t21 + (t137 * t49 + t138 * t50) * t120) + m(6) * (-t122 * t12 + t152 * t120); 0.2e1 * t164 * (t174 * t117 + t118); t139 * t200; 0; m(7) * (t120 * t5 + t153 * t122) + m(6) * (t120 * t12 + t152 * t122); t173 * (-0.1e1 + t174) * t120 * t200; 0.2e1 * t173 * (t174 * t118 + t117); t16 + m(7) * (t17 * t19 + t18 * t20) + (t170 * t137 + t171 * t138) * t122; m(7) * (t19 * t137 - t138 * t20); m(7) * (t15 * t5 + t19 * t25 + t20 * t24) + t2 * t207 + t1 * t196 + (t10 * t137 - t11 * t138) * t205 + (t4 * t196 + t3 * t206) * t122; m(7) * (t154 * t120 - t15 * t122); m(7) * (t15 * t120 + t154 * t122); t120 * t16 + m(7) * (t15 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t138 * t1 + t137 * t2 + t120 * (t10 * t138 + t11 * t137)) * t122;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t47(1) t47(2) t47(4) t47(7) t47(11) t47(16); t47(2) t47(3) t47(5) t47(8) t47(12) t47(17); t47(4) t47(5) t47(6) t47(9) t47(13) t47(18); t47(7) t47(8) t47(9) t47(10) t47(14) t47(19); t47(11) t47(12) t47(13) t47(14) t47(15) t47(20); t47(16) t47(17) t47(18) t47(19) t47(20) t47(21);];
Mq  = res;
