% Calculate joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:20
% DurationCPUTime: 1.78s
% Computational Cost: add. (3926->289), mult. (3277->416), div. (0->0), fcn. (3291->10), ass. (0->139)
t199 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t118 = qJ(3) + pkin(10);
t113 = sin(t118);
t115 = cos(t118);
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t198 = Icges(4,5) * t125 - Icges(4,6) * t122 + (-Icges(6,4) + Icges(5,5)) * t115 + (Icges(6,5) - Icges(5,6)) * t113;
t119 = qJ(1) + pkin(9);
t114 = sin(t119);
t197 = -t114 / 0.2e1;
t184 = t114 / 0.2e1;
t116 = cos(t119);
t196 = -t116 / 0.2e1;
t195 = t116 / 0.2e1;
t194 = pkin(7) * t114;
t193 = t113 / 0.2e1;
t192 = t122 / 0.2e1;
t191 = t125 / 0.2e1;
t163 = t115 * t116;
t124 = cos(qJ(6));
t161 = t116 * t124;
t121 = sin(qJ(6));
t165 = t114 * t121;
t73 = t113 * t161 - t165;
t162 = t116 * t121;
t164 = t114 * t124;
t74 = t113 * t162 + t164;
t33 = t74 * rSges(7,1) + t73 * rSges(7,2) + rSges(7,3) * t163;
t190 = t114 * pkin(5) + pkin(8) * t163 + t33;
t189 = -t198 * t114 + t116 * t199;
t188 = t114 * t199 + t198 * t116;
t111 = t114 ^ 2;
t112 = t116 ^ 2;
t187 = m(6) / 0.2e1;
t186 = m(7) / 0.2e1;
t185 = -m(6) - m(7);
t100 = rSges(4,1) * t122 + rSges(4,2) * t125;
t183 = m(4) * t100;
t123 = sin(qJ(1));
t182 = pkin(1) * t123;
t181 = pkin(3) * t122;
t108 = t116 * pkin(7);
t109 = pkin(3) * t125 + pkin(2);
t95 = t116 * t109;
t180 = t114 * (t108 + (-pkin(2) + t109) * t114) + t116 * (-pkin(2) * t116 - t194 + t95);
t167 = t113 * t116;
t179 = pkin(4) * t163 + qJ(5) * t167;
t178 = rSges(4,1) * t125;
t177 = rSges(4,2) * t122;
t176 = t116 * rSges(4,3);
t174 = Icges(4,4) * t122;
t173 = Icges(4,4) * t125;
t172 = Icges(5,4) * t113;
t171 = Icges(5,4) * t115;
t170 = Icges(6,6) * t113;
t169 = Icges(6,6) * t115;
t168 = qJ(5) * t113;
t166 = t114 * t115;
t160 = t114 * rSges(4,3) + t116 * t178;
t159 = t111 + t112;
t158 = t187 + t186;
t27 = Icges(7,5) * t74 + Icges(7,6) * t73 + Icges(7,3) * t163;
t29 = Icges(7,4) * t74 + Icges(7,2) * t73 + Icges(7,6) * t163;
t31 = Icges(7,1) * t74 + Icges(7,4) * t73 + Icges(7,5) * t163;
t11 = t113 * t27 + (-t121 * t31 - t124 * t29) * t115;
t61 = Icges(7,3) * t113 + (-Icges(7,5) * t121 - Icges(7,6) * t124) * t115;
t64 = Icges(7,6) * t113 + (-Icges(7,4) * t121 - Icges(7,2) * t124) * t115;
t67 = Icges(7,5) * t113 + (-Icges(7,1) * t121 - Icges(7,4) * t124) * t115;
t15 = t61 * t163 + t64 * t73 + t67 * t74;
t157 = t11 / 0.2e1 + t15 / 0.2e1;
t75 = t113 * t164 + t162;
t76 = t113 * t165 - t161;
t28 = Icges(7,5) * t76 + Icges(7,6) * t75 + Icges(7,3) * t166;
t30 = Icges(7,4) * t76 + Icges(7,2) * t75 + Icges(7,6) * t166;
t32 = Icges(7,1) * t76 + Icges(7,4) * t75 + Icges(7,5) * t166;
t12 = t113 * t28 + (-t121 * t32 - t124 * t30) * t115;
t16 = t61 * t166 + t64 * t75 + t67 * t76;
t156 = t12 / 0.2e1 + t16 / 0.2e1;
t155 = -pkin(4) * t113 + qJ(5) * t115 - t181;
t154 = -rSges(5,1) * t113 - rSges(5,2) * t115 - t181;
t153 = t111 * (pkin(4) * t115 + t168) + t116 * t179 + t180;
t152 = rSges(6,2) * t113 + rSges(6,3) * t115 + t155;
t151 = Icges(4,5) * t192 + Icges(4,6) * t191 + Icges(5,5) * t193 - Icges(6,4) * t113 / 0.2e1 + (Icges(5,6) / 0.2e1 - Icges(6,5) / 0.2e1) * t115;
t126 = cos(qJ(1));
t117 = t126 * pkin(1);
t120 = -qJ(4) - pkin(7);
t150 = -t114 * t120 + t117 + t95;
t149 = -rSges(7,1) * t76 - rSges(7,2) * t75;
t148 = -t177 + t178;
t147 = rSges(5,1) * t115 - rSges(5,2) * t113;
t142 = -t121 * t67 - t124 * t64;
t139 = Icges(4,1) * t125 - t174;
t138 = Icges(5,1) * t115 - t172;
t137 = -Icges(4,2) * t122 + t173;
t136 = -Icges(5,2) * t113 + t171;
t132 = -Icges(6,2) * t115 + t170;
t131 = Icges(6,3) * t113 - t169;
t130 = rSges(5,1) * t163 - rSges(5,2) * t167 + t114 * rSges(5,3);
t129 = t114 * rSges(6,1) - rSges(6,2) * t163 + rSges(6,3) * t167;
t128 = t150 + t179;
t70 = rSges(7,3) * t113 + (-rSges(7,1) * t121 - rSges(7,2) * t124) * t115;
t127 = -pkin(8) * t113 + t155 - t70;
t102 = rSges(2,1) * t126 - t123 * rSges(2,2);
t101 = -t123 * rSges(2,1) - rSges(2,2) * t126;
t78 = rSges(3,1) * t116 - rSges(3,2) * t114 + t117;
t77 = -rSges(3,1) * t114 - rSges(3,2) * t116 - t182;
t57 = t154 * t116;
t56 = t154 * t114;
t43 = t113 * t61;
t40 = t194 + t117 + (pkin(2) - t177) * t116 + t160;
t39 = t176 - t182 + t108 + (-pkin(2) - t148) * t114;
t38 = t152 * t116;
t37 = t152 * t114;
t36 = t130 + t150;
t35 = -t182 + (rSges(5,3) - t120) * t116 + (-t109 - t147) * t114;
t34 = rSges(7,3) * t166 - t149;
t26 = t116 * (-t116 * t177 + t160) + (t114 * t148 - t176) * t114;
t25 = t127 * t116;
t24 = t127 * t114;
t23 = t128 + t129;
t22 = -t182 + (rSges(6,1) - t120) * t116 + (-t109 + (rSges(6,2) - pkin(4)) * t115 + (-rSges(6,3) - qJ(5)) * t113) * t114;
t21 = t113 * t33 - t70 * t163;
t20 = -t113 * t34 + t70 * t166;
t19 = (t115 * t142 + t43) * t113;
t18 = t128 + t190;
t17 = -t182 + (pkin(5) - t120) * t116 + (-t168 - t109 + (-rSges(7,3) - pkin(4) - pkin(8)) * t115) * t114 + t149;
t14 = (-t114 * t33 + t116 * t34) * t115;
t13 = t116 * t130 + (-t116 * rSges(5,3) + t114 * t147) * t114 + t180;
t10 = t116 * t129 + (-t116 * rSges(6,1) + (-rSges(6,2) * t115 + rSges(6,3) * t113) * t114) * t114 + t153;
t9 = t28 * t166 + t30 * t75 + t32 * t76;
t8 = t27 * t166 + t29 * t75 + t31 * t76;
t7 = t28 * t163 + t30 * t73 + t32 * t74;
t6 = t27 * t163 + t29 * t73 + t31 * t74;
t5 = t190 * t116 + (-t116 * pkin(5) + pkin(8) * t166 + t34) * t114 + t153;
t4 = t114 * t8 - t116 * t9;
t3 = t114 * t6 - t116 * t7;
t2 = t113 * t16 + (t114 * t9 + t116 * t8) * t115;
t1 = t113 * t15 + (t114 * t7 + t116 * t6) * t115;
t41 = [t122 * (Icges(4,1) * t122 + t173) + t125 * (Icges(4,2) * t125 + t174) + Icges(2,3) + Icges(3,3) + t43 + m(7) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t77 ^ 2 + t78 ^ 2) + m(2) * (t101 ^ 2 + t102 ^ 2) + (t142 + t170 + t172 + (Icges(5,2) + Icges(6,3)) * t115) * t115 + (t169 + t171 + (Icges(5,1) + Icges(6,2)) * t113) * t113; 0; m(3) + m(4) + m(5) - t185; m(7) * (t17 * t25 + t18 * t24) + m(5) * (t35 * t57 + t36 * t56) + m(6) * (t22 * t38 + t23 * t37) + (-t39 * t183 - t122 * (-Icges(4,5) * t116 + t114 * t139) / 0.2e1 - t125 * (-Icges(4,6) * t116 + t114 * t137) / 0.2e1 + t151 * t116 + (Icges(6,5) * t196 + Icges(5,6) * t195 + t131 * t184 + t136 * t197) * t115 - t156) * t116 + (-t40 * t183 + (Icges(4,5) * t114 + t116 * t139) * t192 + (Icges(4,6) * t114 + t116 * t137) * t191 + t151 * t114 + (Icges(6,5) * t197 + Icges(5,6) * t184 + t131 * t196 + t136 * t195) * t115 + t157) * t114 + ((Icges(6,4) * t196 + Icges(5,5) * t195 + t132 * t184 + t138 * t197) * t116 + (Icges(6,4) * t197 + Icges(5,5) * t184 + t132 * t196 + t138 * t195) * t114) * t113; m(4) * t26 + m(5) * t13 + m(6) * t10 + m(7) * t5; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t10 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(5) * (t13 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(4) * (t159 * t100 ^ 2 + t26 ^ 2) + (t189 * t112 - t4) * t116 + (t3 + t188 * t111 + (t189 * t114 + t188 * t116) * t116) * t114; m(7) * (t114 * t17 - t116 * t18) + m(5) * (t114 * t35 - t116 * t36) + m(6) * (t114 * t22 - t116 * t23); 0; m(7) * (t114 * t25 - t116 * t24) + m(6) * (t114 * t38 - t116 * t37) + m(5) * (t114 * t57 - t116 * t56); 0.2e1 * (m(5) / 0.2e1 + t158) * t159; 0.2e1 * ((t114 * t18 + t116 * t17) * t186 + (t114 * t23 + t116 * t22) * t187) * t113; t185 * t115; m(7) * (-t115 * t5 + (t114 * t24 + t116 * t25) * t113) + m(6) * (-t10 * t115 + (t114 * t37 + t116 * t38) * t113); 0; 0.2e1 * t158 * (t159 * t113 ^ 2 + t115 ^ 2); m(7) * (t17 * t20 + t18 * t21) + t19 + (t156 * t114 + t157 * t116) * t115; m(7) * t14; m(7) * (t14 * t5 + t20 * t25 + t21 * t24) + t1 * t184 + t2 * t196 + (t11 * t114 - t12 * t116) * t193 + (t4 * t184 + t3 * t195) * t115; m(7) * (t114 * t20 - t116 * t21); m(7) * (-t115 * t14 + (t114 * t21 + t116 * t20) * t113); t113 * t19 + m(7) * (t14 ^ 2 + t20 ^ 2 + t21 ^ 2) + (t116 * t1 + t114 * t2 + t113 * (t11 * t116 + t114 * t12)) * t115;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t41(1) t41(2) t41(4) t41(7) t41(11) t41(16); t41(2) t41(3) t41(5) t41(8) t41(12) t41(17); t41(4) t41(5) t41(6) t41(9) t41(13) t41(18); t41(7) t41(8) t41(9) t41(10) t41(14) t41(19); t41(11) t41(12) t41(13) t41(14) t41(15) t41(20); t41(16) t41(17) t41(18) t41(19) t41(20) t41(21);];
Mq  = res;
