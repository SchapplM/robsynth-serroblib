% Calculate joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:31:42
% DurationCPUTime: 1.97s
% Computational Cost: add. (3249->308), mult. (3397->462), div. (0->0), fcn. (3513->10), ass. (0->147)
t198 = Icges(3,3) + Icges(4,3);
t121 = qJ(2) + pkin(8);
t114 = sin(t121);
t116 = cos(t121);
t128 = sin(qJ(2));
t130 = cos(qJ(2));
t197 = Icges(3,5) * t130 + Icges(4,5) * t116 - Icges(3,6) * t128 - Icges(4,6) * t114;
t129 = sin(qJ(1));
t186 = t129 / 0.2e1;
t131 = cos(qJ(1));
t196 = t131 / 0.2e1;
t195 = t128 / 0.2e1;
t194 = t130 / 0.2e1;
t193 = -t197 * t129 + t198 * t131;
t192 = t198 * t129 + t197 * t131;
t125 = cos(pkin(9));
t110 = pkin(4) * t125 + pkin(3);
t127 = -pkin(7) - qJ(4);
t124 = sin(pkin(9));
t164 = t129 * t124;
t170 = t116 * t131;
t171 = t114 * t131;
t120 = pkin(9) + qJ(5);
t113 = sin(t120);
t115 = cos(t120);
t165 = t129 * t115;
t74 = -t113 * t170 + t165;
t166 = t129 * t113;
t75 = t115 * t170 + t166;
t33 = t75 * rSges(6,1) + t74 * rSges(6,2) + rSges(6,3) * t171;
t191 = pkin(4) * t164 + t110 * t170 - t127 * t171 + t33;
t190 = t116 ^ 2;
t122 = t129 ^ 2;
t123 = t131 ^ 2;
t189 = m(5) / 0.2e1;
t188 = m(6) / 0.2e1;
t98 = rSges(3,1) * t128 + rSges(3,2) * t130;
t187 = m(3) * t98;
t185 = -t131 / 0.2e1;
t184 = pkin(2) * t128;
t183 = pkin(3) * t116;
t111 = pkin(2) * t130 + pkin(1);
t106 = t131 * t111;
t119 = t131 * pkin(6);
t126 = -qJ(3) - pkin(6);
t167 = t126 * t131;
t182 = t129 * (t167 + t119 + (-pkin(1) + t111) * t129) + t131 * (-pkin(1) * t131 + t106 + (-pkin(6) - t126) * t129);
t181 = rSges(3,1) * t130;
t180 = rSges(3,2) * t128;
t48 = -Icges(6,6) * t116 + (Icges(6,4) * t115 - Icges(6,2) * t113) * t114;
t179 = t113 * t48;
t178 = t131 * rSges(3,3);
t177 = Icges(3,4) * t128;
t176 = Icges(3,4) * t130;
t175 = Icges(4,4) * t114;
t174 = Icges(4,4) * t116;
t173 = t110 * t116;
t172 = t114 * t129;
t169 = t124 * t131;
t168 = t125 * t131;
t163 = t129 * t125;
t162 = pkin(3) * t170 + qJ(4) * t171;
t161 = t129 * rSges(3,3) + t131 * t181;
t160 = t122 + t123;
t159 = t189 + t188;
t86 = -t116 * t169 + t163;
t87 = t116 * t168 + t164;
t158 = t87 * rSges(5,1) + t86 * rSges(5,2) + rSges(5,3) * t171;
t72 = -t115 * t131 - t116 * t166;
t73 = -t113 * t131 + t116 * t165;
t26 = Icges(6,5) * t73 + Icges(6,6) * t72 + Icges(6,3) * t172;
t28 = Icges(6,4) * t73 + Icges(6,2) * t72 + Icges(6,6) * t172;
t30 = Icges(6,1) * t73 + Icges(6,4) * t72 + Icges(6,5) * t172;
t11 = -t116 * t26 + (-t113 * t28 + t115 * t30) * t114;
t47 = -Icges(6,3) * t116 + (Icges(6,5) * t115 - Icges(6,6) * t113) * t114;
t49 = -Icges(6,5) * t116 + (Icges(6,1) * t115 - Icges(6,4) * t113) * t114;
t13 = t47 * t172 + t48 * t72 + t49 * t73;
t157 = t13 / 0.2e1 + t11 / 0.2e1;
t27 = Icges(6,5) * t75 + Icges(6,6) * t74 + Icges(6,3) * t171;
t29 = Icges(6,4) * t75 + Icges(6,2) * t74 + Icges(6,6) * t171;
t31 = Icges(6,1) * t75 + Icges(6,4) * t74 + Icges(6,5) * t171;
t12 = -t116 * t27 + (-t113 * t29 + t115 * t31) * t114;
t14 = t47 * t171 + t74 * t48 + t75 * t49;
t156 = t14 / 0.2e1 + t12 / 0.2e1;
t155 = Icges(3,5) * t195 + Icges(3,6) * t194 + Icges(4,5) * t114 / 0.2e1 + Icges(4,6) * t116 / 0.2e1;
t154 = -pkin(3) * t114 + qJ(4) * t116 - t184;
t153 = -rSges(4,1) * t114 - rSges(4,2) * t116 - t184;
t145 = qJ(4) * t114 + t183;
t152 = t122 * t145 + t131 * t162 + t182;
t151 = -t129 * t126 + t106;
t150 = t154 + rSges(5,3) * t116 - (rSges(5,1) * t125 - rSges(5,2) * t124) * t114;
t84 = -t116 * t164 - t168;
t85 = t116 * t163 - t169;
t149 = -t85 * rSges(5,1) - t84 * rSges(5,2);
t148 = -t73 * rSges(6,1) - t72 * rSges(6,2);
t147 = -t180 + t181;
t146 = rSges(4,1) * t116 - rSges(4,2) * t114;
t50 = -rSges(6,3) * t116 + (rSges(6,1) * t115 - rSges(6,2) * t113) * t114;
t140 = t154 - (qJ(4) + t127) * t116 - (-pkin(3) + t110) * t114 - t50;
t139 = Icges(3,1) * t130 - t177;
t138 = Icges(4,1) * t116 - t175;
t137 = -Icges(3,2) * t128 + t176;
t136 = -Icges(4,2) * t114 + t174;
t133 = rSges(4,1) * t170 - rSges(4,2) * t171 + t129 * rSges(4,3);
t100 = rSges(2,1) * t131 - t129 * rSges(2,2);
t99 = -t129 * rSges(2,1) - rSges(2,2) * t131;
t60 = t153 * t131;
t59 = t153 * t129;
t57 = -Icges(5,5) * t116 + (Icges(5,1) * t125 - Icges(5,4) * t124) * t114;
t56 = -Icges(5,6) * t116 + (Icges(5,4) * t125 - Icges(5,2) * t124) * t114;
t54 = t129 * pkin(6) + (pkin(1) - t180) * t131 + t161;
t53 = t178 + t119 + (-pkin(1) - t147) * t129;
t45 = t133 + t151;
t44 = (rSges(4,3) - t126) * t131 + (-t111 - t146) * t129;
t43 = t114 * t115 * t49;
t42 = t131 * (-t131 * t180 + t161) + (t147 * t129 - t178) * t129;
t41 = Icges(5,1) * t87 + Icges(5,4) * t86 + Icges(5,5) * t171;
t40 = Icges(5,1) * t85 + Icges(5,4) * t84 + Icges(5,5) * t172;
t39 = Icges(5,4) * t87 + Icges(5,2) * t86 + Icges(5,6) * t171;
t38 = Icges(5,4) * t85 + Icges(5,2) * t84 + Icges(5,6) * t172;
t37 = Icges(5,5) * t87 + Icges(5,6) * t86 + Icges(5,3) * t171;
t36 = Icges(5,5) * t85 + Icges(5,6) * t84 + Icges(5,3) * t172;
t35 = t150 * t131;
t34 = t150 * t129;
t32 = rSges(6,3) * t172 - t148;
t25 = t151 + t158 + t162;
t24 = -t167 + (-t183 - t111 + (-rSges(5,3) - qJ(4)) * t114) * t129 + t149;
t23 = t140 * t131;
t22 = t140 * t129;
t21 = -t116 * t33 - t50 * t171;
t20 = t116 * t32 + t50 * t172;
t19 = t151 + t191;
t18 = (pkin(4) * t124 - t126) * t131 + (-t173 - t111 + (-rSges(6,3) + t127) * t114) * t129 + t148;
t17 = t131 * t133 + (-t131 * rSges(4,3) + t146 * t129) * t129 + t182;
t16 = -t114 * t179 - t116 * t47 + t43;
t15 = (-t129 * t33 + t131 * t32) * t114;
t10 = t129 * (rSges(5,3) * t172 - t149) + t131 * t158 + t152;
t9 = t27 * t171 + t74 * t29 + t75 * t31;
t8 = t26 * t171 + t74 * t28 + t75 * t30;
t7 = t27 * t172 + t29 * t72 + t31 * t73;
t6 = t26 * t172 + t28 * t72 + t30 * t73;
t5 = (-t162 + t191) * t131 + (-pkin(4) * t169 + t32 + (-t114 * t127 - t145 + t173) * t129) * t129 + t152;
t4 = t9 * t129 - t131 * t8;
t3 = t7 * t129 - t131 * t6;
t2 = -t14 * t116 + (t129 * t8 + t131 * t9) * t114;
t1 = -t13 * t116 + (t129 * t6 + t131 * t7) * t114;
t46 = [t128 * (Icges(3,1) * t128 + t176) + t130 * (Icges(3,2) * t130 + t177) + Icges(2,3) + t43 + (-t47 + t175 - (Icges(5,5) * t125 - Icges(5,6) * t124) * t114 + (Icges(4,2) + Icges(5,3)) * t116) * t116 + (Icges(4,1) * t114 - t124 * t56 + t125 * t57 + t174 - t179) * t114 + m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(3) * (t53 ^ 2 + t54 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(2) * (t100 ^ 2 + t99 ^ 2); m(6) * (t18 * t23 + t22 * t19) + m(5) * (t24 * t35 + t25 * t34) + m(4) * (t44 * t60 + t45 * t59) + (-t128 * (-Icges(3,5) * t131 + t139 * t129) / 0.2e1 - t130 * (-Icges(3,6) * t131 + t137 * t129) / 0.2e1 - t53 * t187 - t84 * t56 / 0.2e1 - t85 * t57 / 0.2e1 + t155 * t131 - t157) * t131 + (-t54 * t187 + t86 * t56 / 0.2e1 + t87 * t57 / 0.2e1 + (Icges(3,5) * t129 + t139 * t131) * t195 + (Icges(3,6) * t129 + t137 * t131) * t194 + t155 * t129 + t156) * t129 + ((Icges(4,6) * t196 - t136 * t129 / 0.2e1 + t36 / 0.2e1) * t131 + (Icges(4,6) * t186 + t136 * t196 - t37 / 0.2e1) * t129) * t116 + ((Icges(4,5) * t129 - t124 * t39 + t125 * t41 + t138 * t131) * t186 + (-Icges(4,5) * t131 - t124 * t38 + t125 * t40 + t138 * t129) * t185) * t114; m(6) * (t22 ^ 2 + t23 ^ 2 + t5 ^ 2) + m(5) * (t10 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(4) * (t17 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(3) * (t160 * t98 ^ 2 + t42 ^ 2) + (-t3 + (t36 * t172 + t84 * t38 + t85 * t40) * t131 + t193 * t123) * t131 + (t4 + (t37 * t171 + t86 * t39 + t87 * t41) * t129 + t192 * t122 + (t193 * t129 + t192 * t131 - t36 * t171 - t37 * t172 - t86 * t38 - t84 * t39 - t87 * t40 - t85 * t41) * t131) * t129; m(6) * (t129 * t18 - t131 * t19) + m(5) * (t129 * t24 - t131 * t25) + m(4) * (t129 * t44 - t131 * t45); m(6) * (t129 * t23 - t131 * t22) + m(5) * (t129 * t35 - t131 * t34) + m(4) * (t129 * t60 - t131 * t59); 0.2e1 * (m(4) / 0.2e1 + t159) * t160; 0.2e1 * ((t129 * t19 + t131 * t18) * t188 + (t129 * t25 + t131 * t24) * t189) * t114; m(6) * (-t116 * t5 + (t129 * t22 + t131 * t23) * t114) + m(5) * (-t116 * t10 + (t129 * t34 + t131 * t35) * t114); 0; 0.2e1 * t159 * (t160 * t114 ^ 2 + t190); m(6) * (t18 * t20 + t19 * t21) - t16 * t116 + (t157 * t129 + t156 * t131) * t114; m(6) * (t15 * t5 + t20 * t23 + t21 * t22) + t1 * t185 - t116 * (-t11 * t131 + t12 * t129) / 0.2e1 + t2 * t186 + (t3 * t186 + t4 * t196) * t114; m(6) * (t20 * t129 - t131 * t21); m(6) * (-t15 * t116 + (t129 * t21 + t131 * t20) * t114); m(6) * (t15 ^ 2 + t20 ^ 2 + t21 ^ 2) + t190 * t16 + (t131 * t2 + t129 * t1 - t116 * (t11 * t129 + t12 * t131)) * t114;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t46(1), t46(2), t46(4), t46(7), t46(11); t46(2), t46(3), t46(5), t46(8), t46(12); t46(4), t46(5), t46(6), t46(9), t46(13); t46(7), t46(8), t46(9), t46(10), t46(14); t46(11), t46(12), t46(13), t46(14), t46(15);];
Mq = res;
