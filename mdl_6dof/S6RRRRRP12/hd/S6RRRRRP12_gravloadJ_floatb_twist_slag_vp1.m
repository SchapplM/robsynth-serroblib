% Calculate Gravitation load on the joints for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:41
% EndTime: 2019-03-10 02:59:47
% DurationCPUTime: 2.33s
% Computational Cost: add. (1773->266), mult. (4828->395), div. (0->0), fcn. (6195->14), ass. (0->120)
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t111 = sin(pkin(7));
t117 = cos(qJ(1));
t163 = sin(pkin(6));
t164 = cos(pkin(7));
t141 = t164 * t163;
t165 = cos(pkin(6));
t178 = cos(qJ(2));
t151 = t165 * t178;
t175 = sin(qJ(2));
t176 = sin(qJ(1));
t94 = -t117 * t151 + t176 * t175;
t191 = -t94 * t111 + t117 * t141;
t114 = sin(qJ(3));
t154 = t117 * t163;
t155 = t114 * t164;
t177 = cos(qJ(3));
t150 = t165 * t175;
t95 = t117 * t150 + t176 * t178;
t54 = -t111 * t114 * t154 - t94 * t155 + t95 * t177;
t22 = -t113 * t191 + t116 * t54;
t147 = t177 * t163;
t149 = t164 * t177;
t53 = t117 * t111 * t147 + t95 * t114 + t94 * t149;
t1 = t112 * t22 - t53 * t115;
t2 = t112 * t53 + t115 * t22;
t125 = t117 * t175 + t176 * t151;
t119 = t125 * t111 + t176 * t141;
t190 = -t113 * t54 - t116 * t191;
t182 = rSges(7,2) + pkin(12);
t146 = t163 * t176;
t189 = -t111 * t146 + t125 * t164;
t188 = t165 * t111 + t178 * t141;
t183 = rSges(7,1) + pkin(5);
t180 = rSges(5,3) + pkin(11);
t179 = rSges(6,3) + pkin(12);
t174 = pkin(4) * t116;
t173 = pkin(10) * t111;
t172 = t113 * t53;
t96 = t117 * t178 - t176 * t150;
t57 = t114 * t96 + t177 * t189;
t171 = t113 * t57;
t145 = t163 * t175;
t77 = t114 * t145 - t177 * t188;
t170 = t113 * t77;
t166 = rSges(7,3) + qJ(6);
t162 = t111 * t113;
t161 = t111 * t116;
t160 = t112 * t116;
t159 = t115 * t116;
t136 = t111 * t145;
t148 = t178 * t163;
t158 = pkin(2) * t148 + pkin(10) * t136;
t157 = t117 * pkin(1) + pkin(9) * t146;
t131 = t175 * t141;
t88 = -t114 * t131 + t178 * t147;
t156 = t88 * pkin(3) + t158;
t152 = -t176 * pkin(1) + pkin(9) * t154;
t47 = t53 * pkin(3);
t144 = pkin(11) * t54 - pkin(12) * t172 - t53 * t174 - t47;
t49 = t57 * pkin(3);
t58 = -t114 * t189 + t96 * t177;
t143 = pkin(11) * t58 - pkin(12) * t171 - t57 * t174 - t49;
t76 = t77 * pkin(3);
t78 = t114 * t188 + t177 * t145;
t142 = pkin(11) * t78 - pkin(12) * t170 - t77 * t174 - t76;
t140 = -rSges(5,1) * t116 + rSges(5,2) * t113;
t139 = rSges(6,1) * t115 - rSges(6,2) * t112;
t64 = -t95 * t155 - t94 * t177;
t89 = t94 * pkin(2);
t138 = t64 * pkin(3) + t95 * t173 - t89;
t66 = -t125 * t177 - t96 * t155;
t91 = t125 * pkin(2);
t137 = t66 * pkin(3) + t96 * t173 - t91;
t69 = t113 * t136 + t88 * t116;
t87 = t114 * t148 + t177 * t131;
t133 = t69 * pkin(4) + pkin(11) * t87 + t156;
t32 = t116 * t64 + t95 * t162;
t63 = -t114 * t94 + t95 * t149;
t129 = t32 * pkin(4) + pkin(11) * t63 + t138;
t34 = t116 * t66 + t96 * t162;
t65 = -t125 * t114 + t96 * t149;
t128 = t34 * pkin(4) + pkin(11) * t65 + t137;
t127 = -t95 * pkin(2) + t191 * pkin(10) + t152;
t126 = -pkin(3) * t54 + t127;
t123 = -pkin(4) * t22 - pkin(11) * t53 + t126;
t121 = t96 * pkin(2) + t119 * pkin(10) + t157;
t120 = t58 * pkin(3) + t121;
t26 = t119 * t113 + t58 * t116;
t118 = t26 * pkin(4) + t57 * pkin(11) + t120;
t93 = -t111 * t148 + t165 * t164;
t68 = t113 * t88 - t116 * t136;
t52 = t113 * t93 + t116 * t78;
t51 = -t113 * t78 + t116 * t93;
t46 = t51 * pkin(4);
t36 = t112 * t87 + t115 * t69;
t35 = t112 * t69 - t87 * t115;
t33 = t113 * t66 - t96 * t161;
t31 = t113 * t64 - t95 * t161;
t28 = t112 * t78 - t77 * t159;
t27 = -t78 * t115 - t77 * t160;
t25 = t113 * t58 - t119 * t116;
t19 = t25 * pkin(4);
t17 = t190 * pkin(4);
t16 = t112 * t77 + t115 * t52;
t15 = t112 * t52 - t77 * t115;
t14 = t112 * t65 + t115 * t34;
t13 = t112 * t34 - t65 * t115;
t12 = t112 * t63 + t115 * t32;
t11 = t112 * t32 - t63 * t115;
t10 = t112 * t58 - t57 * t159;
t9 = -t58 * t115 - t57 * t160;
t8 = t112 * t54 - t53 * t159;
t7 = -t54 * t115 - t53 * t160;
t6 = t112 * t57 + t115 * t26;
t5 = t112 * t26 - t57 * t115;
t3 = [-m(2) * (g(1) * (-t176 * rSges(2,1) - t117 * rSges(2,2)) + g(2) * (t117 * rSges(2,1) - t176 * rSges(2,2))) - m(3) * (g(1) * (-t95 * rSges(3,1) + t94 * rSges(3,2) + rSges(3,3) * t154 + t152) + g(2) * (t96 * rSges(3,1) - t125 * rSges(3,2) + rSges(3,3) * t146 + t157)) - m(4) * (g(1) * (-rSges(4,1) * t54 + rSges(4,2) * t53 + rSges(4,3) * t191 + t127) + g(2) * (t58 * rSges(4,1) - t57 * rSges(4,2) + t119 * rSges(4,3) + t121)) - m(5) * (g(1) * (-rSges(5,1) * t22 - rSges(5,2) * t190 - t180 * t53 + t126) + g(2) * (t26 * rSges(5,1) - t25 * rSges(5,2) + t180 * t57 + t120)) - m(6) * (g(1) * (-rSges(6,1) * t2 + rSges(6,2) * t1 + t179 * t190 + t123) + g(2) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t179 * t25 + t118)) - m(7) * (g(1) * (-t1 * t166 + t182 * t190 - t183 * t2 + t123) + g(2) * (t166 * t5 + t182 * t25 + t183 * t6 + t118)) -m(3) * (g(1) * (-t125 * rSges(3,1) - t96 * rSges(3,2)) + g(2) * (-rSges(3,1) * t94 - rSges(3,2) * t95) + g(3) * (rSges(3,1) * t148 - rSges(3,2) * t145)) - m(4) * (g(1) * (rSges(4,1) * t66 - rSges(4,2) * t65 - t91) + g(2) * (rSges(4,1) * t64 - rSges(4,2) * t63 - t89) + g(3) * (t88 * rSges(4,1) - t87 * rSges(4,2) + t158) + (g(3) * rSges(4,3) * t145 + (g(1) * t96 + g(2) * t95) * (rSges(4,3) + pkin(10))) * t111) - m(5) * (g(1) * (rSges(5,1) * t34 - rSges(5,2) * t33 + t180 * t65 + t137) + g(2) * (rSges(5,1) * t32 - rSges(5,2) * t31 + t180 * t63 + t138) + g(3) * (rSges(5,1) * t69 - rSges(5,2) * t68 + t180 * t87 + t156)) - m(6) * (g(1) * (rSges(6,1) * t14 - rSges(6,2) * t13 + t179 * t33 + t128) + g(2) * (rSges(6,1) * t12 - rSges(6,2) * t11 + t179 * t31 + t129) + g(3) * (rSges(6,1) * t36 - rSges(6,2) * t35 + t179 * t68 + t133)) - m(7) * (g(1) * (t166 * t13 + t183 * t14 + t182 * t33 + t128) + g(2) * (t166 * t11 + t183 * t12 + t182 * t31 + t129) + g(3) * (t166 * t35 + t182 * t68 + t183 * t36 + t133)) -m(4) * (g(1) * (-rSges(4,1) * t57 - rSges(4,2) * t58) + g(2) * (-rSges(4,1) * t53 - rSges(4,2) * t54) + g(3) * (-rSges(4,1) * t77 - rSges(4,2) * t78)) - m(5) * (g(1) * (t140 * t57 + t180 * t58 - t49) + g(2) * (t140 * t53 + t180 * t54 - t47) + g(3) * (t140 * t77 + t180 * t78 - t76)) - m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t9 - rSges(6,3) * t171 + t143) + g(2) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t172 + t144) + g(3) * (rSges(6,1) * t28 - rSges(6,2) * t27 - rSges(6,3) * t170 + t142)) - m(7) * (g(1) * (-rSges(7,2) * t171 + t183 * t10 + t166 * t9 + t143) + g(2) * (-rSges(7,2) * t172 + t166 * t7 + t183 * t8 + t144) + g(3) * (-rSges(7,2) * t170 + t166 * t27 + t183 * t28 + t142)) -m(5) * (g(1) * (-rSges(5,1) * t25 - rSges(5,2) * t26) + g(2) * (rSges(5,1) * t190 - rSges(5,2) * t22) + g(3) * (rSges(5,1) * t51 - rSges(5,2) * t52)) - m(6) * (g(1) * (-t139 * t25 + t179 * t26 - t19) + g(2) * (t139 * t190 + t179 * t22 + t17) + g(3) * (t139 * t51 + t179 * t52 + t46)) + (-g(1) * (t182 * t26 - t19) - g(2) * (t182 * t22 + t17) - g(3) * (t182 * t52 + t46) + (-t112 * t166 - t115 * t183) * (-g(1) * t25 + g(2) * t190 + g(3) * t51)) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t5 - rSges(6,2) * t6) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t15 - rSges(6,2) * t16)) - m(7) * (g(1) * (t166 * t6 - t183 * t5) + g(2) * (-t183 * t1 + t166 * t2) + g(3) * (-t183 * t15 + t166 * t16)) -m(7) * (g(1) * t5 + g(2) * t1 + g(3) * t15)];
taug  = t3(:);
