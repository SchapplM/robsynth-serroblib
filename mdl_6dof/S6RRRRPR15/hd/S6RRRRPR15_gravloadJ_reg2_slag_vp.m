% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR15_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:39:26
% EndTime: 2019-05-08 03:39:34
% DurationCPUTime: 1.82s
% Computational Cost: add. (1387->211), mult. (3872->312), div. (0->0), fcn. (4968->14), ass. (0->121)
t100 = cos(qJ(6));
t101 = cos(qJ(4));
t102 = cos(qJ(1));
t154 = cos(pkin(7));
t95 = sin(pkin(6));
t145 = t95 * t154;
t155 = cos(pkin(6));
t139 = t102 * t155;
t169 = sin(qJ(1));
t171 = cos(qJ(2));
t99 = sin(qJ(2));
t77 = -t171 * t139 + t169 * t99;
t94 = sin(pkin(7));
t185 = t102 * t145 - t77 * t94;
t98 = sin(qJ(3));
t144 = t98 * t154;
t156 = t102 * t95;
t151 = t94 * t156;
t170 = cos(qJ(3));
t78 = t99 * t139 + t169 * t171;
t40 = -t77 * t144 - t98 * t151 + t78 * t170;
t97 = sin(qJ(4));
t16 = t101 * t185 + t40 * t97;
t129 = t154 * t170;
t39 = t77 * t129 + t170 * t151 + t78 * t98;
t96 = sin(qJ(6));
t190 = t100 * t39 + t16 * t96;
t189 = t100 * t16 - t39 * t96;
t130 = t155 * t169;
t112 = t102 * t99 + t171 * t130;
t186 = t112 * t94 + t169 * t145;
t17 = t40 * t101 - t185 * t97;
t160 = t101 * t39;
t162 = qJ(5) * t97;
t184 = -pkin(4) * t160 - t39 * t162;
t149 = t95 * t169;
t179 = t112 * t154 - t94 * t149;
t79 = t102 * t171 - t99 * t130;
t43 = t179 * t170 + t79 * t98;
t159 = t101 * t43;
t183 = -pkin(4) * t159 - t43 * t162;
t143 = t155 * t94;
t150 = t95 * t171;
t166 = t95 * t99;
t60 = -t129 * t150 - t170 * t143 + t98 * t166;
t158 = t101 * t60;
t182 = -pkin(4) * t158 - t60 * t162;
t178 = pkin(5) + pkin(11);
t177 = pkin(10) * t94;
t47 = t78 * t129 - t77 * t98;
t176 = pkin(11) * t47;
t49 = -t112 * t98 + t79 * t129;
t175 = pkin(11) * t49;
t70 = (t99 * t129 + t171 * t98) * t95;
t174 = pkin(11) * t70;
t173 = t39 * pkin(11);
t172 = t43 * pkin(11);
t167 = t94 * t97;
t165 = t96 * t97;
t153 = t94 * t166;
t164 = pkin(2) * t150 + pkin(10) * t153;
t163 = t102 * pkin(1) + pkin(9) * t149;
t161 = t100 * t97;
t157 = t101 * t94;
t71 = (-t99 * t144 + t171 * t170) * t95;
t152 = t71 * pkin(3) + t164;
t33 = t39 * pkin(3);
t148 = pkin(11) * t40 - t33;
t35 = t43 * pkin(3);
t44 = t79 * t170 - t179 * t98;
t147 = pkin(11) * t44 - t35;
t59 = t60 * pkin(3);
t61 = t98 * t143 + (t171 * t144 + t170 * t99) * t95;
t146 = pkin(11) * t61 - t59;
t142 = -t16 * pkin(4) + qJ(5) * t17;
t20 = -t101 * t186 + t44 * t97;
t21 = t44 * t101 + t186 * t97;
t141 = -t20 * pkin(4) + qJ(5) * t21;
t111 = -t94 * t150 + t155 * t154;
t37 = -t111 * t101 + t61 * t97;
t38 = t61 * t101 + t111 * t97;
t140 = -t37 * pkin(4) + qJ(5) * t38;
t137 = -t77 * pkin(2) + t78 * t177;
t136 = -t112 * pkin(2) + t79 * t177;
t134 = -t169 * pkin(1) + pkin(9) * t156;
t133 = -g(1) * t16 + g(2) * t20;
t132 = -g(1) * t17 + g(2) * t21;
t131 = -g(1) * t39 + g(2) * t43;
t48 = -t78 * t144 - t77 * t170;
t128 = t48 * pkin(3) + t137;
t50 = -t112 * t170 - t79 * t144;
t127 = t50 * pkin(3) + t136;
t125 = -g(1) * t102 - g(2) * t169;
t52 = -t101 * t153 + t71 * t97;
t53 = t101 * t71 + t97 * t153;
t124 = t53 * pkin(4) + qJ(5) * t52 + t152;
t2 = g(1) * t20 + g(2) * t16 + g(3) * t37;
t123 = g(1) * t21 + g(2) * t17 + g(3) * t38;
t24 = -t78 * t157 + t48 * t97;
t26 = -t79 * t157 + t50 * t97;
t122 = g(1) * t26 + g(2) * t24 + g(3) * t52;
t25 = t101 * t48 + t78 * t167;
t27 = t101 * t50 + t79 * t167;
t121 = g(1) * t27 + g(2) * t25 + g(3) * t53;
t120 = g(1) * t43 + g(2) * t39 + g(3) * t60;
t119 = g(1) * t44 + g(2) * t40 + g(3) * t61;
t118 = g(1) * t49 + g(2) * t47 + g(3) * t70;
t116 = g(1) * t79 + g(2) * t78 + g(3) * t166;
t115 = t25 * pkin(4) + qJ(5) * t24 + t128;
t114 = t27 * pkin(4) + qJ(5) * t26 + t127;
t113 = -t78 * pkin(2) + t185 * pkin(10) + t134;
t110 = -pkin(3) * t40 + t113;
t107 = -pkin(4) * t17 - qJ(5) * t16 + t110;
t106 = t79 * pkin(2) + t186 * pkin(10) + t163;
t105 = t44 * pkin(3) + t106;
t103 = t21 * pkin(4) + t20 * qJ(5) + t105;
t8 = t100 * t43 + t20 * t96;
t7 = t100 * t20 - t43 * t96;
t6 = t120 * t101;
t5 = t120 * t97;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t169 - g(2) * t102, -t125, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t78 - g(2) * t79, -g(1) * t77 + g(2) * t112, t125 * t95, -g(1) * t134 - g(2) * t163, 0, 0, 0, 0, 0, 0, g(1) * t40 - g(2) * t44, t131, -g(1) * t185 - g(2) * t186, -g(1) * t113 - g(2) * t106, 0, 0, 0, 0, 0, 0, -t132, t133, -t131, -g(1) * (t110 - t173) - g(2) * (t105 + t172) 0, 0, 0, 0, 0, 0, -t131, t132, -t133, -g(1) * (t107 - t173) - g(2) * (t103 + t172) 0, 0, 0, 0, 0, 0, g(1) * t190 - g(2) * t8, g(1) * t189 - g(2) * t7, -t132, -g(1) * (-pkin(12) * t17 - t178 * t39 + t107) - g(2) * (t21 * pkin(12) + t178 * t43 + t103); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t112 + g(2) * t77 - g(3) * t150, t116, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t48 - g(3) * t71, t118, -t116 * t94, -g(1) * t136 - g(2) * t137 - g(3) * t164, 0, 0, 0, 0, 0, 0, -t121, t122, -t118, -g(1) * (t127 + t175) - g(2) * (t128 + t176) - g(3) * (t152 + t174) 0, 0, 0, 0, 0, 0, -t118, t121, -t122, -g(1) * (t114 + t175) - g(2) * (t115 + t176) - g(3) * (t124 + t174) 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t49 + t26 * t96) - g(2) * (t100 * t47 + t24 * t96) - g(3) * (t100 * t70 + t52 * t96) -g(1) * (t100 * t26 - t49 * t96) - g(2) * (t100 * t24 - t47 * t96) - g(3) * (t100 * t52 - t70 * t96) -t121, -g(1) * (pkin(12) * t27 + t178 * t49 + t114) - g(2) * (pkin(12) * t25 + t178 * t47 + t115) - g(3) * (pkin(12) * t53 + t178 * t70 + t124); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t119, -g(1) * t147 - g(2) * t148 - g(3) * t146, 0, 0, 0, 0, 0, 0, -t119, -t6, t5, -g(1) * (t147 + t183) - g(2) * (t148 + t184) - g(3) * (t146 + t182) 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t44 - t43 * t165) - g(2) * (t100 * t40 - t39 * t165) - g(3) * (t100 * t61 - t60 * t165) -g(1) * (-t43 * t161 - t44 * t96) - g(2) * (-t39 * t161 - t40 * t96) - g(3) * (-t60 * t161 - t61 * t96) t6, -g(1) * (-pkin(12) * t159 + t178 * t44 + t183 - t35) - g(2) * (-pkin(12) * t160 + t178 * t40 + t184 - t33) - g(3) * (-pkin(12) * t158 + t178 * t61 + t182 - t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t123, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t123, -g(1) * t141 - g(2) * t142 - g(3) * t140, 0, 0, 0, 0, 0, 0, -t123 * t96, -t123 * t100, t2, -g(1) * (-pkin(12) * t20 + t141) - g(2) * (-pkin(12) * t16 + t142) - g(3) * (-pkin(12) * t37 + t140); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t189 - g(3) * (t100 * t37 - t60 * t96) g(1) * t8 + g(2) * t190 - g(3) * (-t100 * t60 - t37 * t96) 0, 0;];
taug_reg  = t1;
