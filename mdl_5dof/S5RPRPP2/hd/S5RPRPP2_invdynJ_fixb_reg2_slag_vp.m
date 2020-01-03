% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:17
% EndTime: 2019-12-31 18:11:18
% DurationCPUTime: 1.15s
% Computational Cost: add. (982->241), mult. (1874->266), div. (0->0), fcn. (994->8), ass. (0->142)
t99 = cos(qJ(3));
t163 = qJ(4) * t99;
t173 = pkin(3) + pkin(4);
t97 = sin(qJ(3));
t176 = t173 * t97;
t115 = t163 - t176;
t144 = qJD(1) * qJD(3);
t134 = t99 * t144;
t77 = t97 * qJDD(1);
t180 = t134 + t77;
t95 = sin(pkin(7));
t65 = t95 * pkin(1) + pkin(6);
t179 = qJD(2) * qJD(3) + qJDD(1) * t65;
t44 = t65 * qJD(1);
t35 = t97 * t44;
t23 = t99 * qJD(2) - t35;
t178 = qJD(4) - t23;
t136 = t173 * qJD(3);
t150 = qJ(5) * qJD(1);
t12 = t97 * t150 + t23;
t152 = qJD(4) - t12;
t9 = -t136 + t152;
t96 = cos(pkin(7));
t171 = t96 * pkin(1);
t139 = pkin(2) + t171;
t45 = t139 * qJD(1);
t132 = t173 * qJDD(3);
t154 = t97 * qJD(2);
t24 = t99 * t44 + t154;
t91 = qJD(3) * qJ(4);
t19 = t24 + t91;
t88 = qJ(1) + pkin(7);
t74 = sin(t88);
t75 = cos(t88);
t166 = g(1) * t75 + g(2) * t74;
t13 = t154 + (t44 - t150) * t99;
t11 = t91 + t13;
t81 = t97 * qJ(4);
t133 = pkin(2) + t81;
t10 = qJD(5) + (t173 * t99 + t133 + t171) * qJD(1);
t151 = qJD(5) + t10;
t177 = t151 * t97;
t15 = -qJD(3) * pkin(3) + t178;
t161 = pkin(1) * qJDD(1);
t168 = t75 * t97;
t169 = t74 * t97;
t175 = g(1) * t168 + g(2) * t169 - g(3) * t99;
t148 = qJDD(3) * t65;
t85 = t99 * pkin(3);
t116 = -t133 - t85;
t20 = (t116 - t171) * qJD(1);
t165 = t85 + t81;
t27 = -t139 - t165;
t174 = (qJD(1) * t27 + t20) * qJD(3) - t148;
t92 = t97 ^ 2;
t93 = t99 ^ 2;
t164 = -t92 + t93;
t79 = t99 * qJDD(1);
t14 = 0.2e1 * t164 * t144 + 0.2e1 * t97 * t79;
t63 = g(1) * t74;
t172 = g(2) * t75;
t84 = t99 * pkin(4);
t170 = t11 * t99;
t167 = t75 * t99;
t69 = t96 * t161;
t43 = -qJDD(1) * pkin(2) - t69;
t162 = qJ(5) - t65;
t160 = qJD(1) * t97;
t29 = t162 * t99;
t159 = qJD(3) * t29;
t158 = qJD(3) * t44;
t157 = qJD(3) * t97;
t156 = qJDD(3) * pkin(3);
t155 = t11 * qJD(3);
t153 = t97 * qJD(4);
t149 = qJ(5) * qJD(3);
t147 = t92 * qJDD(1);
t146 = t93 * qJDD(1);
t145 = qJ(5) * qJDD(1);
t143 = qJD(1) * qJD(5);
t141 = t97 * qJDD(2) + t179 * t99;
t100 = cos(qJ(1));
t140 = t100 * pkin(1) + t75 * pkin(2) + t74 * pkin(6);
t98 = sin(qJ(1));
t138 = -t98 * pkin(1) + t75 * pkin(6);
t137 = t63 - t172;
t135 = t97 * t144;
t131 = -g(3) - t158;
t130 = (-qJDD(2) + t158) * t99 + t179 * t97;
t129 = -t166 + (t146 + t147) * t65;
t22 = -t27 + t84;
t128 = qJD(1) * t22 + t10;
t125 = t97 * t134;
t89 = qJDD(3) * qJ(4);
t90 = qJD(3) * qJD(4);
t124 = 0.2e1 * t89 + 0.2e1 * t90 + t141;
t123 = pkin(3) * t167 + t75 * t81 + t140;
t121 = g(1) * t98 - g(2) * t100;
t120 = pkin(3) * t97 - t163;
t102 = qJD(3) ^ 2;
t119 = t102 * t65 + t172;
t118 = -qJDD(4) - t130;
t117 = pkin(3) * t79 + t180 * qJ(4) + qJD(1) * t153 - t43;
t114 = -t130 + t175;
t7 = -t44 * t157 + t141;
t113 = pkin(4) * t79 + qJDD(5) + t117;
t5 = -t118 - t156;
t4 = t7 + t89 + t90;
t16 = t115 * qJD(3) + t153;
t3 = -t173 * t135 + t113;
t112 = qJD(1) * t16 + qJDD(1) * t22 - t172 + t3;
t111 = -qJDD(1) * t139 + t119 + t43;
t110 = -0.2e1 * t45 * qJD(3) - t148;
t109 = t24 * qJD(3) + t114;
t26 = t120 * qJD(3) - t153;
t6 = pkin(3) * t135 - t117;
t108 = -qJD(1) * t26 - qJDD(1) * t27 - t119 - t6;
t107 = -t97 * t145 + qJDD(4) - t114;
t106 = t4 * t99 + t5 * t97 + (t15 * t99 - t19 * t97) * qJD(3);
t105 = t7 * t99 + t130 * t97 + (-t23 * t99 - t24 * t97) * qJD(3);
t103 = qJD(1) ^ 2;
t58 = t97 * t103 * t99;
t56 = qJ(5) * t135;
t55 = -t92 * t103 - t102;
t53 = t99 * t63;
t52 = g(1) * t169;
t49 = t75 * t163;
t47 = t74 * t163;
t46 = -qJDD(3) - t58;
t41 = t164 * t103;
t40 = qJDD(3) * t99 - t102 * t97;
t39 = qJDD(3) * t97 + t102 * t99;
t36 = t120 * qJD(1);
t31 = -0.2e1 * t125 + t146;
t30 = 0.2e1 * t125 + t147;
t28 = t162 * t97;
t25 = t115 * qJD(1);
t18 = -t97 * qJD(5) - t159;
t17 = -t99 * qJD(5) + t162 * t157;
t2 = t56 + (-t143 - t145) * t99 + t4;
t1 = -t180 * qJ(5) - t97 * t143 - t118 - t132;
t8 = [0, 0, 0, 0, 0, qJDD(1), t121, g(1) * t100 + g(2) * t98, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t137 + 0.2e1 * t69, -0.2e1 * t95 * t161 + t166, 0, (t121 + (t95 ^ 2 + t96 ^ 2) * t161) * pkin(1), t30, t14, t39, t31, t40, 0, t110 * t97 - t111 * t99 + t53, t110 * t99 + t111 * t97 - t52, t105 + t129, -t43 * t139 - g(1) * (-t74 * pkin(2) + t138) - g(2) * t140 + t105 * t65, t30, t39, -t14, 0, -t40, t31, t108 * t99 + t174 * t97 + t53, t106 + t129, t108 * t97 - t174 * t99 + t52, -g(1) * t138 - g(2) * t123 + t106 * t65 - t116 * t63 + t20 * t26 + t6 * t27, t30, -t14, -t39, t31, t40, 0, t28 * qJDD(3) + t53 + (-t128 * t97 - t18) * qJD(3) + t112 * t99, -t29 * qJDD(3) + t52 + (t128 * t99 + t17) * qJD(3) + t112 * t97, (-qJD(3) * t9 + qJDD(1) * t29 - t2 + (qJD(3) * t28 - t17) * qJD(1)) * t99 + (t155 + qJDD(1) * t28 - t1 + (-t18 - t159) * qJD(1)) * t97 + t166, -t2 * t29 + t11 * t17 - t1 * t28 + t9 * t18 + t3 * t22 + t10 * t16 - g(1) * (-t75 * qJ(5) + t138) - g(2) * (pkin(4) * t167 + t123) + (-g(1) * (t116 - t84) + g(2) * qJ(5)) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t40, -t39, 0, t7 * t97 - t130 * t99 - g(3) + (-t23 * t97 + t24 * t99) * qJD(3), 0, 0, 0, 0, 0, 0, t40, 0, t39, t4 * t97 - t5 * t99 - g(3) + (t15 * t97 + t19 * t99) * qJD(3), 0, 0, 0, 0, 0, 0, t40, t39, 0, -t1 * t99 + t2 * t97 - g(3) + (t9 * t97 + t170) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t41, t77, t58, t79, qJDD(3), t160 * t45 + t109, g(3) * t97 + (t23 + t35) * qJD(3) + (qJD(1) * t45 + t166) * t99 - t141, 0, 0, -t58, t77, t41, qJDD(3), -t79, t58, 0.2e1 * t156 - qJDD(4) + (-t20 * t97 + t36 * t99) * qJD(1) + t109, -t120 * qJDD(1), -t23 * qJD(3) + (qJD(1) * t20 - t166) * t99 + (qJD(1) * t36 + t131) * t97 + t124, t4 * qJ(4) - t5 * pkin(3) - t20 * t36 - t15 * t24 - g(1) * (-pkin(3) * t168 + t49) - g(2) * (-pkin(3) * t169 + t47) - g(3) * t165 + t178 * t19, -t58, t41, -t77, t58, t79, qJDD(3), t13 * qJD(3) + 0.2e1 * t132 + ((-t25 + t149) * t99 + t177) * qJD(1) - t107, -t12 * qJD(3) + t56 + (-qJD(1) * t25 + t131) * t97 + (-qJD(1) * t151 - t145 - t166) * t99 + t124, -t115 * qJDD(1), t2 * qJ(4) - t1 * t173 - t9 * t13 - t10 * t25 - g(1) * t49 - g(2) * t47 - g(3) * (t84 + t165) + t152 * t11 + t166 * t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t77, t55, -t19 * qJD(3) + t160 * t20 - t175 + t5, 0, 0, 0, 0, 0, 0, t46, t55, -t77, -t155 - t132 + (-t149 * t99 - t177) * qJD(1) + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 - 0.2e1 * t135, t77 + 0.2e1 * t134, (-t92 - t93) * t103, (t170 + (t9 - t136) * t97) * qJD(1) + t113 + t137;];
tau_reg = t8;
