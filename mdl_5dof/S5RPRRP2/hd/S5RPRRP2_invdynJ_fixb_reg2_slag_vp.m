% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:35
% EndTime: 2020-01-03 11:45:39
% DurationCPUTime: 1.54s
% Computational Cost: add. (1782->249), mult. (3137->293), div. (0->0), fcn. (1761->12), ass. (0->152)
t102 = qJDD(1) + qJDD(3);
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t108 = sin(pkin(8));
t171 = pkin(1) * qJDD(1);
t155 = t108 * t171;
t109 = cos(pkin(8));
t87 = pkin(1) * t109 + pkin(2);
t64 = t87 * qJD(1);
t196 = qJD(3) * t64 + t155;
t193 = pkin(1) * t108;
t156 = qJD(1) * t193;
t197 = qJD(3) * t156 - t87 * qJDD(1);
t147 = -t112 * t197 + t115 * t196;
t13 = pkin(7) * t102 + t147;
t198 = qJD(2) * qJD(4) + t13;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t103 = qJD(1) + qJD(3);
t37 = t112 * t64 + t115 * t156;
t29 = pkin(7) * t103 + t37;
t96 = t114 * qJD(2);
t21 = -t111 * t29 + t96;
t159 = t111 * qJD(2);
t177 = t114 * t29;
t22 = t159 + t177;
t161 = qJD(4) * t111;
t5 = t111 * qJDD(2) + t114 * t198 - t29 * t161;
t93 = t114 * qJDD(2);
t6 = -t22 * qJD(4) - t111 * t13 + t93;
t120 = -t6 * t111 + t5 * t114 + (-t111 * t22 - t114 * t21) * qJD(4);
t104 = qJ(1) + pkin(8);
t94 = qJ(3) + t104;
t86 = cos(t94);
t77 = g(3) * t86;
t85 = sin(t94);
t78 = g(2) * t85;
t183 = t77 - t78;
t154 = t103 * t161;
t189 = t114 * pkin(4);
t89 = pkin(3) + t189;
t179 = t102 * t89;
t195 = -pkin(4) * t154 + t179;
t44 = -t112 * t193 + t115 * t87;
t145 = qJ(5) * t103 + t29;
t131 = t145 * t114;
t194 = g(2) * t86;
t105 = t111 ^ 2;
t192 = pkin(4) * t105;
t191 = g(1) * t114;
t190 = t102 * pkin(3);
t110 = -qJ(5) - pkin(7);
t16 = -t111 * t145 + t96;
t180 = qJD(4) * pkin(4);
t15 = t16 + t180;
t188 = t15 - t16;
t165 = t103 * t114;
t36 = -t112 * t156 + t115 * t64;
t187 = t161 * t36 + t165 * t37;
t186 = t110 * t86 + t85 * t89;
t178 = t111 * t85;
t185 = g(3) * t178 + t111 * t194;
t184 = pkin(3) * t86 + pkin(7) * t85;
t45 = t112 * t87 + t115 * t193;
t113 = sin(qJ(1));
t90 = sin(t104);
t182 = pkin(1) * t113 + pkin(2) * t90;
t116 = cos(qJ(1));
t91 = cos(t104);
t181 = pkin(1) * t116 + pkin(2) * t91;
t175 = t37 * t103;
t38 = t44 * qJD(3);
t174 = t38 * t103;
t39 = t45 * qJD(3);
t173 = t39 * t103;
t43 = pkin(7) + t45;
t172 = -qJ(5) - t43;
t170 = qJ(5) * t102;
t168 = qJDD(4) * pkin(4);
t101 = t103 ^ 2;
t167 = t101 * t114;
t166 = t103 * t111;
t81 = t111 * t102;
t82 = t114 * t102;
t106 = t114 ^ 2;
t164 = t105 - t106;
t163 = t105 + t106;
t162 = qJD(4) * t103;
t160 = qJD(4) * t114;
t95 = t114 * qJD(5);
t3 = t103 * t95 + (-t154 + t82) * qJ(5) + t5;
t157 = t114 * t3 + t183;
t153 = t103 * t160;
t20 = -t103 * t89 + qJD(5) - t36;
t146 = t112 * t196 + t115 * t197;
t8 = qJDD(5) + t146 - t195;
t152 = t111 * t8 + t160 * t20 + t185;
t151 = -t110 * t85 + t86 * t89;
t150 = qJD(4) * t110;
t149 = t163 * t36;
t14 = t146 - t190;
t28 = -pkin(3) * t103 - t36;
t148 = t111 * t14 + t160 * t28 + t185;
t144 = qJD(4) * t172;
t143 = t163 * t102;
t42 = -pkin(3) - t44;
t140 = g(2) * t178 - t191 + t93;
t139 = t111 * t153;
t138 = g(3) * t85 + t194;
t117 = qJD(4) ^ 2;
t137 = pkin(7) * t117 - t190;
t136 = -g(2) * t116 - g(3) * t113;
t135 = -t103 * t28 - t77;
t34 = pkin(4) * t161 + t39;
t35 = t42 - t189;
t134 = t102 * t35 + t103 * t34;
t17 = t159 + t131;
t133 = t111 * t15 - t114 * t17;
t132 = t111 * t21 - t114 * t22;
t130 = -t138 - t8;
t129 = -t147 - t183;
t128 = -t138 - t14;
t127 = g(1) * t111 + t114 * t78 - t5;
t126 = -pkin(3) * t162 - pkin(7) * qJDD(4);
t125 = t102 * t42 + t117 * t43 + t173;
t124 = t138 + t146;
t123 = -qJDD(4) * t43 + (t103 * t42 - t38) * qJD(4);
t122 = -t77 - t170 + (-qJD(5) - t20) * t103;
t119 = t183 + t120;
t107 = qJDD(2) - g(1);
t97 = t114 * qJ(5);
t75 = t85 * pkin(3);
t68 = t111 * t167;
t67 = pkin(7) * t114 + t97;
t66 = t110 * t111;
t56 = qJDD(4) * t114 - t111 * t117;
t55 = qJDD(4) * t111 + t114 * t117;
t48 = t164 * t101;
t47 = -t111 * qJD(5) + t114 * t150;
t46 = t111 * t150 + t95;
t41 = t102 * t106 - 0.2e1 * t139;
t40 = t102 * t105 + 0.2e1 * t139;
t33 = t114 * t43 + t97;
t32 = t172 * t111;
t31 = t36 * t160;
t26 = 0.2e1 * t111 * t82 - 0.2e1 * t162 * t164;
t23 = t28 * t161;
t18 = t20 * t161;
t10 = (-qJD(5) - t38) * t111 + t114 * t144;
t9 = t111 * t144 + t114 * t38 + t95;
t2 = t168 + t93 - qJD(4) * t131 + (-qJD(5) * t103 - t170 - t198) * t111;
t1 = [0, 0, 0, 0, 0, qJDD(1), t136, g(2) * t113 - g(3) * t116, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -g(2) * t91 - g(3) * t90 + 0.2e1 * t109 * t171, g(2) * t90 - g(3) * t91 - 0.2e1 * t155, 0, (t136 + (t108 ^ 2 + t109 ^ 2) * t171) * pkin(1), 0, 0, 0, 0, 0, t102, t102 * t44 - t124 - t173, -t102 * t45 + t129 - t174, 0, -g(2) * t181 - g(3) * t182 - t146 * t44 + t147 * t45 - t36 * t39 + t37 * t38, t40, t26, t55, t41, t56, 0, t23 + t123 * t111 + (-t125 + t128) * t114, t111 * t125 + t114 * t123 + t148, t143 * t43 + t163 * t174 + t119, t14 * t42 + t28 * t39 - g(2) * (t181 + t184) - g(3) * (-pkin(7) * t86 + t182 + t75) - t132 * t38 + t120 * t43, t40, t26, t55, t41, t56, 0, t32 * qJDD(4) + t18 + (t166 * t35 + t10) * qJD(4) + (t130 - t134) * t114, -t33 * qJDD(4) + t134 * t111 + (t165 * t35 - t9) * qJD(4) + t152, (t102 * t33 + t103 * t9 + (-t103 * t32 - t15) * qJD(4)) * t114 + (-t10 * t103 - t102 * t32 - t2 + (-t103 * t33 - t17) * qJD(4)) * t111 + t157, t3 * t33 + t17 * t9 + t2 * t32 + t15 * t10 + t8 * t35 + t20 * t34 - g(2) * (t151 + t181) - g(3) * (t182 + t186); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, t56, -t55, 0, -qJD(4) * t132 + t5 * t111 + t6 * t114 - g(1), 0, 0, 0, 0, 0, 0, t56, -t55, 0, -qJD(4) * t133 + t3 * t111 + t2 * t114 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t124 + t175, t103 * t36 + t129, 0, 0, t40, t26, t55, t41, t56, 0, t23 + t126 * t111 + (t128 - t137) * t114 + t187, t31 + t126 * t114 + (t137 - t175) * t111 + t148, pkin(7) * t143 - t103 * t149 + t119, -t14 * pkin(3) - t28 * t37 - g(2) * t184 - g(3) * t75 + t132 * t36 + (t120 + t77) * pkin(7), t40, t26, t55, t41, t56, 0, t66 * qJDD(4) + t18 + (-t166 * t89 + t47) * qJD(4) + (t130 + t195) * t114 + t187, -t67 * qJDD(4) + t31 + (-t175 - t179) * t111 + (-t46 + (-t114 * t89 + t192) * t103) * qJD(4) + t152, (-qJD(4) * t15 + t102 * t67) * t114 + (-t17 * qJD(4) - t102 * t66 - t2) * t111 + (-t111 * t47 + t114 * t46 - t149 + (-t111 * t67 - t114 * t66) * qJD(4)) * t103 + t157, t3 * t67 + t2 * t66 + t15 * t47 - t8 * t89 - t20 * t37 - g(2) * t151 - g(3) * t186 + (-t114 * t36 + t46) * t17 + (t15 * t36 + t180 * t20) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t48, t81, t68, t82, qJDD(4), (t22 - t177) * qJD(4) + (t135 - t198) * t111 + t140, t21 * qJD(4) + t114 * t135 + t127, 0, 0, -t68, t48, t81, t68, t82, qJDD(4), 0.2e1 * t168 + (t17 - t131) * qJD(4) + (pkin(4) * t167 + t122 - t198) * t111 + t140, -t101 * t192 + (qJ(5) * t166 + t16) * qJD(4) + t122 * t114 + t127, -pkin(4) * t81 + (-t180 + t188) * t165, t188 * t17 + (-t191 + t2 + (-t103 * t20 - t183) * t111) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 + 0.2e1 * t154, t81 + 0.2e1 * t153, -t163 * t101, t103 * t133 - t130;];
tau_reg = t1;
