% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPP3
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:59
% EndTime: 2019-12-31 18:13:01
% DurationCPUTime: 1.44s
% Computational Cost: add. (1757->275), mult. (4145->291), div. (0->0), fcn. (2986->8), ass. (0->149)
t116 = sin(qJ(1));
t117 = cos(qJ(1));
t158 = g(1) * t116 - g(2) * t117;
t179 = qJDD(1) * pkin(1);
t139 = -qJDD(2) + t179 + t158;
t209 = (qJDD(3) * qJ(4) + qJD(3) * qJD(4));
t108 = pkin(7) + qJ(3);
t101 = sin(t108);
t102 = cos(t108);
t146 = g(1) * t117 + g(2) * t116;
t115 = sin(qJ(3));
t193 = cos(qJ(3));
t159 = qJD(3) * t193;
t166 = qJD(3) * t115;
t111 = sin(pkin(7));
t165 = qJD(1) * qJD(2);
t186 = pkin(6) + qJ(2);
t199 = t186 * qJDD(1) + t165;
t49 = t199 * t111;
t112 = cos(pkin(7));
t50 = t199 * t112;
t77 = t186 * t111;
t73 = qJD(1) * t77;
t78 = t186 * t112;
t74 = qJD(1) * t78;
t156 = t115 * t49 + t73 * t159 + t74 * t166 - t193 * t50;
t124 = -g(3) * t101 - t102 * t146 - t156;
t41 = t115 * t74 + t193 * t73;
t208 = -t41 * qJD(3) - t124;
t188 = pkin(3) + qJ(5);
t154 = t188 * qJDD(3);
t157 = qJDD(1) * t193;
t164 = t111 * qJDD(1);
t144 = -t112 * t157 + t115 * t164;
t72 = t193 * t111 + t115 * t112;
t68 = t72 * qJD(3);
t38 = qJD(1) * t68 + t144;
t161 = t193 * t112;
t147 = qJD(1) * t161;
t173 = t115 * t111;
t160 = qJD(1) * t173;
t63 = -t147 + t160;
t207 = t38 * qJ(5) + t63 * qJD(5);
t184 = t102 * pkin(3) + t101 * qJ(4);
t206 = t72 * qJD(1);
t205 = 2 * t209;
t197 = t206 ^ 2;
t198 = t63 ^ 2;
t203 = -t198 - t197;
t202 = -t198 + t197;
t201 = qJ(2) * qJDD(1);
t200 = -t38 * pkin(4) + qJDD(5);
t28 = (qJD(2) * t111 + qJD(3) * t78) * t115 - qJD(2) * t161 + t77 * t159;
t163 = t112 * qJDD(1);
t162 = qJD(3) * t147 + t111 * t157 + t115 * t163;
t37 = qJD(3) * t160 - t162;
t196 = t37 * pkin(4);
t194 = t63 * pkin(4);
t189 = t63 * t206;
t187 = pkin(4) + t186;
t183 = qJ(4) * t38;
t182 = t63 * qJ(4);
t181 = qJ(5) * t102;
t42 = -t115 * t73 + t193 * t74;
t180 = qJD(3) * t42;
t178 = qJDD(3) * pkin(3);
t177 = t101 * t116;
t176 = t101 * t117;
t175 = t102 * t116;
t174 = t102 * t117;
t172 = t117 * t186;
t170 = -qJD(4) - t41;
t25 = -pkin(4) * t206 - t41;
t169 = qJD(4) - t25;
t26 = t42 - t194;
t168 = -qJD(5) - t26;
t106 = t111 ^ 2;
t107 = t112 ^ 2;
t167 = t106 + t107;
t99 = t112 * pkin(2) + pkin(1);
t43 = t115 * t78 + t193 * t77;
t155 = t115 * t50 + t74 * t159 - t73 * t166 + t193 * t49;
t153 = t167 * qJD(1) ^ 2;
t152 = qJDD(3) - t189;
t82 = t117 * t99;
t151 = g(2) * (pkin(3) * t174 + qJ(4) * t176 + t82);
t150 = 0.2e1 * t167;
t149 = -g(1) * t177 + g(2) * t176;
t148 = g(1) * t175 - g(2) * t174;
t67 = t111 * t166 - t112 * t159;
t10 = -t206 * t67 - t37 * t72;
t71 = -t161 + t173;
t11 = t38 * t71 + t63 * t68;
t145 = -qJDD(4) - t155;
t143 = t67 * qJ(4) - t72 * qJD(4);
t141 = qJD(3) * t67 - t72 * qJDD(3);
t40 = qJD(3) * t68 + qJDD(3) * t71;
t5 = t156 - t209;
t140 = -t72 * qJ(4) - t99;
t138 = -t99 - t184;
t44 = -t115 * t77 + t193 * t78;
t136 = g(1) * t176 + g(2) * t177 - g(3) * t102 - t155;
t76 = -t99 * qJD(1) + qJD(2);
t75 = -t99 * qJDD(1) + qJDD(2);
t135 = -qJDD(4) + t136;
t133 = t139 + t179;
t33 = -qJD(3) * qJ(4) - t42;
t132 = -qJ(4) * t206 + t76;
t131 = -t206 * t68 + t37 * t71 - t38 * t72 + t63 * t67;
t130 = qJD(3) * t28 - qJDD(3) * t44 + t149;
t29 = qJD(2) * t72 + qJD(3) * t44;
t129 = -qJD(3) * t29 - qJDD(3) * t43 + t148;
t27 = t63 * pkin(3) + t132;
t128 = t206 * t27 - t135;
t127 = t38 * pkin(3) + t37 * qJ(4) + t75;
t126 = t150 * t165 - t146;
t12 = t188 * t63 + t132;
t125 = t12 * t206 - t135 - t196;
t123 = t206 * t29 + t28 * t63 - t37 * t43 - t38 * t44 - t146;
t4 = -qJD(4) * t206 + t127;
t122 = t127 - t158;
t21 = 0.2e1 * t206 * qJD(3) + t144;
t121 = -t12 * t63 + t124 + t200;
t118 = qJD(3) ^ 2;
t81 = qJ(4) * t174;
t79 = qJ(4) * t175;
t53 = qJD(3) * t63;
t51 = -t197 - t118;
t36 = pkin(3) * t71 + t140;
t35 = pkin(3) * t206 + t182;
t32 = -qJD(3) * pkin(3) - t170;
t31 = -t71 * pkin(4) + t44;
t30 = pkin(4) * t72 + t43;
t24 = t188 * t71 + t140;
t23 = pkin(3) * t68 + t143;
t20 = (t63 - t160) * qJD(3) + t162;
t19 = t37 - t53;
t18 = (t63 + t160) * qJD(3) - t162;
t17 = t188 * t206 + t182;
t16 = qJD(5) - t33 - t194;
t15 = -t188 * qJD(3) + t169;
t14 = -t67 * pkin(4) + t29;
t13 = -t68 * pkin(4) - t28;
t7 = t71 * qJD(5) + t188 * t68 + t143;
t6 = -t145 - t178;
t3 = -t5 + t200;
t2 = -qJD(3) * qJD(5) - t145 - t154 - t196;
t1 = t4 + t207;
t8 = [0, 0, 0, 0, 0, qJDD(1), t158, t146, 0, 0, t106 * qJDD(1), 0.2e1 * t111 * t163, 0, t107 * qJDD(1), 0, 0, t133 * t112, -t133 * t111, t150 * t201 + t126, t139 * pkin(1) + (t167 * t201 + t126) * qJ(2), t10, t131, -t141, t11, -t40, 0, -t38 * t99 + t68 * t76 + t71 * t75 + t129, t37 * t99 - t67 * t76 + t72 * t75 + t130, t155 * t72 + t156 * t71 - t41 * t67 - t42 * t68 + t123, -t156 * t44 - t42 * t28 + t155 * t43 + t41 * t29 - t75 * t99 - g(1) * (-t116 * t99 + t172) - g(2) * (t116 * t186 + t82), 0, t141, t40, t10, t131, t11, -t32 * t67 + t33 * t68 + t5 * t71 + t6 * t72 + t123, -t23 * t63 - t27 * t68 - t36 * t38 - t4 * t71 - t129, -t206 * t23 + t27 * t67 + t36 * t37 - t4 * t72 - t130, t4 * t36 + t27 * t23 - t5 * t44 + t33 * t28 + t6 * t43 + t32 * t29 - g(1) * t172 - t151 + (-g(1) * t138 - g(2) * t186) * t116, 0, t40, -t141, t11, -t131, t10, -t13 * t63 + t14 * t206 - t15 * t67 - t16 * t68 + t2 * t72 - t3 * t71 - t30 * t37 - t31 * t38 - t146, qJD(3) * t13 + qJDD(3) * t31 - t1 * t72 + t12 * t67 - t206 * t7 + t24 * t37 - t149, -qJD(3) * t14 - qJDD(3) * t30 + t1 * t71 + t12 * t68 + t24 * t38 + t63 * t7 + t148, t1 * t24 + t12 * t7 + t2 * t30 + t15 * t14 + t3 * t31 + t16 * t13 - t151 + (-g(1) * t187 - g(2) * t181) * t117 + (-g(1) * (t138 - t181) - g(2) * t187) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163, t164, -t153, -qJ(2) * t153 - t139, 0, 0, 0, 0, 0, 0, t21, -t18, t203, -t206 * t41 + t42 * t63 - t158 + t75, 0, 0, 0, 0, 0, 0, t203, -t21, t37 + t53, -t33 * t63 + (-qJD(4) - t32) * t206 + t122, 0, 0, 0, 0, 0, 0, t203, t18, t21, t16 * t63 + (-qJD(4) - t15) * t206 + t122 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t202, t20, -t189, -t144, qJDD(3), -t206 * t76 + t136 + t180, t76 * t63 + t208, 0, 0, qJDD(3), t19, t144, t189, t202, -t189, pkin(3) * t37 - t183 + (-t33 - t42) * t206 + (t32 + t170) * t63, t35 * t63 + t128 - 0.2e1 * t178 - t180, t206 * t35 - t27 * t63 + t205 - t208, -t5 * qJ(4) - t6 * pkin(3) - t27 * t35 - t32 * t42 - g(1) * (-pkin(3) * t176 + t81) - g(2) * (-pkin(3) * t177 + t79) - g(3) * t184 + t170 * t33, qJDD(3), t144, t20, -t189, -t202, t189, -t183 + t188 * t37 + (t16 + t168) * t206 + (t15 - t169) * t63, -t25 * qJD(3) + t17 * t206 + t121 + t205, -t17 * t63 + (0.2e1 * qJD(5) + t26) * qJD(3) + 0.2e1 * t154 - t125, t3 * qJ(4) - t12 * t17 - g(1) * t81 - g(2) * t79 - g(3) * (t181 + t184) + t169 * t16 + t168 * t15 + (t146 * t101 - t2) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t152, t51, qJD(3) * t33 + t128 - t178, 0, 0, 0, 0, 0, 0, t20, t51, -t152, -t154 + (-qJD(5) - t16) * qJD(3) + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, qJDD(3) + t189, -t198 - t118, t15 * qJD(3) + t121 + t209;];
tau_reg = t8;
