% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:40:04
% EndTime: 2019-12-05 17:40:08
% DurationCPUTime: 1.72s
% Computational Cost: add. (3129->268), mult. (6141->325), div. (0->0), fcn. (4374->12), ass. (0->150)
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t173 = qJD(5) * t134;
t174 = qJD(5) * t131;
t132 = sin(qJ(4));
t135 = cos(qJ(4));
t128 = cos(pkin(8));
t170 = t128 * qJDD(1);
t127 = sin(pkin(8));
t171 = t127 * qJDD(1);
t153 = -t132 * t171 + t135 * t170;
t85 = t135 * t127 + t132 * t128;
t78 = t85 * qJD(4);
t43 = qJD(1) * t78 - t153;
t177 = qJD(1) * t127;
t168 = t132 * t177;
t146 = -qJD(4) * t168 + t85 * qJDD(1);
t181 = t135 * t128;
t167 = qJD(1) * t181;
t44 = qJD(4) * t167 + t146;
t74 = t85 * qJD(1);
t76 = t167 - t168;
t12 = t131 * t44 + t134 * t43 + t74 * t173 + t76 * t174;
t123 = qJD(4) + qJD(5);
t36 = t131 * t76 + t134 * t74;
t184 = t36 * t123;
t211 = -t12 + t184;
t193 = t36 ^ 2;
t152 = -t131 * t74 + t134 * t76;
t194 = t152 ^ 2;
t210 = -t193 + t194;
t192 = t36 * t152;
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t204 = g(1) * t133 - g(2) * t136;
t13 = t152 * qJD(5) - t131 * t43 + t134 * t44;
t185 = t152 * t123;
t209 = -t13 + t185;
t86 = -t132 * t127 + t181;
t130 = -pkin(1) - qJ(3);
t96 = t130 * qJD(1) + qJD(2);
t165 = -pkin(6) * qJD(1) + t96;
t67 = t165 * t127;
t187 = t132 * t67;
t68 = t165 * t128;
t33 = t135 * t68 - t187;
t26 = -t76 * pkin(7) + t33;
t25 = qJD(4) * pkin(4) + t26;
t34 = t132 * t68 + t135 * t67;
t27 = -t74 * pkin(7) + t34;
t201 = -qJD(1) * qJD(3) + qJDD(1) * t130;
t87 = qJDD(2) + t201;
t159 = -pkin(6) * qJDD(1) + t87;
t56 = t159 * t127;
t57 = t159 * t128;
t161 = -t132 * t56 + t135 * t57;
t15 = -t34 * qJD(4) + t161;
t6 = qJDD(4) * pkin(4) + t43 * pkin(7) + t15;
t175 = qJD(4) * t135;
t169 = -t132 * t57 - t135 * t56 - t68 * t175;
t176 = qJD(4) * t132;
t14 = -t67 * t176 - t169;
t7 = -t44 * pkin(7) + t14;
t1 = (qJD(5) * t25 + t7) * t134 + t131 * t6 - t27 * t174;
t122 = pkin(8) + qJ(4);
t111 = qJ(5) + t122;
t104 = sin(t111);
t105 = cos(t111);
t108 = qJD(1) * qJ(2) + qJD(3);
t91 = pkin(3) * t177 + t108;
t51 = t74 * pkin(4) + t91;
t208 = g(3) * t105 + t204 * t104 + t51 * t36 - t1;
t150 = -t131 * t85 + t134 * t86;
t207 = t12 * t150;
t120 = t127 ^ 2;
t121 = t128 ^ 2;
t178 = t120 + t121;
t206 = t178 * t96;
t157 = g(1) * t136 + g(2) * t133;
t124 = qJDD(1) * qJ(2);
t125 = qJD(1) * qJD(2);
t203 = t124 + t125;
t94 = qJDD(3) + t203;
t147 = -t157 + t94;
t118 = qJDD(4) + qJDD(5);
t205 = t118 * t150;
t186 = t134 * t27;
t9 = t131 * t25 + t186;
t2 = -t9 * qJD(5) - t131 * t7 + t134 * t6;
t202 = g(3) * t104 - t204 * t105 - t51 * t152 + t2;
t79 = -t127 * t176 + t128 * t175;
t140 = t150 * qJD(5) - t131 * t78 + t134 * t79;
t151 = t131 * t86 + t134 * t85;
t200 = t13 * t151 + t140 * t36;
t199 = -t118 * t151 - t123 * t140;
t198 = t1 * t151 + t140 * t9 + t150 * t2 - t204;
t197 = t76 ^ 2;
t196 = 0.2e1 * t125;
t195 = t44 * pkin(4);
t112 = t127 * pkin(3);
t191 = t76 * t74;
t129 = -pkin(6) - qJ(3);
t190 = -pkin(6) + t130;
t189 = -t78 * qJD(4) + t86 * qJDD(4);
t88 = t190 * t127;
t89 = t190 * t128;
t50 = t132 * t89 + t135 * t88;
t188 = t131 * t27;
t183 = pkin(1) * qJDD(1);
t103 = qJ(2) + t112;
t180 = t136 * pkin(1) + t133 * qJ(2);
t162 = -t131 * t79 - t134 * t78;
t49 = -t132 * t88 + t135 * t89;
t160 = t178 * t87;
t158 = qJDD(2) - t183;
t106 = pkin(3) * t171;
t81 = t106 + t94;
t155 = -t86 * t43 - t76 * t78;
t154 = t85 * t44 + t79 * t74;
t31 = -t86 * pkin(7) + t49;
t32 = -t85 * pkin(7) + t50;
t16 = -t131 * t32 + t134 * t31;
t17 = t131 * t31 + t134 * t32;
t149 = -t79 * qJD(4) - t85 * qJDD(4);
t145 = t106 + t147;
t109 = sin(t122);
t110 = cos(t122);
t144 = g(3) * t109 - t110 * t204;
t142 = t203 + t147;
t141 = t14 * t85 + t15 * t86 - t33 * t78 + t34 * t79 - t204;
t30 = -t86 * qJD(3) - t50 * qJD(4);
t29 = -t85 * qJD(3) + t89 * t175 - t88 * t176;
t137 = qJD(1) ^ 2;
t119 = -pkin(7) + t129;
t114 = t136 * qJ(2);
t90 = pkin(4) * t109 + t112;
t71 = t74 ^ 2;
t61 = t79 * pkin(4) + qJD(2);
t59 = t85 * pkin(4) + t103;
t28 = t81 + t195;
t24 = t78 * pkin(7) + t30;
t23 = -t79 * pkin(7) + t29;
t20 = -t151 * qJD(5) + t162;
t19 = t85 * t173 + t86 * t174 - t162;
t11 = t134 * t26 - t188;
t10 = -t131 * t26 - t186;
t8 = t134 * t25 - t188;
t4 = -t17 * qJD(5) - t131 * t23 + t134 * t24;
t3 = t16 * qJD(5) + t131 * t24 + t134 * t23;
t5 = [0, 0, 0, 0, 0, qJDD(1), t204, t157, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t204 - 0.2e1 * t183, 0.2e1 * t124 + t196 - t157, -t158 * pkin(1) - g(1) * (-t133 * pkin(1) + t114) - g(2) * t180 + (t124 + t196) * qJ(2), t121 * qJDD(1), -0.2e1 * t127 * t170, 0, t120 * qJDD(1), 0, 0, t142 * t127, t142 * t128, t204 + t178 * (-t201 - t87), t94 * qJ(2) + t108 * qJD(2) - g(1) * (t130 * t133 + t114) - g(2) * (t136 * qJ(3) + t180) + t130 * t160 - qJD(3) * t206, t155, t43 * t85 - t86 * t44 + t78 * t74 - t76 * t79, t189, t154, t149, 0, qJD(2) * t74 + t30 * qJD(4) + t49 * qJDD(4) + t103 * t44 - t157 * t109 + t91 * t79 + t81 * t85, qJD(2) * t76 - t29 * qJD(4) - t50 * qJDD(4) - t103 * t43 - t157 * t110 - t91 * t78 + t81 * t86, -t29 * t74 - t30 * t76 + t49 * t43 - t50 * t44 - t141, t14 * t50 + t34 * t29 + t15 * t49 + t33 * t30 + t81 * t103 + t91 * qJD(2) - g(1) * (t136 * t112 + t114 + (-pkin(1) + t129) * t133) - g(2) * (t133 * t112 - t136 * t129 + t180), -t152 * t19 - t207, t12 * t151 - t13 * t150 - t140 * t152 + t19 * t36, -t19 * t123 + t205, t200, t199, 0, -t157 * t104 + t16 * t118 + t4 * t123 + t59 * t13 + t140 * t51 + t151 * t28 + t61 * t36, -t105 * t157 - t17 * t118 - t59 * t12 - t3 * t123 + t150 * t28 + t152 * t61 - t51 * t19, t16 * t12 - t17 * t13 - t152 * t4 + t8 * t19 - t3 * t36 - t198, t1 * t17 + t9 * t3 + t2 * t16 + t8 * t4 + t28 * t59 + t51 * t61 - g(1) * (t136 * t90 + t114 + (-pkin(1) + t119) * t133) - g(2) * (-t136 * t119 + t133 * t90 + t180); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t137, -t137 * qJ(2) + t158 - t204, 0, 0, 0, 0, 0, 0, -t137 * t127, -t137 * t128, -t178 * qJDD(1), -t108 * qJD(1) + t160 - t204, 0, 0, 0, 0, 0, 0, -qJD(1) * t74 + t189, -qJD(1) * t76 + t149, -t154 - t155, -t91 * qJD(1) + t141, 0, 0, 0, 0, 0, 0, -qJD(1) * t36 + t20 * t123 + t205, -qJD(1) * t152 + t199, -t152 * t20 - t200 + t207, -t51 * qJD(1) + t8 * t20 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t170, -t178 * t137, qJD(1) * t206 + t147, 0, 0, 0, 0, 0, 0, (t76 + t167) * qJD(4) + t146, -0.2e1 * t74 * qJD(4) + t153, -t71 - t197, t33 * t76 + t34 * t74 + t145, 0, 0, 0, 0, 0, 0, t13 + t185, -t12 - t184, -t193 - t194, t152 * t8 + t36 * t9 + t145 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t71 + t197, t153, -t191, (t76 - t167) * qJD(4) - t146, qJDD(4), -t91 * t76 + t144 + t161, g(3) * t110 + t91 * t74 + t204 * t109 + (t33 + t187) * qJD(4) + t169, 0, 0, t192, t210, t211, -t192, t209, t118, -t10 * t123 + (t118 * t134 - t123 * t174 - t36 * t76) * pkin(4) + t202, t11 * t123 + (-t118 * t131 - t123 * t173 - t152 * t76) * pkin(4) + t208, t10 * t152 + t11 * t36 + t9 * t152 - t8 * t36 + (t12 * t134 - t13 * t131 + (t131 * t152 - t134 * t36) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t131 + t134 * t2 - t51 * t76 + (-t131 * t8 + t134 * t9) * qJD(5) + t144) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t210, t211, -t192, t209, t118, t9 * t123 + t202, t8 * t123 + t208, 0, 0;];
tau_reg = t5;
