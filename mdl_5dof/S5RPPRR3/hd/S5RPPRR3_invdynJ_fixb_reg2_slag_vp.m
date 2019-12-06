% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:42:12
% EndTime: 2019-12-05 17:42:17
% DurationCPUTime: 1.74s
% Computational Cost: add. (3153->267), mult. (6878->334), div. (0->0), fcn. (5142->16), ass. (0->150)
t140 = sin(qJ(5));
t200 = cos(qJ(5));
t143 = cos(qJ(4));
t137 = cos(pkin(9));
t181 = qJD(1) * t137;
t169 = t143 * t181;
t135 = sin(pkin(9));
t141 = sin(qJ(4));
t184 = t141 * t135;
t170 = qJD(1) * t184;
t85 = -t169 + t170;
t96 = t143 * t135 + t141 * t137;
t87 = t96 * qJD(1);
t154 = t140 * t85 - t200 * t87;
t174 = t137 * qJDD(1);
t175 = t135 * qJDD(1);
t171 = qJD(4) * t169 + t141 * t174 + t143 * t175;
t48 = qJD(4) * t170 - t171;
t159 = t141 * t175 - t143 * t174;
t90 = t96 * qJD(4);
t49 = qJD(1) * t90 + t159;
t150 = t154 * qJD(5) + t140 * t48 - t200 * t49;
t133 = qJD(4) + qJD(5);
t187 = t154 * t133;
t211 = t150 - t187;
t168 = qJD(5) * t200;
t178 = qJD(5) * t140;
t153 = -t140 * t49 - t85 * t168 - t87 * t178 - t200 * t48;
t40 = -t140 * t87 - t200 * t85;
t186 = t40 * t133;
t210 = t153 - t186;
t195 = t154 ^ 2;
t196 = t40 ^ 2;
t209 = t195 - t196;
t194 = t40 * t154;
t134 = qJ(1) + pkin(8);
t123 = sin(t134);
t125 = cos(t134);
t206 = -g(2) * t123 + g(3) * t125;
t136 = sin(pkin(8));
t111 = t136 * pkin(1) + qJ(3);
t102 = t111 * qJD(1);
t70 = t135 * qJD(2) + t137 * t102;
t64 = pkin(6) * t181 + t70;
t188 = t141 * t64;
t121 = t137 * qJD(2);
t63 = t121 + (-pkin(6) * qJD(1) - t102) * t135;
t29 = t143 * t63 - t188;
t23 = -t87 * pkin(7) + t29;
t22 = qJD(4) * pkin(4) + t23;
t30 = t141 * t63 + t143 * t64;
t24 = -t85 * pkin(7) + t30;
t119 = t137 * qJDD(2);
t93 = qJD(1) * qJD(3) + t111 * qJDD(1);
t59 = t119 + (-pkin(6) * qJDD(1) - t93) * t135;
t66 = t135 * qJDD(2) + t137 * t93;
t60 = pkin(6) * t174 + t66;
t165 = -t141 * t60 + t143 * t59;
t14 = -t30 * qJD(4) + t165;
t6 = qJDD(4) * pkin(4) + t48 * pkin(7) + t14;
t179 = qJD(4) * t143;
t173 = -t141 * t59 - t143 * t60 - t63 * t179;
t180 = qJD(4) * t141;
t13 = -t64 * t180 - t173;
t7 = -t49 * pkin(7) + t13;
t1 = t140 * t6 + t22 * t168 - t24 * t178 + t200 * t7;
t132 = pkin(9) + qJ(4);
t126 = qJ(5) + t132;
t114 = sin(t126);
t115 = cos(t126);
t138 = cos(pkin(8));
t117 = -t138 * pkin(1) - pkin(2);
t127 = t137 * pkin(3);
t205 = t117 - t127;
t82 = qJD(1) * t205 + qJD(3);
t47 = t85 * pkin(4) + t82;
t208 = g(1) * t114 + t206 * t115 - t47 * t40 - t1;
t172 = t200 * t24;
t9 = t140 * t22 + t172;
t2 = -t9 * qJD(5) - t140 * t7 + t200 * t6;
t207 = -g(1) * t115 + t206 * t114 + t47 * t154 + t2;
t185 = pkin(1) * qJDD(1);
t176 = qJDD(1) * t117;
t100 = qJDD(3) + t176;
t164 = g(2) * t125 + g(3) * t123;
t204 = t164 - t100;
t203 = t87 ^ 2;
t202 = pkin(4) * t90;
t201 = t49 * pkin(4);
t142 = sin(qJ(1));
t198 = t142 * pkin(1);
t144 = cos(qJ(1));
t197 = t144 * pkin(1);
t193 = t87 * t85;
t139 = -pkin(6) - qJ(3);
t192 = pkin(6) + t111;
t116 = t127 + pkin(2);
t89 = t135 * t180 - t137 * t179;
t183 = t143 * t137;
t95 = -t183 + t184;
t25 = t140 * t90 + t95 * t168 + t96 * t178 + t200 * t89;
t53 = -t140 * t95 + t200 * t96;
t191 = t150 * t53 - t25 * t40;
t190 = -t96 * t49 + t89 * t85;
t91 = t192 * t135;
t92 = t192 * t137;
t43 = -t141 * t91 + t143 * t92;
t189 = t140 * t24;
t130 = t135 ^ 2;
t131 = t137 ^ 2;
t182 = t130 + t131;
t42 = -t141 * t92 - t143 * t91;
t162 = g(2) * t144 + g(3) * t142;
t26 = t53 * qJD(5) - t140 * t89 + t200 * t90;
t52 = t140 * t96 + t200 * t95;
t161 = t153 * t52 - t154 * t26;
t160 = -t95 * t48 + t87 * t90;
t128 = qJDD(4) + qJDD(5);
t158 = t53 * t128 - t25 * t133;
t65 = -t135 * t93 + t119;
t157 = -t65 * t135 + t66 * t137;
t156 = (-t135 * t102 + t121) * t135 - t70 * t137;
t34 = -t96 * pkin(7) + t42;
t35 = -t95 * pkin(7) + t43;
t19 = -t140 * t35 + t200 * t34;
t20 = t140 * t34 + t200 * t35;
t152 = -t176 + t204;
t80 = qJDD(1) * t205 + qJDD(3);
t122 = sin(t132);
t124 = cos(t132);
t151 = -g(1) * t124 + t122 * t206;
t31 = qJD(3) * t183 - t91 * t179 + (-qJD(3) * t135 - qJD(4) * t92) * t141;
t149 = -t164 + t80;
t32 = -t96 * qJD(3) - t43 * qJD(4);
t129 = -pkin(7) + t139;
t97 = pkin(4) * t124 + t116;
t83 = t85 ^ 2;
t62 = t95 * pkin(4) + t205;
t51 = -t90 * qJD(4) - t95 * qJDD(4);
t50 = -t89 * qJD(4) + t96 * qJDD(4);
t33 = t80 + t201;
t28 = t89 * pkin(7) + t32;
t27 = -t90 * pkin(7) + t31;
t15 = -t52 * t128 - t26 * t133;
t11 = t200 * t23 - t189;
t10 = -t140 * t23 - t172;
t8 = t200 * t22 - t189;
t4 = -t20 * qJD(5) - t140 * t27 + t200 * t28;
t3 = t19 * qJD(5) + t140 * t28 + t200 * t27;
t5 = [0, 0, 0, 0, 0, qJDD(1), t162, -g(2) * t142 + g(3) * t144, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t138 * t185 + t164, -0.2e1 * t136 * t185 + t206, 0, (t162 + (t136 ^ 2 + t138 ^ 2) * t185) * pkin(1), t130 * qJDD(1), 0.2e1 * t135 * t174, 0, t131 * qJDD(1), 0, 0, t152 * t137, -t152 * t135, t93 * t182 + t157 - t206, t100 * t117 - g(2) * (-t125 * pkin(2) - t123 * qJ(3) - t197) - g(3) * (-t123 * pkin(2) + t125 * qJ(3) - t198) + t157 * t111 - t156 * qJD(3), -t48 * t96 - t87 * t89, -t160 + t190, t50, t49 * t95 + t85 * t90, t51, 0, t32 * qJD(4) + t42 * qJDD(4) + t124 * t164 + t205 * t49 + t80 * t95 + t82 * t90, -t31 * qJD(4) - t43 * qJDD(4) - t122 * t164 - t205 * t48 + t80 * t96 - t82 * t89, -t13 * t95 - t14 * t96 + t29 * t89 - t30 * t90 - t31 * t85 - t32 * t87 + t42 * t48 - t43 * t49 - t206, t13 * t43 + t30 * t31 + t14 * t42 + t29 * t32 + t80 * t205 - g(2) * (-t125 * t116 + t123 * t139 - t197) - g(3) * (-t123 * t116 - t125 * t139 - t198), t153 * t53 + t154 * t25, -t161 + t191, t158, -t150 * t52 - t26 * t40, t15, 0, t115 * t164 + t19 * t128 + t4 * t133 - t150 * t62 - t202 * t40 + t47 * t26 + t33 * t52, -t114 * t164 - t20 * t128 - t3 * t133 + t153 * t62 - t154 * t202 - t47 * t25 + t33 * t53, -t1 * t52 + t150 * t20 - t153 * t19 + t154 * t4 - t2 * t53 + t8 * t25 - t9 * t26 + t3 * t40 - t206, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 + t33 * t62 + t47 * t202 - g(2) * (t123 * t129 - t125 * t97 - t197) - g(3) * (-t123 * t97 - t125 * t129 - t198); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t135 + t65 * t137 - g(1), 0, 0, 0, 0, 0, 0, t51, -t50, t160 + t190, t13 * t96 - t14 * t95 - t29 * t90 - t30 * t89 - g(1), 0, 0, 0, 0, 0, 0, t15, -t158, t161 + t191, t1 * t53 - t2 * t52 - t9 * t25 - t8 * t26 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t175, -t182 * qJD(1) ^ 2, qJD(1) * t156 - t204, 0, 0, 0, 0, 0, 0, 0.2e1 * t87 * qJD(4) + t159, (-t85 - t170) * qJD(4) + t171, -t83 - t203, t29 * t87 + t30 * t85 + t149, 0, 0, 0, 0, 0, 0, -t150 - t187, t153 + t186, -t195 - t196, -t154 * t8 - t9 * t40 + t149 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t83 + t203, (t85 - t170) * qJD(4) + t171, -t193, -t159, qJDD(4), -t82 * t87 + t151 + t165, g(1) * t122 + t82 * t85 + t206 * t124 + (t29 + t188) * qJD(4) + t173, 0, 0, t194, t209, t210, -t194, t211, t128, -t10 * t133 + (t200 * t128 - t133 * t178 + t40 * t87) * pkin(4) + t207, t11 * t133 + (-t128 * t140 - t133 * t168 + t154 * t87) * pkin(4) + t208, -t10 * t154 - t11 * t40 - t9 * t154 + t8 * t40 + (-t200 * t153 + t140 * t150 + (-t140 * t154 + t200 * t40) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t200 * t2 + t1 * t140 - t47 * t87 + (-t140 * t8 + t200 * t9) * qJD(5) + t151) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t209, t210, -t194, t211, t128, t9 * t133 + t207, t8 * t133 + t208, 0, 0;];
tau_reg = t5;
