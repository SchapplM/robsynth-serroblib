% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:18
% EndTime: 2019-12-05 16:44:21
% DurationCPUTime: 1.59s
% Computational Cost: add. (2183->287), mult. (4840->332), div. (0->0), fcn. (3293->8), ass. (0->161)
t128 = sin(qJ(4));
t129 = sin(qJ(3));
t130 = cos(qJ(3));
t203 = cos(qJ(4));
t73 = t128 * t130 + t203 * t129;
t66 = t73 * qJD(2);
t190 = t66 * qJ(5);
t179 = qJD(1) * t129;
t131 = -pkin(7) - pkin(6);
t93 = t131 * t130;
t61 = -qJD(2) * t93 + t179;
t55 = t128 * t61;
t194 = qJD(3) * pkin(3);
t117 = t130 * qJD(1);
t164 = qJD(2) * t131;
t60 = t129 * t164 + t117;
t58 = t60 + t194;
t32 = t203 * t58 - t55;
t16 = t32 - t190;
t123 = qJD(3) + qJD(4);
t120 = qJDD(3) + qJDD(4);
t160 = qJD(4) * t203;
t153 = pkin(3) * t160;
t202 = pkin(3) * t128;
t207 = -t120 * t202 - t123 * t153;
t112 = t120 * pkin(4);
t182 = t128 * t129;
t145 = t123 * t182;
t162 = t203 * t130;
t150 = qJD(2) * t162;
t155 = qJDD(2) * t203;
t173 = t130 * qJDD(2);
t154 = -t123 * t150 - t128 * t173 - t129 * t155;
t26 = qJD(2) * t145 + t154;
t192 = t26 * qJ(5);
t206 = t112 + t192;
t122 = pkin(8) + qJ(2);
t114 = cos(t122);
t127 = qJ(3) + qJ(4);
t118 = sin(t127);
t184 = t114 * t118;
t113 = sin(t122);
t186 = t113 * t118;
t119 = cos(t127);
t201 = g(3) * t119;
t205 = g(1) * t184 + g(2) * t186 - t201;
t204 = t66 ^ 2;
t200 = g(3) * t130;
t178 = qJD(2) * t129;
t161 = t128 * t178;
t64 = -t150 + t161;
t199 = t66 * t64;
t14 = pkin(4) * t123 + t16;
t198 = t14 - t16;
t174 = t129 * qJDD(2);
t144 = t128 * t174 - t130 * t155;
t43 = t123 * t73;
t27 = qJD(2) * t43 + t144;
t197 = -t64 * t153 - t27 * t202;
t42 = -qJD(3) * t162 - t130 * t160 + t145;
t196 = -t73 * t27 + t42 * t64;
t35 = t203 * t60 - t55;
t92 = t131 * t129;
t46 = t128 * t92 - t203 * t93;
t195 = qJ(5) * t27;
t193 = t123 * t64;
t191 = t64 * qJ(5);
t188 = pkin(6) * qJDD(2);
t187 = qJDD(2) * pkin(2);
t185 = t113 * t119;
t183 = t114 * t119;
t157 = pkin(4) * t64 + qJD(5);
t110 = pkin(3) * t130 + pkin(2);
t91 = t110 * qJD(2);
t44 = t157 - t91;
t181 = qJD(5) + t44;
t124 = t129 ^ 2;
t125 = t130 ^ 2;
t180 = t124 - t125;
t177 = qJD(4) * t128;
t176 = qJD(1) * qJD(3);
t175 = qJD(2) * qJD(3);
t172 = t203 * pkin(3);
t171 = pkin(3) * t178;
t170 = t129 * t194;
t169 = pkin(3) * t177;
t168 = pkin(6) * t178;
t167 = pkin(6) * qJD(2) * t130;
t57 = t203 * t61;
t133 = qJD(2) ^ 2;
t166 = t129 * t133 * t130;
t165 = pkin(6) * t173 + t129 * qJDD(1) + t130 * t176;
t163 = qJD(3) * t131;
t159 = t129 * t175;
t158 = t130 * t175;
t116 = t130 * qJDD(1);
t37 = qJDD(3) * pkin(3) + t116 + qJDD(2) * t92 + (t130 * t164 - t179) * qJD(3);
t47 = -pkin(6) * t159 + t165;
t41 = (-t159 + t173) * pkin(7) + t47;
t156 = -t128 * t41 + t203 * t37;
t34 = -t128 * t60 - t57;
t45 = t128 * t93 + t203 * t92;
t5 = t128 * t37 + t58 * t160 - t61 * t177 + t203 * t41;
t152 = -g(1) * t186 + g(2) * t184;
t151 = g(1) * t185 - g(2) * t183;
t149 = t129 * t158;
t148 = g(1) * t114 + g(2) * t113;
t147 = g(1) * t113 - g(2) * t114;
t72 = -t162 + t182;
t146 = -t26 * t72 + t43 * t66;
t23 = t120 * t73 - t123 * t42;
t33 = t128 * t58 + t57;
t143 = g(1) * t183 + g(2) * t185 + g(3) * t118 - t5;
t142 = -0.2e1 * pkin(2) * t175 - pkin(6) * qJDD(3);
t78 = t129 * t163;
t79 = t130 * t163;
t18 = t128 * t79 + t92 * t160 + t93 * t177 + t203 * t78;
t62 = pkin(3) * t159 - t110 * qJDD(2);
t141 = pkin(2) * t133 + t148;
t140 = -t64 * t91 + t143;
t132 = qJD(3) ^ 2;
t139 = -pkin(6) * t132 + t147 + 0.2e1 * t187;
t11 = pkin(4) * t27 + qJDD(5) + t62;
t138 = -t123 * t161 - t154;
t6 = -t33 * qJD(4) + t156;
t19 = -t46 * qJD(4) - t128 * t78 + t203 * t79;
t137 = t181 * t64 + t143 + t195;
t136 = t6 + t205;
t48 = -t129 * t176 + t116 + (-t158 - t174) * pkin(6);
t80 = t117 - t168;
t81 = t167 + t179;
t135 = -t48 * t129 + t47 * t130 + (-t129 * t81 - t130 * t80) * qJD(3) - t148;
t134 = t91 * t66 + t136;
t126 = qJDD(1) - g(3);
t121 = -qJ(5) + t131;
t109 = t172 + pkin(4);
t107 = pkin(4) * t119;
t90 = qJDD(3) * t130 - t129 * t132;
t89 = qJDD(3) * t129 + t130 * t132;
t77 = t107 + t110;
t63 = t64 ^ 2;
t54 = pkin(4) * t72 - t110;
t50 = pkin(4) * t66 + t171;
t40 = pkin(4) * t43 + t170;
t39 = -qJ(5) * t72 + t46;
t38 = -qJ(5) * t73 + t45;
t28 = -t63 + t204;
t25 = -t120 * t72 - t123 * t43;
t21 = -t190 + t35;
t20 = t34 + t191;
t17 = t33 - t191;
t12 = t138 + t193;
t10 = t42 * qJ(5) - t73 * qJD(5) + t19;
t9 = -qJ(5) * t43 - qJD(5) * t72 + t18;
t8 = t27 * t72 + t43 * t64;
t7 = -t26 * t73 - t42 * t66;
t4 = -qJD(5) * t64 - t195 + t5;
t3 = -t66 * qJD(5) + t206 + t6;
t2 = t146 + t196;
t1 = -t146 + t196;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, 0, t90, -t89, 0, t129 * t47 + t130 * t48 - g(3) + (-t129 * t80 + t130 * t81) * qJD(3), 0, 0, 0, 0, 0, 0, t25, -t23, t2, -t32 * t43 - t33 * t42 + t5 * t73 - t6 * t72 - g(3), 0, 0, 0, 0, 0, 0, t25, -t23, t2, -t14 * t43 - t17 * t42 - t3 * t72 + t4 * t73 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t147, t148, 0, 0, qJDD(2) * t124 + 0.2e1 * t149, 0.2e1 * t129 * t173 - 0.2e1 * t180 * t175, t89, qJDD(2) * t125 - 0.2e1 * t149, t90, 0, t129 * t142 + t130 * t139, -t129 * t139 + t130 * t142, (t124 + t125) * t188 + t135, (t147 + t187) * pkin(2) + t135 * pkin(6), t7, t1, t23, t8, t25, 0, -t110 * t27 + t120 * t45 + t123 * t19 + t64 * t170 - t43 * t91 + t62 * t72 + t151, t110 * t26 - t120 * t46 - t123 * t18 + t66 * t170 + t42 * t91 + t62 * t73 + t152, -t18 * t64 - t19 * t66 + t26 * t45 - t27 * t46 + t32 * t42 - t33 * t43 - t5 * t72 - t6 * t73 - t148, t5 * t46 + t33 * t18 + t6 * t45 + t32 * t19 - t62 * t110 - t91 * t170 - g(1) * (-t110 * t113 - t114 * t131) - g(2) * (t110 * t114 - t113 * t131), t7, t1, t23, t8, t25, 0, t10 * t123 + t11 * t72 + t120 * t38 + t27 * t54 + t40 * t64 + t43 * t44 + t151, t11 * t73 - t120 * t39 - t123 * t9 - t26 * t54 + t40 * t66 - t42 * t44 + t152, -t10 * t66 + t14 * t42 - t17 * t43 + t26 * t38 - t27 * t39 - t3 * t73 - t4 * t72 - t64 * t9 - t148, t4 * t39 + t17 * t9 + t3 * t38 + t14 * t10 + t11 * t54 + t44 * t40 - g(1) * (-t113 * t77 - t114 * t121) - g(2) * (-t113 * t121 + t114 * t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t180 * t133, t174, t166, t173, qJDD(3), -t200 + t116 + (t81 - t167) * qJD(3) + (t141 - t176 - t188) * t129, g(3) * t129 + (t80 + t168) * qJD(3) + t141 * t130 - t165, 0, 0, t199, t28, t12, -t199, -t144, t120, -t34 * t123 + (t203 * t120 - t123 * t177 - t64 * t178) * pkin(3) + t134, t123 * t35 - t66 * t171 + t140 + t207, t26 * t172 + (-t32 + t35) * t64 + (t33 + t34 + t169) * t66 + t197, -t32 * t34 - t33 * t35 + (t203 * t6 - t200 + t128 * t5 + (-t128 * t32 + t203 * t33) * qJD(4) + (qJD(2) * t91 + t148) * t129) * pkin(3), t199, t28, t12, -t199, -t144, t120, t109 * t120 - t20 * t123 - t50 * t64 - t181 * t66 + (-t57 + (-pkin(3) * t123 - t58) * t128) * qJD(4) + t156 + t205 + t206, t123 * t21 - t50 * t66 + t137 + t207, t109 * t26 + (-t14 + t21) * t64 + (t17 + t20 + t169) * t66 + t197, -g(3) * t107 + t3 * t109 - t14 * t20 - t17 * t21 - t44 * t50 - t148 * (-pkin(3) * t129 - pkin(4) * t118) + (-t200 + t4 * t128 + (-t128 * t14 + t203 * t17) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t28, t12, -t199, -t144, t120, t33 * t123 + t134, t123 * t32 + t140, 0, 0, t199, t28, t12, -t199, -t144, t120, t192 + t17 * t123 + 0.2e1 * t112 + (-t157 - t44) * t66 + t136, -pkin(4) * t204 + t123 * t16 + t137, pkin(4) * t26 - t198 * t64, t198 * t17 + (t118 * t148 - t44 * t66 - t201 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t123 + t27, t138 - t193, -t63 - t204, t14 * t66 + t17 * t64 + t11 - t147;];
tau_reg = t13;
