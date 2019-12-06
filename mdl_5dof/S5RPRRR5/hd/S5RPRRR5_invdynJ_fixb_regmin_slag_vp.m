% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:31
% EndTime: 2019-12-05 18:16:34
% DurationCPUTime: 0.88s
% Computational Cost: add. (1295->172), mult. (2151->235), div. (0->0), fcn. (1414->14), ass. (0->131)
t189 = pkin(7) + pkin(8);
t92 = qJ(1) + pkin(9) + qJ(3);
t85 = sin(t92);
t86 = cos(t92);
t188 = g(2) * t86 + g(3) * t85;
t101 = qJD(1) + qJD(3);
t107 = sin(qJ(5));
t108 = sin(qJ(4));
t111 = cos(qJ(5));
t112 = cos(qJ(4));
t52 = t107 * t112 + t111 * t108;
t45 = t52 * t101;
t105 = sin(pkin(9));
t179 = pkin(1) * t105;
t152 = qJD(1) * t179;
t106 = cos(pkin(9));
t87 = pkin(1) * t106 + pkin(2);
t187 = -qJD(3) * t152 + qJDD(1) * t87;
t185 = -g(2) * t85 + g(3) * t86;
t71 = t87 * qJD(1);
t186 = qJD(3) * t71 + qJDD(1) * t179;
t160 = t111 * t112;
t164 = t107 * t108;
t51 = -t160 + t164;
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t140 = -t109 * t179 + t113 * t87;
t40 = t109 * t71 + t113 * t152;
t144 = t101 * t189 + t40;
t25 = t112 * qJD(2) - t108 * t144;
t100 = qJD(4) + qJD(5);
t26 = t108 * qJD(2) + t112 * t144;
t184 = t109 * t187 + t113 * t186;
t99 = qJDD(1) + qJDD(3);
t181 = t99 * pkin(3);
t172 = t109 * t87 + t113 * t179;
t48 = pkin(7) + t172;
t180 = -pkin(8) - t48;
t178 = t101 * pkin(3);
t177 = t112 * pkin(4);
t148 = t101 * t160;
t149 = t101 * t164;
t43 = -t148 + t149;
t176 = t45 * t43;
t104 = qJ(4) + qJ(5);
t95 = cos(t104);
t175 = t85 * t95;
t174 = t86 * t95;
t154 = t112 * qJD(4);
t139 = -t109 * t186 + t113 * t187;
t20 = -t139 - t181;
t39 = -t109 * t152 + t113 * t71;
t33 = -t39 - t178;
t173 = t108 * t20 + t154 * t33;
t171 = t111 * t26;
t170 = t112 * t99;
t168 = t40 * t101;
t42 = t172 * qJD(3);
t167 = t42 * t101;
t165 = t101 * t108;
t162 = t108 * t112;
t159 = qJDD(2) - g(1);
t102 = t108 ^ 2;
t158 = -t112 ^ 2 + t102;
t157 = qJD(5) * t107;
t155 = t108 * qJD(4);
t153 = t112 * t188 + t33 * t155;
t151 = pkin(4) * t155;
t90 = -pkin(3) - t177;
t19 = t99 * pkin(7) + t184;
t147 = pkin(8) * t99 + t19;
t146 = qJD(4) * t189;
t145 = t101 * t154;
t143 = qJD(4) * t180;
t10 = (t101 * t155 - t170) * pkin(4) + t20;
t27 = t101 * t90 - t39;
t31 = t100 * t52;
t142 = g(2) * t174 + g(3) * t175 + t10 * t51 + t27 * t31;
t47 = -pkin(3) - t140;
t137 = -t40 + t151;
t135 = t51 * t99;
t110 = sin(qJ(1));
t114 = cos(qJ(1));
t134 = g(2) * t114 + g(3) * t110;
t30 = t100 * t51;
t98 = qJDD(4) + qJDD(5);
t21 = -t100 * t30 + t52 * t98;
t24 = qJD(4) * pkin(4) + t25;
t132 = -t107 * t24 - t171;
t35 = t180 * t108;
t96 = t112 * pkin(8);
t36 = t112 * t48 + t96;
t131 = -t107 * t36 + t111 * t35;
t130 = t107 * t35 + t111 * t36;
t129 = t139 + t188;
t77 = pkin(7) * t112 + t96;
t128 = qJD(5) * t77 - t108 * t39 + t112 * t146;
t76 = t189 * t108;
t127 = qJD(5) * t76 + t108 * t146 + t112 * t39;
t115 = qJD(4) ^ 2;
t125 = -pkin(7) * t115 + t168 + t181;
t124 = -t115 * t48 - t47 * t99 - t167;
t123 = -t101 * t33 + t185 - t19;
t94 = sin(t104);
t122 = t10 * t52 - t188 * t94 - t27 * t30;
t121 = -pkin(7) * qJDD(4) + (t39 - t178) * qJD(4);
t41 = t140 * qJD(3);
t120 = -qJDD(4) * t48 + (t101 * t47 - t41) * qJD(4);
t13 = qJD(5) * t148 - t100 * t149 + t111 * t145 + t52 * t99;
t91 = t112 * qJDD(2);
t4 = qJDD(4) * pkin(4) - qJD(4) * t26 - t108 * t147 + t91;
t119 = -g(2) * t175 + t27 * t43 + t26 * t157 + g(3) * t174 + g(1) * t94 + (-t100 * t26 - t4) * t107;
t118 = -t184 + t185;
t5 = qJD(4) * t25 + t108 * qJDD(2) + t112 * t147;
t117 = -g(1) * t95 + t132 * qJD(5) - t107 * t5 + t111 * t4 + t185 * t94 - t27 * t45;
t97 = t101 ^ 2;
t62 = qJDD(4) * t112 - t108 * t115;
t61 = qJDD(4) * t108 + t112 * t115;
t46 = t102 * t99 + 0.2e1 * t108 * t145;
t38 = t47 - t177;
t37 = t42 + t151;
t32 = -0.2e1 * qJD(4) * t101 * t158 + 0.2e1 * t162 * t99;
t22 = -t100 * t31 - t51 * t98;
t18 = -t108 * t41 + t112 * t143;
t17 = t108 * t143 + t112 * t41;
t15 = -t43 ^ 2 + t45 ^ 2;
t14 = t101 * t31 + t135;
t6 = t43 * t100 + t13;
t2 = t13 * t52 - t30 * t45;
t1 = -t13 * t51 - t14 * t52 + t30 * t43 - t31 * t45;
t3 = [qJDD(1), t134, -g(2) * t110 + g(3) * t114, (t134 + (t105 ^ 2 + t106 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t99, t140 * t99 + t129 - t167, -t41 * t101 - t172 * t99 + t118, t46, t32, t61, t62, 0, t120 * t108 + (t124 - t20) * t112 + t153, t120 * t112 + (-t124 - t188) * t108 + t173, t2, t1, t21, t22, 0, t37 * t43 + t38 * t14 + (-qJD(5) * t130 - t107 * t17 + t111 * t18) * t100 + t131 * t98 + t142, t37 * t45 + t38 * t13 - (qJD(5) * t131 + t107 * t18 + t111 * t17) * t100 - t130 * t98 + t122; 0, 0, 0, t159, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, 0, 0, 0, 0, 0, t22, -t21; 0, 0, 0, 0, t99, t129 + t168, t39 * t101 + t118, t46, t32, t61, t62, 0, t121 * t108 + (t125 - t20) * t112 + t153, t121 * t112 + (-t125 - t188) * t108 + t173, t2, t1, t21, t22, 0, t90 * t14 + (-t107 * t77 - t111 * t76) * t98 + t137 * t43 + (t107 * t127 - t111 * t128) * t100 + t142, t90 * t13 - (-t107 * t76 + t111 * t77) * t98 + t137 * t45 + (t107 * t128 + t111 * t127) * t100 + t122; 0, 0, 0, 0, 0, 0, 0, -t97 * t162, t158 * t97, t108 * t99, t170, qJDD(4), -g(1) * t112 + t108 * t123 + t91, -t108 * t159 + t112 * t123, t176, t15, t6, -t135, t98, -(-t107 * t25 - t171) * t100 + (-t100 * t157 + t111 * t98 - t165 * t43) * pkin(4) + t117, (-qJD(5) * t24 + t100 * t25 - t5) * t111 + (-qJD(5) * t100 * t111 - t107 * t98 - t165 * t45) * pkin(4) + t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t15, t6, -t135, t98, -t100 * t132 + t117, (-t5 + (-qJD(5) + t100) * t24) * t111 + t119;];
tau_reg = t3;
