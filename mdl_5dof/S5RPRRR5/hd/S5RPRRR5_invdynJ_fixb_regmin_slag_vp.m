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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:54:10
% EndTime: 2020-01-03 11:54:14
% DurationCPUTime: 0.99s
% Computational Cost: add. (1295->170), mult. (2151->235), div. (0->0), fcn. (1414->14), ass. (0->135)
t189 = pkin(7) + pkin(8);
t105 = sin(qJ(5));
t106 = sin(qJ(4));
t109 = cos(qJ(5));
t110 = cos(qJ(4));
t52 = t105 * t110 + t109 * t106;
t99 = qJD(1) + qJD(3);
t45 = t52 * t99;
t90 = qJ(1) + pkin(9) + qJ(3);
t83 = sin(t90);
t84 = cos(t90);
t135 = -g(2) * t84 - g(3) * t83;
t107 = sin(qJ(3));
t111 = cos(qJ(3));
t103 = sin(pkin(9));
t177 = pkin(1) * t103;
t104 = cos(pkin(9));
t85 = t104 * pkin(1) + pkin(2);
t71 = t85 * qJD(1);
t187 = qJD(3) * t71 + qJDD(1) * t177;
t152 = qJD(1) * t177;
t188 = -qJD(3) * t152 + t85 * qJDD(1);
t138 = -t187 * t107 + t188 * t111;
t97 = qJDD(1) + qJDD(3);
t180 = t97 * pkin(3);
t20 = -t138 - t180;
t126 = t135 - t20;
t186 = g(2) * t83 - g(3) * t84;
t159 = t109 * t110;
t163 = t105 * t106;
t51 = -t159 + t163;
t140 = -t107 * t177 + t111 * t85;
t40 = t107 * t71 + t111 * t152;
t145 = t189 * t99 + t40;
t25 = t110 * qJD(2) - t145 * t106;
t98 = qJD(4) + qJD(5);
t26 = t106 * qJD(2) + t145 * t110;
t185 = t188 * t107 + t187 * t111;
t179 = t99 * pkin(3);
t169 = t107 * t85 + t111 * t177;
t48 = pkin(7) + t169;
t178 = -pkin(8) - t48;
t176 = t110 * pkin(4);
t39 = -t107 * t152 + t111 * t71;
t175 = t39 * t98;
t174 = t40 * t99;
t42 = t169 * qJD(3);
t173 = t42 * t99;
t149 = t99 * t159;
t150 = t99 * t163;
t43 = -t149 + t150;
t172 = t45 * t43;
t102 = qJ(4) + qJ(5);
t92 = sin(t102);
t171 = t83 * t92;
t170 = t84 * t92;
t168 = t106 * t99;
t167 = t109 * t26;
t166 = t110 * t97;
t161 = t106 * t110;
t158 = qJDD(2) - g(1);
t100 = t106 ^ 2;
t157 = -t110 ^ 2 + t100;
t156 = qJD(5) * t105;
t154 = t106 * qJD(4);
t153 = t110 * qJD(4);
t151 = pkin(4) * t154;
t148 = t99 * t153;
t88 = -pkin(3) - t176;
t19 = t97 * pkin(7) + t185;
t146 = pkin(8) * t97 + t19;
t144 = qJD(4) * t189;
t143 = qJD(4) * t178;
t10 = (t99 * t154 - t166) * pkin(4) + t20;
t27 = t88 * t99 - t39;
t30 = t98 * t51;
t142 = g(2) * t170 + g(3) * t171 + t10 * t52 - t27 * t30;
t33 = -t39 - t179;
t139 = -t126 * t106 + t33 * t153;
t47 = -pkin(3) - t140;
t136 = -t40 + t151;
t134 = t51 * t97;
t108 = sin(qJ(1));
t112 = cos(qJ(1));
t132 = -g(2) * t112 - g(3) * t108;
t96 = qJDD(4) + qJDD(5);
t21 = -t30 * t98 + t52 * t96;
t24 = qJD(4) * pkin(4) + t25;
t131 = -t105 * t24 - t167;
t35 = t178 * t106;
t94 = t110 * pkin(8);
t36 = t110 * t48 + t94;
t130 = -t105 * t36 + t109 * t35;
t129 = t105 * t35 + t109 * t36;
t76 = t189 * t106;
t77 = t110 * pkin(7) + t94;
t128 = -t105 * t77 - t109 * t76;
t127 = -t105 * t76 + t109 * t77;
t113 = qJD(4) ^ 2;
t124 = pkin(7) * t113 - t174 - t180;
t123 = t113 * t48 + t47 * t97 + t173;
t122 = -t33 * t99 + t186 - t19;
t121 = t135 + t138;
t120 = -pkin(7) * qJDD(4) + (t39 - t179) * qJD(4);
t31 = t98 * t52;
t93 = cos(t102);
t119 = t10 * t51 + t135 * t93 + t27 * t31;
t41 = t140 * qJD(3);
t118 = -qJDD(4) * t48 + (t47 * t99 - t41) * qJD(4);
t13 = qJD(5) * t149 + t109 * t148 - t98 * t150 + t52 * t97;
t89 = t110 * qJDD(2);
t4 = qJDD(4) * pkin(4) - t26 * qJD(4) - t146 * t106 + t89;
t117 = t27 * t43 + t26 * t156 + g(1) * t92 + (-t26 * t98 - t4) * t105 + t186 * t93;
t116 = -t185 + t186;
t5 = t25 * qJD(4) + t106 * qJDD(2) + t146 * t110;
t115 = -g(1) * t93 + g(2) * t171 - g(3) * t170 + t131 * qJD(5) - t105 * t5 + t109 * t4 - t27 * t45;
t95 = t99 ^ 2;
t62 = qJDD(4) * t110 - t113 * t106;
t61 = qJDD(4) * t106 + t113 * t110;
t58 = t110 * t144;
t57 = t106 * t144;
t46 = t100 * t97 + 0.2e1 * t106 * t148;
t38 = t47 - t176;
t37 = t42 + t151;
t32 = -0.2e1 * t157 * t99 * qJD(4) + 0.2e1 * t97 * t161;
t28 = t33 * t154;
t22 = -t31 * t98 - t51 * t96;
t18 = -t106 * t41 + t110 * t143;
t17 = t106 * t143 + t110 * t41;
t15 = -t43 ^ 2 + t45 ^ 2;
t14 = t31 * t99 + t134;
t6 = t43 * t98 + t13;
t2 = t13 * t52 - t45 * t30;
t1 = -t13 * t51 - t52 * t14 + t30 * t43 - t45 * t31;
t3 = [qJDD(1), t132, g(2) * t108 - g(3) * t112, (t132 + (t103 ^ 2 + t104 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t97, t140 * t97 + t121 - t173, -t169 * t97 - t41 * t99 + t116, t46, t32, t61, t62, 0, t28 + t118 * t106 + (-t123 + t126) * t110, t123 * t106 + t118 * t110 + t139, t2, t1, t21, t22, 0, t37 * t43 + t38 * t14 + (-t129 * qJD(5) - t105 * t17 + t109 * t18) * t98 + t130 * t96 + t119, t37 * t45 + t38 * t13 - (t130 * qJD(5) + t105 * t18 + t109 * t17) * t98 - t129 * t96 + t142; 0, 0, 0, t158, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, 0, 0, 0, 0, 0, t22, -t21; 0, 0, 0, 0, t97, t121 + t174, t39 * t99 + t116, t46, t32, t61, t62, 0, t28 + t120 * t106 + (-t124 + t126) * t110, t124 * t106 + t120 * t110 + t139, t2, t1, t21, t22, 0, t88 * t14 + (-t127 * qJD(5) + t105 * t57 - t109 * t58) * t98 + t128 * t96 + t136 * t43 + t52 * t175 + t119, t88 * t13 - (t128 * qJD(5) - t105 * t58 - t109 * t57) * t98 - t127 * t96 + t136 * t45 - t51 * t175 + t142; 0, 0, 0, 0, 0, 0, 0, -t95 * t161, t157 * t95, t106 * t97, t166, qJDD(4), -g(1) * t110 + t122 * t106 + t89, -t158 * t106 + t122 * t110, t172, t15, t6, -t134, t96, -(-t105 * t25 - t167) * t98 + (t109 * t96 - t98 * t156 - t43 * t168) * pkin(4) + t115, (-qJD(5) * t24 + t25 * t98 - t5) * t109 + (-qJD(5) * t109 * t98 - t105 * t96 - t45 * t168) * pkin(4) + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t15, t6, -t134, t96, -t131 * t98 + t115, (-t5 + (-qJD(5) + t98) * t24) * t109 + t117;];
tau_reg = t3;
