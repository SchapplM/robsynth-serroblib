% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:58:26
% EndTime: 2019-05-05 14:58:32
% DurationCPUTime: 1.46s
% Computational Cost: add. (3985->230), mult. (7611->257), div. (0->0), fcn. (4291->6), ass. (0->160)
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t132 = cos(qJ(4));
t169 = qJD(1) * t132;
t104 = -qJD(4) * t131 + t128 * t169;
t129 = sin(qJ(4));
t165 = qJD(1) * qJD(4);
t159 = t129 * t165;
t163 = t132 * qJDD(1);
t109 = -t159 + t163;
t73 = -qJD(5) * t104 + qJDD(4) * t128 + t109 * t131;
t115 = qJD(1) * t129 + qJD(5);
t93 = t115 * t104;
t62 = t73 + t93;
t203 = qJ(6) * t62;
t158 = t132 * t165;
t164 = t129 * qJDD(1);
t108 = -t158 - t164;
t102 = qJDD(5) - t108;
t106 = qJD(4) * t128 + t131 * t169;
t83 = t106 * t104;
t196 = t102 - t83;
t202 = pkin(5) * t196;
t188 = pkin(1) + qJ(3);
t187 = qJ(2) - pkin(7);
t201 = t128 * t196;
t200 = t131 * t196;
t135 = qJD(1) ^ 2;
t127 = t135 * pkin(7);
t130 = sin(qJ(1));
t133 = cos(qJ(1));
t157 = g(1) * t130 - g(2) * t133;
t151 = -qJDD(2) + t157;
t140 = qJ(2) * t135 + t151;
t156 = t188 * qJDD(1);
t138 = t156 + t140;
t195 = 2 * qJD(3);
t48 = -pkin(4) * t108 - pkin(8) * t109 - t127 + (t195 + (pkin(4) * t132 + pkin(8) * t129) * qJD(4)) * qJD(1) + t138;
t134 = qJD(4) ^ 2;
t153 = pkin(4) * t129 - pkin(8) * t132;
t146 = t135 * t153;
t124 = qJDD(1) * qJ(2);
t152 = g(1) * t133 + g(2) * t130;
t148 = 0.2e1 * qJD(2) * qJD(1) - t152;
t141 = qJDD(3) + t148;
t87 = -t135 * t188 + t124 + t141;
t81 = -qJDD(1) * pkin(7) + t87;
t75 = t132 * g(3) - t129 * t81;
t54 = -t134 * pkin(4) + qJDD(4) * pkin(8) - t129 * t146 - t75;
t24 = t128 * t48 + t131 * t54;
t155 = -qJDD(4) * t131 + t109 * t128;
t72 = -qJD(5) * t106 - t155;
t88 = pkin(5) * t115 - qJ(6) * t106;
t147 = qJ(6) * t72 - 0.2e1 * qJD(6) * t104 - t115 * t88 + t24;
t100 = t104 ^ 2;
t101 = t106 ^ 2;
t114 = t115 ^ 2;
t78 = -t101 - t114;
t179 = t100 + t78;
t199 = pkin(5) * t179 - t147;
t197 = t73 - t93;
t58 = (qJD(5) - t115) * t106 + t155;
t76 = -t114 - t100;
t37 = t128 * t76 + t200;
t194 = pkin(4) * t37;
t67 = t102 + t83;
t182 = t128 * t67;
t42 = t131 * t78 - t182;
t193 = pkin(4) * t42;
t23 = t128 * t54 - t131 * t48;
t142 = t202 - t23 - t203;
t167 = qJD(6) * t106;
t97 = -0.2e1 * t167;
t10 = t142 + t97;
t192 = pkin(5) * t10;
t191 = pkin(5) * t62;
t190 = pkin(8) * t37;
t189 = pkin(8) * t42;
t31 = t128 * t62 - t131 * t58;
t65 = -t100 - t101;
t186 = -pkin(4) * t65 + pkin(8) * t31;
t38 = t131 * t76 - t201;
t57 = (qJD(5) + t115) * t106 + t155;
t185 = -pkin(4) * t57 + pkin(8) * t38;
t180 = t131 * t67;
t43 = -t128 * t78 - t180;
t184 = -pkin(4) * t197 + pkin(8) * t43;
t74 = g(3) * t129 + t132 * t81;
t53 = qJDD(4) * pkin(4) + pkin(8) * t134 - t132 * t146 + t74;
t183 = t128 * t53;
t181 = t131 * t53;
t178 = qJ(6) * t128;
t177 = qJ(6) * t131;
t176 = qJDD(1) * pkin(1);
t175 = t115 * t128;
t174 = t115 * t131;
t125 = t129 ^ 2;
t173 = t125 * t135;
t126 = t132 ^ 2;
t172 = t126 * t135;
t160 = t129 * t135 * t132;
t171 = t129 * (qJDD(4) + t160);
t170 = t125 + t126;
t162 = qJD(1) * t195;
t161 = t129 * t83;
t9 = t128 * t23 + t131 * t24;
t15 = t129 * t31 - t132 * t65;
t30 = -t128 * t58 - t131 * t62;
t154 = pkin(1) * t30 + t15 * t187;
t8 = t128 * t24 - t131 * t23;
t39 = -t129 * t75 + t132 * t74;
t149 = qJ(3) + t153;
t145 = pkin(1) + t149;
t19 = t129 * t38 - t132 * t57;
t144 = t187 * t19 + t188 * t37;
t21 = t129 * t43 - t132 * t197;
t143 = t187 * t21 + t188 * t42;
t139 = t142 + t202;
t137 = -t106 * t88 - qJDD(6) + t53;
t86 = t138 + t162;
t136 = -pkin(5) * t72 - t137;
t120 = 0.2e1 * t124;
t112 = t170 * t135;
t111 = t170 * qJDD(1);
t110 = -0.2e1 * t159 + t163;
t107 = 0.2e1 * t158 + t164;
t103 = t132 * (qJDD(4) - t160);
t98 = 0.2e1 * t167;
t94 = t140 + t176;
t90 = -t101 + t114;
t89 = t100 - t114;
t85 = -t171 + t132 * (-t134 - t172);
t84 = t129 * (-t134 - t173) + t103;
t80 = -t127 + t86;
t79 = t101 - t100;
t63 = (-t104 * t128 - t106 * t131) * t115;
t51 = t106 * t174 + t128 * t73;
t50 = t104 * t175 + t131 * t72;
t49 = t129 * t102 + t132 * (-t104 * t131 + t106 * t128) * t115;
t47 = t128 * t89 + t180;
t46 = t131 * t90 + t201;
t34 = t161 + t132 * (-t106 * t175 + t131 * t73);
t33 = -t161 + t132 * (t104 * t174 - t128 * t72);
t32 = -pkin(5) * t197 - qJ(6) * t67;
t29 = -t128 * t57 + t131 * t197;
t26 = t129 * t62 + t132 * (-t128 * t90 + t200);
t25 = -t129 * t58 + t132 * (t131 * t89 - t182);
t17 = qJ(6) * t100 - t136;
t16 = t129 * t79 + t132 * (-t128 * t197 - t131 * t57);
t13 = -qJ(6) * t179 + t136;
t12 = -pkin(5) * t100 + t147;
t11 = (t100 + t76) * qJ(6) + (-t57 + t72) * pkin(5) + t137;
t7 = -t142 + t98 + t203;
t6 = -qJ(6) * t58 + (-t100 - t65) * pkin(5) + t147;
t5 = t129 * t9 + t132 * t53;
t4 = pkin(5) * t17 + qJ(6) * t12;
t3 = -t10 * t128 + t12 * t131;
t2 = t10 * t131 + t12 * t128;
t1 = t129 * t3 + t132 * t17;
t14 = [0, 0, 0, 0, 0, qJDD(1), t157, t152, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t151 - 0.2e1 * t176, t120 + t148, pkin(1) * t94 + qJ(2) * (-pkin(1) * t135 + t124 + t148), qJDD(1), 0, 0, 0, 0, 0, 0, t120 + t141, t151 + 0.2e1 * t156 + t162, qJ(2) * t87 + t188 * t86, (t109 - t159) * t132, -t107 * t132 - t110 * t129, -t129 * (t134 - t172) + t103, (-t108 + t158) * t129, -t171 + t132 * (-t134 + t173), 0, t107 * t188 + t129 * t80 + t187 * t84, t110 * t188 + t132 * t80 + t187 * t85, -t111 * t187 - t112 * t188 - t39, t187 * t39 + t188 * t80, t34, t16, t26, t33, t25, t49, -t129 * (t23 - t194) + t132 * (-t183 - t190) + t144, -t129 * (t24 - t193) + t132 * (-t181 - t189) + t143, -t132 * t8 + t149 * t30 + t154, t145 * t8 + t187 * t5, t34, t16, t26, t33, t25, t49, -t129 * (-t139 + t98 - t194) + t132 * (-t11 * t128 - t177 * t196 - t190) + t144, -t129 * (-t193 - t199) + t132 * (-t128 * t32 + t13 * t131 - t189) + t143, qJ(3) * t30 - t129 * (-pkin(4) * t30 + t191) + t132 * (-pkin(8) * t30 - t128 * t6 + t131 * t7) + t154, t129 * t192 + t132 * (-t10 * t177 - t128 * t4) + t145 * t2 + t187 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t135, -t94, 0, 0, 0, 0, 0, 0, 0, -t135, -qJDD(1), -t86, 0, 0, 0, 0, 0, 0, -t107, -t110, t112, -t80, 0, 0, 0, 0, 0, 0, -t37, -t42, -t30, -t8, 0, 0, 0, 0, 0, 0, -t37, -t42, -t30, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t135, t87, 0, 0, 0, 0, 0, 0, t84, t85, -t111, t39, 0, 0, 0, 0, 0, 0, t19, t21, t15, t5, 0, 0, 0, 0, 0, 0, t19, t21, t15, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, (-t125 + t126) * t135, t163, -t160, -t164, qJDD(4), t74, t75, 0, 0, t51, t29, t46, t50, t47, t63, t181 + t185, -t183 + t184, t9 + t186, pkin(4) * t53 + pkin(8) * t9, t51, t29, t46, t50, t47, t63, t11 * t131 - t178 * t196 + t185, t128 * t13 + t131 * t32 + t184, t128 * t7 + t131 * t6 + t186, pkin(4) * t17 + pkin(8) * t3 - t10 * t178 + t131 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t79, t62, -t83, -t58, t102, -t23, -t24, 0, 0, t83, t79, t62, -t83, -t58, t102, t139 + t97, t199, -t191, t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t197, t65, -t17;];
tauJ_reg  = t14;
