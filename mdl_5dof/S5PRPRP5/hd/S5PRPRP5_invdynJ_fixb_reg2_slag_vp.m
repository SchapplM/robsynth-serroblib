% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:48
% EndTime: 2019-12-05 15:38:51
% DurationCPUTime: 1.54s
% Computational Cost: add. (1536->242), mult. (3478->297), div. (0->0), fcn. (2585->10), ass. (0->141)
t95 = cos(pkin(7));
t184 = g(1) * t95;
t93 = sin(pkin(7));
t134 = g(2) * t93 + t184;
t99 = cos(qJ(2));
t88 = g(3) * t99;
t98 = sin(qJ(2));
t114 = t134 * t98 - t88;
t180 = cos(qJ(4));
t94 = cos(pkin(8));
t148 = t180 * t94;
t92 = sin(pkin(8));
t97 = sin(qJ(4));
t173 = t97 * t92;
t120 = t148 - t173;
t117 = t99 * t120;
t172 = pkin(6) + qJ(3);
t69 = t172 * t92;
t70 = t172 * t94;
t121 = -t180 * t69 - t97 * t70;
t171 = -qJD(1) * t117 + t120 * qJD(3) + t121 * qJD(4);
t35 = t180 * t70 - t97 * t69;
t91 = pkin(8) + qJ(4);
t86 = sin(t91);
t196 = t171 * qJD(4) + t35 * qJDD(4) + t114 * t86;
t161 = t99 * qJD(1);
t137 = qJD(3) - t161;
t194 = t134 * t99;
t65 = t180 * t92 + t97 * t94;
t190 = t65 * qJD(2);
t185 = t190 ^ 2;
t135 = qJD(2) * t148;
t150 = qJD(2) * t173;
t56 = -t135 + t150;
t52 = t56 ^ 2;
t192 = -t52 - t185;
t191 = -t52 + t185;
t61 = t65 * qJD(4);
t144 = qJD(4) * t180;
t164 = qJD(4) * t97;
t149 = t92 * t164;
t189 = t94 * t144 - t149;
t142 = qJDD(2) * t180;
t156 = t94 * qJDD(2);
t151 = qJD(4) * t135 + t92 * t142 + t97 * t156;
t25 = qJD(2) * t149 - t151;
t157 = t92 * qJDD(2);
t133 = -t94 * t142 + t97 * t157;
t26 = qJD(2) * t61 + t133;
t188 = t26 * pkin(4) + t25 * qJ(5);
t154 = qJD(1) * qJD(2);
t163 = qJDD(2) * pkin(2);
t155 = t99 * qJDD(1);
t128 = t98 * t154 + qJDD(3) - t155;
t55 = t128 - t163;
t187 = (t134 + t154) * t98 + t163 - t55 - t88;
t182 = g(3) * t98;
t179 = t190 * t56;
t87 = cos(t91);
t178 = t87 * t98;
t177 = t93 * t87;
t176 = t93 * t99;
t175 = t95 * t99;
t167 = qJD(1) * t98;
t75 = qJD(2) * qJ(3) + t167;
t143 = pkin(6) * qJD(2) + t75;
t48 = t143 * t94;
t174 = t97 * t48;
t170 = t35 * qJD(4) + t137 * t65;
t89 = t92 ^ 2;
t90 = t94 ^ 2;
t169 = t89 + t90;
t166 = qJD(2) * t98;
t47 = t143 * t92;
t20 = t180 * t48 - t97 * t47;
t165 = qJD(4) * t20;
t162 = qJDD(4) * pkin(4);
t19 = -t180 * t47 - t174;
t160 = qJD(5) - t19;
t159 = qJDD(1) - g(3);
t153 = qJDD(4) * qJ(5);
t49 = qJDD(2) * qJ(3) + t98 * qJDD(1) + (qJD(3) + t161) * qJD(2);
t140 = pkin(6) * qJDD(2) + t49;
t32 = t140 * t92;
t33 = t140 * t94;
t152 = -t47 * t144 + t180 * t33 - t97 * t32;
t84 = pkin(3) * t94 + pkin(2);
t147 = t169 * t49;
t146 = t169 * t99;
t145 = t19 + t174;
t141 = t48 * t144 - t47 * t164 + t180 * t32 + t97 * t33;
t13 = pkin(4) * t61 - qJ(5) * t189 - qJD(5) * t65;
t139 = -t13 + t167;
t138 = t169 * qJDD(2);
t132 = pkin(4) * t87 + qJ(5) * t86;
t131 = -t120 * t26 + t56 * t61;
t129 = t84 * qJDD(2);
t100 = qJD(2) ^ 2;
t127 = qJDD(2) * t99 - t100 * t98;
t126 = qJD(4) * t61 - qJDD(4) * t120;
t21 = qJD(2) * t117 - t98 * t61;
t22 = t189 * t98 + t190 * t99;
t50 = t65 * t98;
t51 = t120 * t98;
t119 = t190 * t22 - t21 * t56 - t25 * t50 - t51 * t26;
t115 = -t182 - t194;
t113 = -qJD(4) * t22 - qJDD(4) * t50 + t56 * t166 - t26 * t99;
t43 = t86 * t176 + t87 * t95;
t45 = t86 * t175 - t177;
t112 = g(1) * t45 + g(2) * t43 + t86 * t182 - t141;
t111 = t137 * t169;
t62 = -t84 * qJD(2) + t137;
t110 = -t120 * t25 - t189 * t56 - t190 * t61 - t26 * t65;
t44 = t87 * t176 - t95 * t86;
t46 = t87 * t175 + t86 * t93;
t109 = g(1) * t46 + g(2) * t44 + g(3) * t178 - t152;
t39 = -t129 + t128;
t108 = qJD(4) * t21 + qJDD(4) * t51 - t166 * t190 - t25 * t99;
t107 = g(2) * t98 * t177 - t170 * qJD(4) + qJDD(4) * t121 + t178 * t184 - t87 * t88;
t106 = t128 - t114;
t105 = t147 + t115;
t14 = t56 * pkin(4) - qJ(5) * t190 + t62;
t104 = t14 * t190 + qJDD(5) - t112;
t103 = 0.2e1 * t190 * qJD(4) + t133;
t102 = -t129 + t106;
t101 = t25 * t121 + t170 * t190 - t171 * t56 - t35 * t26 + t115;
t72 = t99 * t84;
t71 = -qJD(2) * pkin(2) + t137;
t27 = qJD(4) * t189 + qJDD(4) * t65;
t24 = -pkin(4) * t120 - qJ(5) * t65 - t84;
t23 = pkin(4) * t190 + qJ(5) * t56;
t16 = qJD(4) * qJ(5) + t20;
t15 = -qJD(4) * pkin(4) + t160;
t12 = (t56 - t150) * qJD(4) + t151;
t11 = (t56 + t150) * qJD(4) - t151;
t6 = t189 * t190 - t25 * t65;
t4 = -t48 * t164 + t152;
t3 = -qJD(5) * t190 + t188 + t39;
t2 = qJDD(5) + t141 - t162;
t1 = t153 + (qJD(5) - t174) * qJD(4) + t152;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t159, 0, 0, 0, 0, 0, 0, t127, -qJDD(2) * t98 - t100 * t99, 0, -g(3) + (t98 ^ 2 + t99 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t127 * t94, -t127 * t92, t100 * t146 + t98 * t138, -t55 * t99 - g(3) + t98 * t147 + (t75 * t146 + t71 * t98) * qJD(2), 0, 0, 0, 0, 0, 0, t113, -t108, t119, t141 * t50 + t62 * t166 - t19 * t22 + t20 * t21 - t39 * t99 + t4 * t51 - g(3), 0, 0, 0, 0, 0, 0, t113, t119, t108, t1 * t51 + t14 * t166 + t15 * t22 + t16 * t21 + t2 * t50 - t3 * t99 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t155 + t114, -t159 * t98 + t194, 0, 0, t89 * qJDD(2), 0.2e1 * t92 * t156, 0, t90 * qJDD(2), 0, 0, t187 * t94, -t187 * t92, qJ(3) * t138 + t111 * qJD(2) + t105, -t71 * t167 + (-t55 + t114) * pkin(2) + t105 * qJ(3) + t111 * t75, t6, t110, t27, t131, -t126, 0, -t120 * t39 - t56 * t167 - t84 * t26 + t62 * t61 + t107, -t167 * t190 + t189 * t62 + t84 * t25 + t39 * t65 - t196, t120 * t4 + t141 * t65 - t189 * t19 - t20 * t61 + t101, t4 * t35 - t141 * t121 - t39 * t84 - t62 * t167 - g(3) * (t172 * t98 + t72) + t171 * t20 - t170 * t19 + t134 * (-t172 * t99 + t84 * t98), t6, t27, -t110, 0, t126, t131, -t120 * t3 - t139 * t56 + t14 * t61 + t24 * t26 + t107, t1 * t120 + t15 * t189 - t16 * t61 + t2 * t65 + t101, t139 * t190 - t14 * t189 + t24 * t25 - t3 * t65 + t196, -g(3) * t72 + t1 * t35 + t14 * t13 - t2 * t121 + t3 * t24 + t171 * t16 + t170 * t15 + (-g(3) * t132 - t134 * t172) * t99 + (-g(3) * t172 - t14 * qJD(1) + t134 * (t132 + t84)) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, t157, -t169 * t100, -t169 * t75 * qJD(2) + t106 - t163, 0, 0, 0, 0, 0, 0, t103, -t11, t192, t19 * t190 + t20 * t56 + t102, 0, 0, 0, 0, 0, 0, t103, t192, t11, t16 * t56 + (-qJD(5) - t15) * t190 + t102 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t191, t12, -t179, -t133, qJDD(4), -t190 * t62 + t112 + t165, t145 * qJD(4) + t62 * t56 + t109, 0, 0, t179, t12, -t191, qJDD(4), t133, -t179, -t23 * t56 - t104 + 0.2e1 * t162 + t165, pkin(4) * t25 - t26 * qJ(5) + (t16 - t20) * t190 + (t15 - t160) * t56, 0.2e1 * t153 - t14 * t56 + t23 * t190 + (0.2e1 * qJD(5) - t145) * qJD(4) - t109, t1 * qJ(5) - t2 * pkin(4) - t14 * t23 - t15 * t20 - g(1) * (-pkin(4) * t45 + qJ(5) * t46) - g(2) * (-pkin(4) * t43 + qJ(5) * t44) - (-pkin(4) * t86 + qJ(5) * t87) * t182 + t160 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t179, t12, -qJD(4) ^ 2 - t185, -qJD(4) * t16 + t104 - t162;];
tau_reg = t5;
