% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:47
% EndTime: 2019-12-05 17:47:51
% DurationCPUTime: 1.24s
% Computational Cost: add. (1426->212), mult. (2831->279), div. (0->0), fcn. (1982->10), ass. (0->131)
t130 = qJD(1) ^ 2;
t124 = sin(qJ(1));
t127 = cos(qJ(1));
t167 = g(1) * t124 - g(2) * t127;
t135 = -t130 * qJ(2) - t167;
t113 = qJD(3) + qJD(5);
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t119 = sin(pkin(8));
t120 = cos(pkin(8));
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t141 = t119 * t126 + t120 * t123;
t70 = t141 * qJD(1);
t60 = t125 * t70;
t163 = qJD(1) * t126;
t164 = qJD(1) * t123;
t73 = -t119 * t164 + t120 * t163;
t32 = t122 * t73 + t60;
t175 = t32 * t113;
t160 = qJD(5) * t122;
t140 = t119 * t123 - t120 * t126;
t158 = qJD(1) * qJD(3);
t39 = -qJDD(1) * t141 + t140 * t158;
t154 = t126 * qJDD(1);
t155 = t123 * qJDD(1);
t40 = t119 * t155 - t120 * t154 + t141 * t158;
t4 = -qJD(5) * t60 + t122 * t39 - t125 * t40 - t160 * t73;
t194 = t4 + t175;
t143 = -t122 * t70 + t125 * t73;
t193 = t143 * t32;
t176 = t143 * t113;
t5 = qJD(5) * t143 - t122 * t40 - t125 * t39;
t192 = -t5 + t176;
t191 = t143 ^ 2 - t32 ^ 2;
t183 = t70 * pkin(7);
t128 = -pkin(1) - pkin(6);
t88 = qJD(1) * t128 + qJD(2);
t66 = -qJ(4) * t164 + t123 * t88;
t178 = t120 * t66;
t67 = -qJ(4) * t163 + t126 * t88;
t57 = qJD(3) * pkin(3) + t67;
t23 = t119 * t57 + t178;
t12 = t23 - t183;
t82 = pkin(3) * t164 + qJD(1) * qJ(2) + qJD(4);
t45 = pkin(4) * t70 + t82;
t105 = qJ(3) + pkin(8) + qJ(5);
t97 = sin(t105);
t98 = cos(t105);
t190 = g(3) * t98 + t12 * t160 + t167 * t97 + t45 * t32;
t42 = -t122 * t141 - t125 * t140;
t112 = qJDD(3) + qJDD(5);
t41 = -t122 * t140 + t125 * t141;
t161 = qJD(3) * t126;
t162 = qJD(3) * t123;
t71 = t119 * t162 - t120 * t161;
t72 = t141 * qJD(3);
t8 = -qJD(5) * t41 + t122 * t71 - t125 * t72;
t188 = t112 * t42 + t113 * t8;
t187 = qJD(5) - t113;
t114 = qJDD(1) * qJ(2);
t149 = g(1) * t127 + g(2) * t124;
t115 = qJD(1) * qJD(2);
t153 = 0.2e1 * t115;
t186 = 0.2e1 * t114 + t153 - t149;
t157 = qJD(1) * qJD(4);
t87 = qJDD(1) * t128 + qJDD(2);
t80 = t126 * t87;
t21 = -t126 * t157 - t88 * t162 + qJDD(3) * pkin(3) + t80 + (t123 * t158 - t154) * qJ(4);
t26 = (-qJ(4) * qJD(1) + t88) * t161 + (-qJ(4) * qJDD(1) - t157 + t87) * t123;
t6 = -t119 * t26 + t120 * t21;
t2 = qJDD(3) * pkin(4) + pkin(7) * t40 + t6;
t7 = t119 * t21 + t120 * t26;
t3 = pkin(7) * t39 + t7;
t185 = g(3) * t97 - t122 * t3 + t125 * t2 - t45 * t143 - t167 * t98;
t182 = t73 * pkin(7);
t181 = pkin(3) * t119;
t180 = g(3) * t123;
t108 = t123 * pkin(3);
t169 = qJ(4) - t128;
t64 = -t126 * qJD(4) + t162 * t169;
t84 = t169 * t126;
t65 = -qJD(3) * t84 - t123 * qJD(4);
t25 = t119 * t64 + t120 * t65;
t53 = t119 * t66;
t30 = t120 * t67 - t53;
t83 = t169 * t123;
t44 = -t119 * t84 - t120 * t83;
t173 = pkin(1) * qJDD(1);
t171 = t82 * qJD(1);
t170 = qJ(2) + t108;
t168 = pkin(1) * t127 + qJ(2) * t124;
t118 = t126 ^ 2;
t166 = t123 ^ 2 - t118;
t129 = qJD(3) ^ 2;
t165 = -t129 - t130;
t159 = pkin(3) * t161 + qJD(2);
t156 = qJDD(3) * t123;
t152 = t126 * t158;
t22 = t120 * t57 - t53;
t24 = -t119 * t65 + t120 * t64;
t29 = -t119 * t67 - t178;
t43 = t119 * t83 - t120 * t84;
t150 = qJDD(2) - t173;
t9 = qJD(5) * t42 - t122 * t72 - t125 * t71;
t147 = -t112 * t41 - t113 * t9;
t11 = qJD(3) * pkin(4) - t182 + t22;
t146 = -t122 * t11 - t125 * t12;
t27 = pkin(7) * t140 + t43;
t28 = -pkin(7) * t141 + t44;
t145 = -t122 * t28 + t125 * t27;
t144 = t122 * t27 + t125 * t28;
t142 = qJDD(4) + t114 + t115 + (t152 + t155) * pkin(3);
t99 = pkin(3) * t120 + pkin(4);
t139 = t122 * t99 + t125 * t181;
t138 = -t122 * t181 + t125 * t99;
t136 = 0.2e1 * qJ(2) * t158 + qJDD(3) * t128;
t132 = -t140 * t6 + t141 * t7 - t22 * t72 - t23 * t71 - t167;
t131 = -t128 * t129 + t186;
t121 = -qJ(4) - pkin(6);
t107 = t127 * qJ(2);
t104 = qJDD(3) * t126;
t58 = pkin(4) * t141 + t170;
t49 = pkin(3) * t163 + pkin(4) * t73;
t46 = -pkin(4) * t71 + t159;
t17 = -pkin(4) * t39 + t142;
t16 = t30 - t182;
t15 = t29 + t183;
t14 = pkin(7) * t71 + t25;
t13 = pkin(7) * t72 + t24;
t1 = [qJDD(1), t167, t149, qJDD(2) - t167 - 0.2e1 * t173, t186, -t150 * pkin(1) - g(1) * (-pkin(1) * t124 + t107) - g(2) * t168 + (t153 + t114) * qJ(2), qJDD(1) * t118 - 0.2e1 * t123 * t152, -0.2e1 * t123 * t154 + 0.2e1 * t158 * t166, -t123 * t129 + t104, -t126 * t129 - t156, 0, t123 * t131 + t126 * t136, -t123 * t136 + t126 * t131, -t24 * t73 - t25 * t70 + t39 * t44 + t40 * t43 - t132, t7 * t44 + t23 * t25 + t6 * t43 + t22 * t24 + t142 * t170 + t82 * t159 - g(1) * (t127 * t108 + t107 + (-pkin(1) + t121) * t124) - g(2) * (t108 * t124 - t121 * t127 + t168), t143 * t8 + t4 * t42, -t143 * t9 - t32 * t8 - t4 * t41 - t42 * t5, t188, t147, 0, t46 * t32 + t58 * t5 + t17 * t41 + t45 * t9 + (-qJD(5) * t144 - t122 * t14 + t125 * t13) * t113 + t145 * t112 - t149 * t97, t46 * t143 + t58 * t4 + t17 * t42 + t45 * t8 - (qJD(5) * t145 + t122 * t13 + t125 * t14) * t113 - t144 * t112 - t149 * t98; 0, 0, 0, qJDD(1), -t130, t150 + t135, 0, 0, 0, 0, 0, t123 * t165 + t104, t126 * t165 - t156, -t140 * t40 + t141 * t39 + t70 * t71 + t72 * t73, t132 - t171, 0, 0, 0, 0, 0, -qJD(1) * t32 + t188, -qJD(1) * t143 + t147; 0, 0, 0, 0, 0, 0, t126 * t130 * t123, -t166 * t130, t154, -t155, qJDD(3), t126 * t135 + t180 + t80, g(3) * t126 + (-t135 - t87) * t123, (t23 + t29) * t73 - (t22 - t30) * t70 + (t119 * t39 + t120 * t40) * pkin(3), -t22 * t29 - t23 * t30 + (t180 + t119 * t7 + t120 * t6 + (-t167 - t171) * t126) * pkin(3), t193, t191, t194, t192, t112, t138 * t112 - t49 * t32 - (-t122 * t16 + t125 * t15) * t113 + (-t113 * t139 + t146) * qJD(5) + t185, -t139 * t112 - t125 * t3 - t122 * t2 - t49 * t143 + (t122 * t15 + t125 * t16) * t113 + (-t125 * t11 - t113 * t138) * qJD(5) + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 ^ 2 - t73 ^ 2, t22 * t73 + t23 * t70 + t142 - t149, 0, 0, 0, 0, 0, t5 + t176, t4 - t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t191, t194, t192, t112, t146 * t187 + t185, (-t113 * t12 - t2) * t122 + (-t11 * t187 - t3) * t125 + t190;];
tau_reg = t1;
