% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:23
% EndTime: 2021-01-15 14:30:27
% DurationCPUTime: 0.95s
% Computational Cost: add. (1180->202), mult. (2779->259), div. (0->0), fcn. (1854->8), ass. (0->127)
t159 = cos(qJ(3));
t164 = pkin(5) + pkin(6);
t93 = cos(qJ(2));
t57 = t164 * t93;
t52 = qJD(1) * t57;
t90 = sin(qJ(3));
t36 = t90 * t52;
t145 = qJD(2) * pkin(2);
t91 = sin(qJ(2));
t56 = t164 * t91;
t50 = qJD(1) * t56;
t42 = -t50 + t145;
t126 = t159 * t42 - t36;
t45 = t159 * t91 + t90 * t93;
t34 = t45 * qJD(1);
t141 = t34 * qJ(4);
t9 = t126 - t141;
t86 = qJD(2) + qJD(3);
t83 = t93 * pkin(2);
t160 = pkin(1) + t83;
t148 = t159 * t57 - t56 * t90;
t123 = t159 * qJD(3);
t162 = pkin(2) * t86;
t84 = qJDD(2) + qJDD(3);
t169 = -pkin(2) * t84 * t90 - t123 * t162;
t150 = t90 * t91;
t113 = t86 * t150;
t130 = t159 * t93;
t117 = qJD(1) * t130;
t121 = qJDD(1) * t159;
t136 = t93 * qJDD(1);
t120 = -t117 * t86 - t91 * t121 - t90 * t136;
t11 = qJD(1) * t113 + t120;
t144 = t11 * qJ(4);
t80 = t84 * pkin(3);
t168 = t144 + t80;
t135 = qJD(1) * qJD(2);
t127 = t93 * t135;
t137 = t91 * qJDD(1);
t22 = qJDD(2) * pkin(2) + t164 * (-t127 - t137);
t128 = t91 * t135;
t24 = t164 * (-t128 + t136);
t167 = t159 * t22 - t90 * t24;
t89 = qJ(2) + qJ(3);
t81 = sin(t89);
t94 = cos(qJ(1));
t154 = t81 * t94;
t92 = sin(qJ(1));
t155 = t81 * t92;
t82 = cos(t89);
t161 = g(3) * t82;
t166 = g(1) * t154 + g(2) * t155 - t161;
t165 = t34 ^ 2;
t8 = pkin(3) * t86 + t9;
t163 = t8 - t9;
t140 = qJD(1) * t91;
t132 = t90 * t140;
t32 = -t117 + t132;
t158 = t32 * t86;
t157 = t34 * t32;
t153 = t82 * t92;
t152 = t82 * t94;
t149 = -t159 * t50 - t36;
t147 = pkin(3) * t82 + t83;
t87 = t91 ^ 2;
t146 = -t93 ^ 2 + t87;
t114 = -t121 * t93 + t137 * t90;
t21 = t86 * t45;
t12 = qJD(1) * t21 + t114;
t143 = t12 * qJ(4);
t142 = t32 * qJ(4);
t139 = qJD(3) * t90;
t122 = t32 * pkin(3) + qJD(4);
t55 = t160 * qJD(1);
t23 = t122 - t55;
t138 = qJD(4) + t23;
t134 = t91 * t145;
t133 = pkin(2) * t140;
t40 = t159 * t52;
t129 = qJD(2) * t164;
t125 = t50 * t90 - t40;
t124 = -t159 * t56 - t57 * t90;
t119 = -g(1) * t155 + g(2) * t154;
t118 = g(1) * t153 - g(2) * t152;
t116 = g(1) * t94 + g(2) * t92;
t115 = g(1) * t92 - g(2) * t94;
t112 = -t42 * t90 - t40;
t51 = t91 * t129;
t53 = t93 * t129;
t111 = -t123 * t56 - t139 * t57 - t159 * t51 - t53 * t90;
t29 = pkin(2) * t128 - qJDD(1) * t160;
t110 = -0.2e1 * pkin(1) * t135 - pkin(5) * qJDD(2);
t95 = qJD(2) ^ 2;
t107 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t95 + t115;
t96 = qJD(1) ^ 2;
t106 = pkin(1) * t96 - pkin(5) * qJDD(1) + t116;
t5 = t12 * pkin(3) + qJDD(4) + t29;
t105 = -t132 * t86 - t120;
t104 = qJD(3) * t112 + t167;
t103 = -qJD(3) * t148 - t159 * t53 + t90 * t51;
t102 = t123 * t42 - t139 * t52 + t159 * t24 + t90 * t22;
t101 = g(1) * t152 + g(2) * t153 + g(3) * t81 - t102;
t100 = t104 + t166;
t99 = -t55 * t32 + t101;
t98 = t55 * t34 + t100;
t97 = t138 * t32 + t101 + t143;
t85 = -qJ(4) - t164;
t78 = pkin(2) * t159 + pkin(3);
t49 = pkin(1) + t147;
t44 = -t130 + t150;
t31 = t32 ^ 2;
t26 = pkin(3) * t44 - t160;
t25 = pkin(3) * t34 + t133;
t20 = -qJD(2) * t130 - t123 * t93 + t113;
t18 = pkin(3) * t21 + t134;
t17 = -qJ(4) * t44 + t148;
t16 = -qJ(4) * t45 + t124;
t15 = -t31 + t165;
t14 = -t141 + t149;
t13 = t125 + t142;
t10 = -t112 - t142;
t6 = t105 + t158;
t4 = t20 * qJ(4) - t45 * qJD(4) + t103;
t3 = -qJ(4) * t21 - qJD(4) * t44 + t111;
t2 = -qJD(4) * t32 + t102 - t143;
t1 = -t34 * qJD(4) + t104 + t168;
t7 = [qJDD(1), t115, t116, qJDD(1) * t87 + 0.2e1 * t127 * t91, -0.2e1 * t135 * t146 + 0.2e1 * t136 * t91, qJDD(2) * t91 + t93 * t95, qJDD(2) * t93 - t91 * t95, 0, t107 * t93 + t110 * t91, -t107 * t91 + t110 * t93, -t11 * t45 - t20 * t34, t11 * t44 - t12 * t45 + t20 * t32 - t21 * t34, -t20 * t86 + t45 * t84, -t21 * t86 - t44 * t84, 0, t103 * t86 - t12 * t160 + t124 * t84 + t134 * t32 - t55 * t21 + t29 * t44 + t118, t11 * t160 - t111 * t86 + t134 * t34 - t148 * t84 + t20 * t55 + t29 * t45 + t119, t12 * t26 + t16 * t84 + t18 * t32 + t21 * t23 + t4 * t86 + t44 * t5 + t118, -t11 * t26 - t17 * t84 + t18 * t34 - t20 * t23 - t3 * t86 + t45 * t5 + t119, -t1 * t45 - t10 * t21 + t11 * t16 - t12 * t17 - t2 * t44 + t20 * t8 - t3 * t32 - t34 * t4 - t116, t2 * t17 + t10 * t3 + t1 * t16 + t8 * t4 + t5 * t26 + t23 * t18 - g(1) * (-t49 * t92 - t85 * t94) - g(2) * (t49 * t94 - t85 * t92); 0, 0, 0, -t91 * t96 * t93, t146 * t96, t137, t136, qJDD(2), -g(3) * t93 + t106 * t91, g(3) * t91 + t106 * t93, t157, t15, t6, -t114, t84, -t125 * t86 + (-t139 * t86 - t140 * t32 + t159 * t84) * pkin(2) + t98, -t133 * t34 + t149 * t86 + t169 + t99, -t13 * t86 - t25 * t32 + t78 * t84 - t138 * t34 + (-t40 + (-t42 - t162) * t90) * qJD(3) + t166 + t167 + t168, t14 * t86 - t25 * t34 + t169 + t97, t78 * t11 + (t10 + t13) * t34 + (t14 - t8) * t32 + (-t12 * t90 + (-t159 * t32 + t34 * t90) * qJD(3)) * pkin(2), t1 * t78 - t10 * t14 - t8 * t13 - t23 * t25 - g(3) * t147 - t116 * (-pkin(2) * t91 - pkin(3) * t81) + (t2 * t90 + (t10 * t159 - t8 * t90) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t15, t6, -t114, t84, -t112 * t86 + t98, t126 * t86 + t99, t144 + t10 * t86 + 0.2e1 * t80 + (-t122 - t23) * t34 + t100, -pkin(3) * t165 + t9 * t86 + t97, t11 * pkin(3) - t163 * t32, t163 * t10 + (t116 * t81 - t23 * t34 + t1 - t161) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t86 + t12, t105 - t158, -t31 - t165, t10 * t32 + t8 * t34 - t115 + t5;];
tau_reg = t7;
