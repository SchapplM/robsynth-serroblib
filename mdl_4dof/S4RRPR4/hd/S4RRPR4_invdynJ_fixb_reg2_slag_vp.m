% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:36
% EndTime: 2019-12-31 17:02:38
% DurationCPUTime: 0.80s
% Computational Cost: add. (1386->173), mult. (2126->223), div. (0->0), fcn. (1367->12), ass. (0->123)
t103 = sin(qJ(2));
t162 = t103 * pkin(1);
t139 = qJD(1) * t162;
t106 = cos(qJ(2));
t160 = t106 * pkin(1);
t153 = -qJD(2) * t139 + qJDD(1) * t160;
t135 = qJDD(3) - t153;
t93 = qJDD(1) + qJDD(2);
t163 = t93 * pkin(2);
t42 = t135 - t163;
t98 = qJ(1) + qJ(2);
t89 = cos(t98);
t81 = g(2) * t89;
t167 = t42 + t81;
t88 = sin(t98);
t82 = g(1) * t88;
t165 = t81 - t82;
t99 = sin(pkin(7));
t94 = t99 ^ 2;
t100 = cos(pkin(7));
t95 = t100 ^ 2;
t152 = t94 + t95;
t105 = cos(qJ(4));
t147 = t105 * t100;
t102 = sin(qJ(4));
t150 = t102 * t99;
t114 = t147 - t150;
t52 = t102 * t100 + t105 * t99;
t146 = qJD(1) * t106;
t118 = -pkin(1) * t146 + qJD(3);
t97 = qJD(1) + qJD(2);
t45 = t52 * t97;
t166 = g(1) * t89 + g(2) * t88;
t164 = t45 ^ 2;
t104 = sin(qJ(1));
t161 = t104 * pkin(1);
t141 = t97 * t150;
t43 = -t97 * t147 + t141;
t159 = t45 * t43;
t101 = -pkin(6) - qJ(3);
t61 = t101 * t99;
t90 = t100 * pkin(6);
t62 = t100 * qJ(3) + t90;
t32 = -t102 * t62 + t105 * t61;
t158 = qJD(4) * t32 + t118 * t114;
t33 = t102 * t61 + t105 * t62;
t157 = -qJD(4) * t33 - t118 * t52;
t156 = t167 * t99;
t155 = t89 * pkin(2) + t88 * qJ(3);
t151 = t100 * t93;
t145 = qJD(2) * t103;
t144 = qJD(2) * t106;
t143 = qJDD(1) * t103;
t128 = qJD(4) * t147;
t142 = t97 * t128 + t52 * t93;
t140 = pkin(1) * t145;
t137 = t97 * t145;
t136 = qJD(4) * t150;
t78 = t100 * pkin(3) + pkin(2);
t134 = -t88 * pkin(2) + t89 * qJ(3);
t34 = t93 * qJ(3) + t97 * qJD(3) + (qJD(1) * t144 + t143) * pkin(1);
t133 = pkin(6) * t93 + t34;
t54 = t97 * qJ(3) + t139;
t132 = pkin(6) * t97 + t54;
t69 = pkin(1) * t144 + qJD(3);
t130 = t152 * t69;
t77 = qJ(3) + t162;
t129 = t152 * t77;
t127 = qJ(3) * t152;
t126 = -t88 * t101 + t89 * t78;
t21 = t133 * t99;
t22 = t133 * t100;
t125 = -t102 * t22 - t105 * t21;
t124 = t152 * t34 - t166;
t123 = -t153 + t165;
t122 = t97 * t139;
t121 = t114 * t93;
t35 = t132 * t99;
t36 = t132 * t100;
t12 = -t102 * t36 - t105 * t35;
t13 = -t102 * t35 + t105 * t36;
t115 = -t102 * t21 + t105 * t22;
t4 = t12 * qJD(4) + t115;
t47 = -t128 + t136;
t48 = t52 * qJD(4);
t5 = -t13 * qJD(4) + t125;
t119 = t114 * t4 + t12 * t47 - t13 * t48 - t5 * t52 - t166;
t107 = cos(qJ(1));
t117 = g(1) * t104 - g(2) * t107;
t116 = -t89 * t101 - t88 * t78;
t49 = (-pkin(6) - t77) * t99;
t50 = t100 * t77 + t90;
t25 = -t102 * t50 + t105 * t49;
t26 = t102 * t49 + t105 * t50;
t29 = -t78 * t93 + t135;
t41 = -t78 * t97 + t118;
t96 = pkin(7) + qJ(4);
t86 = sin(t96);
t113 = t165 * t86 + t29 * t52 - t41 * t47;
t87 = cos(t96);
t112 = -t29 * t114 - t165 * t87 + t41 * t48;
t111 = -t122 - t163;
t84 = -pkin(2) - t160;
t110 = pkin(1) * t137 + t84 * t93;
t109 = t118 * t152;
t91 = t107 * pkin(1);
t72 = t95 * t93;
t71 = t94 * t93;
t68 = t100 * t82;
t60 = -t78 - t160;
t56 = 0.2e1 * t99 * t151;
t53 = -t97 * pkin(2) + t118;
t40 = t43 ^ 2;
t28 = -t48 * qJD(4) + qJDD(4) * t114;
t27 = -t47 * qJD(4) + t52 * qJDD(4);
t17 = t48 * t97 - t121;
t16 = t97 * t136 - t142;
t11 = -qJD(4) * t26 - t52 * t69;
t10 = qJD(4) * t25 + t114 * t69;
t7 = -t114 * t17 + t43 * t48;
t6 = -t16 * t52 - t45 * t47;
t1 = -t114 * t16 - t52 * t17 + t47 * t43 - t45 * t48;
t2 = [0, 0, 0, 0, 0, qJDD(1), t117, g(1) * t107 + g(2) * t104, 0, 0, 0, 0, 0, 0, 0, t93, (t106 * t93 - t137) * pkin(1) - t123, ((-qJDD(1) - t93) * t103 + (-qJD(1) - t97) * t144) * pkin(1) + t166, 0, (t117 + (t103 ^ 2 + t106 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t71, t56, 0, t72, 0, 0, t68 + (-t110 - t167) * t100, (t110 - t82) * t99 + t156, t93 * t129 + t97 * t130 + t124, t42 * t84 + t53 * t140 - g(1) * (t134 - t161) - g(2) * (t91 + t155) + t54 * t130 + t34 * t129, t6, t1, t27, t7, t28, 0, t11 * qJD(4) + t25 * qJDD(4) + t43 * t140 + t60 * t17 + t112, -t10 * qJD(4) - t26 * qJDD(4) + t45 * t140 - t60 * t16 + t113, -t10 * t43 - t11 * t45 + t25 * t16 - t26 * t17 + t119, t4 * t26 + t13 * t10 + t5 * t25 + t12 * t11 + t29 * t60 + t41 * t140 - g(1) * (t116 - t161) - g(2) * (t126 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t122 - t123, (-t143 + (-qJD(2) + t97) * t146) * pkin(1) + t166, 0, 0, t71, t56, 0, t72, 0, 0, t68 + (-t111 - t167) * t100, (t111 - t82) * t99 + t156, t109 * t97 + t93 * t127 + t124, -t42 * pkin(2) - g(1) * t134 - g(2) * t155 + t109 * t54 + t34 * t127 - t53 * t139, t6, t1, t27, t7, t28, 0, t157 * qJD(4) + t32 * qJDD(4) - t43 * t139 - t78 * t17 + t112, -t158 * qJD(4) - t33 * qJDD(4) - t45 * t139 + t78 * t16 + t113, -t157 * t45 - t158 * t43 + t32 * t16 - t33 * t17 + t119, -g(1) * t116 - g(2) * t126 + t157 * t12 + t158 * t13 - t41 * t139 - t29 * t78 + t5 * t32 + t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t99 * t93, -t152 * t97 ^ 2, -t152 * t97 * t54 + t165 + t42, 0, 0, 0, 0, 0, 0, 0.2e1 * t45 * qJD(4) - t121, (-t43 - t141) * qJD(4) + t142, -t40 - t164, t12 * t45 + t13 * t43 + t165 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t40 + t164, (t43 - t141) * qJD(4) + t142, -t159, t121, qJDD(4), -g(3) * t87 + t166 * t86 - t41 * t45 + t125, g(3) * t86 + t166 * t87 + t41 * t43 - t115, 0, 0;];
tau_reg = t2;
