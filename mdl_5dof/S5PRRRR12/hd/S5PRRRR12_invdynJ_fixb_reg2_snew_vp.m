% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:21
% EndTime: 2024-09-28 18:08:22
% DurationCPUTime: 0.88s
% Computational Cost: add. (8517->149), mult. (12102->249), div. (0->0), fcn. (9536->14), ass. (0->131)
t109 = -g(3) + qJDD(1);
t112 = sin(pkin(5));
t115 = cos(pkin(5));
t110 = sin(pkin(11));
t113 = cos(pkin(11));
t140 = t110 * g(1) - t113 * g(2);
t169 = t112 * t109 + t115 * t140;
t118 = sin(qJ(3));
t122 = cos(qJ(3));
t106 = qJD(2) + qJD(3);
t102 = qJD(4) + t106;
t100 = t102 ^ 2;
t104 = qJDD(2) + qJDD(3);
t101 = qJDD(4) + t104;
t117 = sin(qJ(4));
t121 = cos(qJ(4));
t134 = -t121 * t100 - t117 * t101;
t87 = t117 * t100 - t121 * t101;
t66 = t118 * t87 + t122 * t134;
t168 = -t118 * t134 + t122 * t87;
t111 = sin(pkin(6));
t114 = cos(pkin(6));
t127 = -t115 * t109 + t112 * t140;
t163 = pkin(10) * t111;
t119 = sin(qJ(2));
t123 = cos(qJ(2));
t98 = -t113 * g(1) - t110 * g(2);
t70 = -t119 * t98 + t169 * t123;
t68 = qJDD(2) * pkin(2) + t70;
t124 = qJD(2) ^ 2;
t71 = t169 * t119 + t123 * t98;
t69 = -t124 * pkin(2) + t71;
t53 = -t118 * t69 + t122 * t68;
t51 = t104 * pkin(3) + t53;
t103 = t106 ^ 2;
t54 = t118 * t68 + t122 * t69;
t52 = -t103 * pkin(3) + t54;
t30 = -t117 * t52 + t121 * t51;
t128 = t101 * pkin(4) + t100 * t163 + t30;
t165 = -t111 * t127 + t114 * t128;
t97 = t114 * t102 + qJD(5);
t96 = t97 ^ 2;
t105 = t111 ^ 2;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t31 = t117 * t51 + t121 * t52;
t29 = -t100 * pkin(4) + t101 * t163 + t31;
t15 = t116 * t29 - t165 * t120;
t16 = t165 * t116 + t120 * t29;
t139 = t116 * t16 - t120 * t15;
t24 = t111 * t128 + t114 * t127;
t6 = t111 * t24 + t139 * t114;
t8 = t116 * t15 + t120 * t16;
t164 = pkin(4) * t6 + t8 * t163;
t162 = t102 * t97;
t156 = t105 * t120;
t93 = t116 * t100 * t156;
t95 = t114 * t101 + qJDD(5);
t74 = t93 + t95;
t161 = t116 * t74;
t75 = -t93 + t95;
t160 = t120 * t75;
t159 = t100 * t105;
t158 = t101 * t116;
t157 = t102 * t120;
t155 = t111 * t116;
t154 = t111 * t120;
t3 = t117 * t8 + t121 * t6;
t152 = pkin(3) * t3 + t164;
t151 = t102 * t155;
t107 = t116 ^ 2;
t150 = t107 * t159;
t108 = t120 ^ 2;
t149 = t108 * t159;
t130 = qJD(5) * t157 + t158;
t129 = t130 * t111;
t82 = t154 * t162;
t62 = -t82 + t129;
t141 = -qJD(5) * t151 + t101 * t154;
t81 = t97 * t151;
t63 = t141 + t81;
t138 = t116 * t63 - t120 * t62;
t80 = (-t107 - t108) * t159;
t41 = -t111 * t80 + t138 * t114;
t49 = t116 * t62 + t120 * t63;
t148 = pkin(4) * t41 + t15 * t155 + t16 * t154 + t49 * t163;
t72 = -t150 - t96;
t137 = -t116 * t75 + t120 * t72;
t61 = (t158 + (qJD(5) + t97) * t157) * t111;
t43 = -t111 * t61 + t137 * t114;
t58 = -t116 * t72 - t160;
t147 = pkin(4) * t43 - t114 * t16 - t24 * t155 + t58 * t163;
t78 = -t96 - t149;
t136 = t116 * t78 + t120 * t74;
t64 = t141 - t81;
t47 = t111 * t64 + t136 * t114;
t60 = t120 * t78 - t161;
t146 = pkin(4) * t47 - t114 * t15 + t24 * t154 + t60 * t163;
t145 = -pkin(3) * t87 + t30;
t26 = t117 * t49 + t121 * t41;
t144 = pkin(3) * t26 + t148;
t33 = t117 * t58 + t121 * t43;
t143 = pkin(3) * t33 + t147;
t37 = t117 * t60 + t121 * t47;
t142 = pkin(3) * t37 + t146;
t135 = t100 * t114 - t162;
t91 = -t122 * t103 - t118 * t104;
t133 = t118 * t103 - t122 * t104;
t132 = pkin(3) * t134 - t31;
t90 = t114 * t95;
t79 = (t107 - t108) * t159;
t77 = t115 * t127;
t56 = (t130 * t105 - t135 * t156) * t116;
t55 = (t135 * t116 * t105 + t141 * t111) * t120;
t46 = t114 * t63 + (t116 * (-t96 + t149) + t160) * t111;
t45 = t114 * t62 + (t161 + t120 * (t96 - t150)) * t111;
t40 = t114 * t79 + (t116 * t64 + (t129 + t82) * t120) * t111;
t38 = -t117 * t47 + t121 * t60;
t35 = t118 * t54 + t122 * t53;
t34 = -t117 * t43 + t121 * t58;
t27 = -t117 * t41 + t121 * t49;
t21 = t118 * t38 + t122 * t37;
t20 = t118 * t34 + t122 * t33;
t19 = -t117 * t30 + t121 * t31;
t18 = t117 * t31 + t121 * t30;
t17 = pkin(3) * t18;
t10 = t118 * t27 + t122 * t26;
t9 = t118 * t19 + t122 * t18;
t4 = -t117 * t6 + t121 * t8;
t1 = t118 * t4 + t122 * t3;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, (qJDD(2) * t123 - t119 * t124) * t112, (-qJDD(2) * t119 - t123 * t124) * t112, 0, -t77 + (t119 * t71 + t123 * t70) * t112, 0, 0, 0, 0, 0, 0, (t119 * t91 - t123 * t133) * t112, (t119 * t133 + t123 * t91) * t112, 0, -t77 + (t119 * (-t118 * t53 + t122 * t54) + t123 * t35) * t112, 0, 0, 0, 0, 0, 0, (t119 * t66 - t123 * t168) * t112, (t119 * t168 + t123 * t66) * t112, 0, -t77 + (t119 * (-t118 * t18 + t122 * t19) + t123 * t9) * t112, 0, 0, 0, 0, 0, 0, t115 * (t136 * t111 - t114 * t64) + (t119 * (-t118 * t37 + t122 * t38) + t123 * t21) * t112, t115 * (t137 * t111 + t114 * t61) + (t119 * (-t118 * t33 + t122 * t34) + t123 * t20) * t112, t115 * (t138 * t111 + t114 * t80) + (t119 * (-t118 * t26 + t122 * t27) + t123 * t10) * t112, t115 * (t139 * t111 - t114 * t24) + (t119 * (-t118 * t3 + t122 * t4) + t123 * t1) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t70, -t71, 0, 0, 0, 0, 0, 0, 0, t104, -pkin(2) * t133 + t53, pkin(2) * t91 - t54, 0, pkin(2) * t35, 0, 0, 0, 0, 0, t101, -pkin(2) * t168 + t145, pkin(2) * t66 + t132, 0, pkin(2) * t9 + t17, t56, t40, t45, t55, t46, t90, pkin(2) * t21 + t142, pkin(2) * t20 + t143, pkin(2) * t10 + t144, pkin(2) * t1 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t53, -t54, 0, 0, 0, 0, 0, 0, 0, t101, t145, t132, 0, t17, t56, t40, t45, t55, t46, t90, t142, t143, t144, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t30, -t31, 0, 0, t56, t40, t45, t55, t46, t90, t146, t147, t148, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t79, t62, t93, t63, t95, -t15, -t16, 0, 0;];
tauJ_reg = t2;
