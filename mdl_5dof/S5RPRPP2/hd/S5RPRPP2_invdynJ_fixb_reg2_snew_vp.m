% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:15
% DurationCPUTime: 0.80s
% Computational Cost: add. (1155->175), mult. (2422->184), div. (0->0), fcn. (1142->6), ass. (0->104)
t86 = sin(qJ(3));
t87 = cos(qJ(3));
t89 = qJD(1) ^ 2;
t122 = t87 * t89 * t86;
t58 = qJDD(3) - t122;
t137 = t87 * t58;
t78 = t86 ^ 2;
t145 = t78 * t89;
t88 = qJD(3) ^ 2;
t62 = t88 + t145;
t29 = -t86 * t62 + t137;
t125 = qJD(1) * qJD(3);
t114 = t87 * t125;
t128 = t86 * qJDD(1);
t47 = 0.2e1 * t114 + t128;
t83 = sin(pkin(7));
t84 = cos(pkin(7));
t120 = pkin(1) * (t83 * t29 + t84 * t47) + pkin(2) * t47 + pkin(6) * t29;
t57 = qJDD(3) + t122;
t154 = pkin(4) * t57;
t130 = qJD(1) * t86;
t104 = -qJD(3) * pkin(4) - qJ(5) * t130;
t127 = t87 * qJDD(1);
t74 = t86 * t125;
t48 = -t74 + t127;
t153 = t48 * qJ(5) - qJD(3) * t104;
t79 = t87 ^ 2;
t144 = t79 * t89;
t65 = -t88 - t144;
t152 = pkin(3) * t57 + qJ(4) * t65;
t24 = t137 + t86 * (-t88 + t144);
t151 = t104 * t130 + qJDD(5);
t150 = t48 * pkin(4) + t151;
t129 = -g(3) + qJDD(2);
t73 = t87 * t129;
t109 = qJDD(3) * pkin(3) + t88 * qJ(4) - qJDD(4) + t73;
t132 = t86 * qJ(4);
t107 = -t87 * pkin(3) - t132;
t44 = t107 * qJD(1);
t131 = qJD(1) * t44;
t146 = sin(qJ(1));
t147 = cos(qJ(1));
t100 = t147 * g(1) + t146 * g(2);
t45 = -t89 * pkin(1) - t100;
t99 = t146 * g(1) - t147 * g(2);
t94 = qJDD(1) * pkin(1) + t99;
t136 = t84 * t45 + t83 * t94;
t17 = -t89 * pkin(2) + qJDD(1) * pkin(6) + t136;
t110 = t17 + t131;
t9 = t110 * t86 - t109;
t149 = pkin(3) + pkin(4);
t142 = t86 * t17;
t141 = t86 * t57;
t138 = t87 * t47;
t12 = t86 * t129 + t87 * t17;
t135 = t78 + t79;
t55 = t135 * t89;
t134 = qJ(4) * t55;
t133 = qJ(4) * t87;
t126 = qJ(5) * qJDD(1);
t124 = qJD(1) * qJD(5);
t123 = qJD(4) * qJD(3);
t28 = t87 * t65 - t141;
t49 = -0.2e1 * t74 + t127;
t121 = pkin(1) * (t83 * t28 + t84 * t49) + pkin(6) * t28 + pkin(2) * t49;
t54 = t135 * qJDD(1);
t119 = pkin(1) * (t83 * t54 + t84 * t55) + pkin(6) * t54 + pkin(2) * t55;
t118 = 0.2e1 * t124;
t117 = -pkin(1) * t84 - pkin(2);
t116 = pkin(1) * t83 + pkin(6);
t115 = qJ(5) * qJD(3) * t87;
t11 = -t73 + t142;
t6 = t86 * t11 + t87 * t12;
t113 = -t83 * t45 + t84 * t94;
t112 = t62 - t144;
t111 = qJDD(3) * qJ(4) + t87 * t131 + t12;
t106 = t86 * t49 + t138;
t26 = -t87 * (-t88 + t145) + t141;
t27 = t86 * t58 + t87 * t62;
t16 = -qJDD(1) * pkin(2) - t89 * pkin(6) - t113;
t105 = t88 * pkin(3) - t111;
t103 = -t87 * t149 + t117;
t101 = t114 + t128;
t102 = -t101 * qJ(5) - t109 - t154;
t75 = 0.2e1 * t123;
t8 = -t105 + t75;
t98 = -t48 * pkin(3) + t16 + (-t101 - t114) * qJ(4);
t97 = pkin(3) * t62 + qJ(4) * t58 + t8;
t96 = 0.2e1 * qJD(4) * t130 - t98;
t95 = t102 + t142;
t93 = pkin(4) * t144 + t105 + t153;
t92 = (pkin(3) * qJD(3) - 0.2e1 * qJD(4)) * t130 + t98;
t91 = (t49 - t74) * pkin(3) + t96;
t90 = -pkin(3) * t74 + qJ(4) * t47 + t96;
t70 = -0.2e1 * t87 * t124;
t69 = t86 * t118;
t56 = (t78 - t79) * t89;
t25 = t87 * t57 + t86 * t65;
t23 = t47 * t86;
t22 = (t48 - t74) * t87;
t5 = (t115 + (-0.2e1 * qJD(5) + t44) * t86) * qJD(1) + t95;
t4 = t70 + t75 - t93;
t3 = qJ(5) * t144 - t150 + t92;
t1 = [0, 0, 0, 0, 0, qJDD(1), t99, t100, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t84 * qJDD(1) - t83 * t89) + t113, pkin(1) * (-t83 * qJDD(1) - t84 * t89) - t136, 0, pkin(1) * (t84 * t113 + t83 * t136), t23, t106, t26, t22, t24, 0, -t87 * t16 + t121, t86 * t16 - t120, t6 + t119, -pkin(2) * t16 + pkin(6) * t6 + pkin(1) * (-t84 * t16 + t83 * t6), t23, t26, -t106, 0, -t24, t22, t49 * t132 + t87 * t91 + t121, t87 * (t75 + (t55 - t88) * pkin(3) + t111) + (t9 + t134) * t86 + t119, pkin(3) * t138 + t86 * t90 + t120, t116 * (t87 * t8 + t86 * t9) + (t107 + t117) * t92, t23, -t106, -t26, t22, t24, 0, t86 * (qJ(4) * t49 + qJ(5) * t57) + t87 * ((-t65 - t144) * qJ(5) + (t48 + t49) * pkin(4) + t91 + t151) + t121, t86 * (t112 * qJ(5) + t150 + t90) + t87 * (-qJ(5) * t58 + t149 * t47) + t120, (-0.2e1 * t123 + (t118 + t126) * t87 + t93) * t87 - t116 * t54 + t103 * t55 + (-qJ(5) * t114 - t134 + t69 + (-t110 + t126) * t86 - t102) * t86, (t103 - t132) * t3 + (t116 - qJ(5)) * (t87 * t4 + t86 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, 0, 0, 0, 0, 0, t25, -t27, 0, -t87 * t11 + t86 * t12, 0, 0, 0, 0, 0, 0, t25, 0, t27, t86 * t8 - t87 * t9, 0, 0, 0, 0, 0, 0, t25, t27, 0, t86 * t4 - t87 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t56, t128, t122, t127, qJDD(3), -t11, -t12, 0, 0, -t122, t128, -t56, qJDD(3), -t127, t122, -t9 + t152, (-pkin(3) * t86 + t133) * qJDD(1), t97, -pkin(3) * t9 + qJ(4) * t8, -t122, -t56, -t128, t122, t127, qJDD(3), t154 + t69 + (-t44 * t86 - t115) * qJD(1) - t95 + t152, t112 * pkin(4) - t153 + t70 + t97, (t149 * t86 - t133) * qJDD(1), qJ(4) * t4 - t149 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t128, -t62, t9, 0, 0, 0, 0, 0, 0, -t57, -t62, -t128, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t47, -t55, -t3;];
tauJ_reg = t1;
