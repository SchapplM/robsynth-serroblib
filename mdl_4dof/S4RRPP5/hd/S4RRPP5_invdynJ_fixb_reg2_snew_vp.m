% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:39
% DurationCPUTime: 0.50s
% Computational Cost: add. (641->139), mult. (1465->135), div. (0->0), fcn. (634->4), ass. (0->81)
t65 = sin(qJ(2));
t93 = qJD(1) * qJD(2);
t53 = t65 * t93;
t67 = cos(qJ(2));
t94 = t67 * qJDD(1);
t31 = -t53 + t94;
t125 = t31 - t53;
t118 = 2 * qJD(3);
t89 = t67 * t93;
t95 = t65 * qJDD(1);
t30 = t89 + t95;
t97 = qJD(1) * t65;
t38 = pkin(3) * t97 - qJD(2) * qJ(4);
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t91 = t66 * g(1) - t68 * g(2);
t79 = qJDD(1) * pkin(1) + t91;
t77 = t125 * pkin(2) + t79;
t124 = t30 * qJ(3) + t31 * qJ(4) + ((qJ(3) * qJD(2) + (2 * qJD(4))) * t67 + (t118 + t38) * t65) * qJD(1) + t77;
t70 = qJD(1) ^ 2;
t92 = t67 * t70 * t65;
t39 = qJDD(2) + t92;
t116 = t67 * g(3);
t69 = qJD(2) ^ 2;
t80 = -qJDD(2) * pkin(2) - t69 * qJ(3) + qJDD(3) + t116;
t75 = -0.2e1 * qJD(4) * qJD(2) + t80 - qJ(4) * t39 + (t30 - t89) * pkin(3);
t86 = t68 * g(1) + t66 * g(2);
t96 = qJDD(1) * pkin(5);
t20 = -t70 * pkin(1) - t86 + t96;
t99 = t65 * qJ(3);
t84 = -t67 * pkin(2) - t99;
t98 = t70 * t84;
t88 = t20 + t98;
t81 = t88 * t65;
t2 = t81 + t75;
t40 = qJDD(2) - t92;
t106 = t67 * t40;
t29 = 0.2e1 * t89 + t95;
t62 = t65 ^ 2;
t113 = t62 * t70;
t44 = -t69 - t113;
t104 = pkin(1) * t29 + pkin(5) * (t65 * t44 + t106);
t123 = -pkin(2) * t44 + qJ(3) * t40;
t63 = t67 ^ 2;
t112 = t63 * t70;
t115 = t70 * pkin(5);
t122 = pkin(3) * t112 + t115 + t124;
t121 = t106 + t65 * (-t69 + t112);
t119 = qJD(2) * t38 + qJDD(4);
t110 = t65 * t39;
t32 = -0.2e1 * t53 + t94;
t46 = -t69 - t112;
t114 = pkin(1) * t32 - pkin(5) * (-t67 * t46 + t110);
t117 = pkin(2) * t32;
t111 = t65 * t32;
t107 = t67 * t29;
t105 = pkin(2) + qJ(4);
t16 = -t65 * g(3) + t67 * t20;
t102 = t62 + t63;
t36 = t102 * t70;
t103 = pkin(1) * t36 + t102 * t96;
t101 = qJ(3) * t36;
t100 = qJ(3) * t46;
t15 = t65 * t20 + t116;
t90 = t65 * t15 + t67 * t16;
t87 = qJDD(2) * qJ(3) + qJD(2) * t118 + t67 * t98 + t16;
t83 = t30 + t89;
t7 = t107 + t111;
t14 = t67 * (t69 - t113) + t110;
t5 = -t69 * pkin(2) + t87;
t78 = (t36 - t69) * pkin(2) + t87;
t6 = t81 + t80;
t76 = t31 * pkin(3) + t119 + t5;
t72 = t97 * t118 + t115 + t77;
t54 = qJ(3) * t94;
t37 = (t62 - t63) * t70;
t19 = t79 + t115;
t12 = t83 * t65;
t11 = t125 * t67;
t3 = -qJ(4) * t112 + t76;
t1 = [0, 0, 0, 0, 0, qJDD(1), t91, t86, 0, 0, t12, t7, t14, t11, t121, 0, t67 * t19 + t114, -t65 * t19 - t104, t90 + t103, pkin(1) * t19 + pkin(5) * t90, 0, -t14, -t121, t12, t7, t11, t67 * t78 + (t6 + t101) * t65 + t103, t67 * (-t72 - t117) + (-t67 * t83 - t111) * qJ(3) - t114, t65 * t72 + pkin(2) * t107 + (t29 + t83) * t99 + t104, pkin(5) * (t67 * t5 + t65 * t6) + (pkin(1) - t84) * (t83 * qJ(3) + t72), 0, -t121, t14, t11, -t7, t12, t67 * ((t36 - t112) * qJ(4) + (t31 + t94) * pkin(3) + t78 + t119) + (t101 + (pkin(3) * qJDD(1) + t88) * t65 + t75) * t65 + t103, t67 * (pkin(3) * t40 + t105 * t29) + t104 + (pkin(3) * t44 + qJ(3) * t29 + t122) * t65, t65 * (-pkin(3) * t39 + qJ(3) * t32) + t114 + (pkin(3) * t46 + qJ(4) * t32 + t117 + t122) * t67, (t67 * t105 + pkin(1) + t99) * ((pkin(3) * t63 + pkin(5)) * t70 + t124) + (pkin(5) + pkin(3)) * (t65 * t2 + t67 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t37, t95, t92, t94, qJDD(2), -t15, -t16, 0, 0, qJDD(2), -t95, -t94, -t92, t37, t92, -pkin(2) * t95 + t54, -pkin(2) * t39 - t100 + t6, t5 + t123, -pkin(2) * t6 + qJ(3) * t5, qJDD(2), -t94, t95, t92, -t37, -t92, -t105 * t95 + t54, (-t44 - t112) * qJ(4) + t76 + t123, t105 * t39 + t100 - t2, qJ(3) * t3 - t105 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t39, t44, t6, 0, 0, 0, 0, 0, 0, t95, t44, -t39, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t40, t46, t3;];
tauJ_reg = t1;
