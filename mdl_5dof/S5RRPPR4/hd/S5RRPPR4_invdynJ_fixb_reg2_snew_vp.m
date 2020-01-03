% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:51
% DurationCPUTime: 0.54s
% Computational Cost: add. (2647->115), mult. (3262->147), div. (0->0), fcn. (1556->8), ass. (0->87)
t129 = pkin(2) + pkin(3);
t87 = (qJD(1) + qJD(2));
t127 = t87 ^ 2;
t86 = qJDD(1) + qJDD(2);
t94 = sin(qJ(2));
t97 = cos(qJ(2));
t128 = pkin(1) * (t127 * t97 + t94 * t86);
t126 = 2 * t87;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t67 = t96 * t127 * t93;
t63 = qJDD(5) + t67;
t125 = t93 * t63;
t124 = t93 * t86;
t64 = qJDD(5) - t67;
t123 = t96 * t64;
t95 = sin(qJ(1));
t98 = cos(qJ(1));
t118 = t95 * g(1) - t98 * g(2);
t61 = qJDD(1) * pkin(1) + t118;
t112 = t98 * g(1) + t95 * g(2);
t62 = -qJD(1) ^ 2 * pkin(1) - t112;
t37 = t97 * t61 - t94 * t62;
t111 = -qJDD(3) + t37;
t79 = t86 * pkin(2);
t31 = -qJ(3) * t127 - t111 - t79;
t103 = -t86 * pkin(3) + t31;
t38 = t94 * t61 + t97 * t62;
t119 = (qJD(3) * t126) + t38;
t78 = t86 * qJ(3);
t115 = t78 + t119;
t21 = -t129 * t127 + t115;
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t15 = t91 * t103 + t92 * t21;
t29 = -pkin(2) * t127 + t115;
t122 = -pkin(2) * t31 + qJ(3) * t29;
t121 = g(3) + qJDD(4);
t120 = qJD(5) * t126;
t117 = -t92 * t103 + t91 * t21;
t116 = -t127 * t92 + t91 * t86;
t114 = 0.2e1 * t78 + t119;
t8 = -t117 * t92 + t91 * t15;
t9 = t117 * t91 + t92 * t15;
t113 = qJ(3) * t9 - t129 * t8;
t71 = t96 * t86;
t51 = t93 * t120 - t71;
t13 = -pkin(4) * t127 - t86 * pkin(7) + t15;
t10 = -t96 * t121 + t93 * t13;
t11 = t93 * t121 + t96 * t13;
t5 = t93 * t10 + t96 * t11;
t109 = t127 * t91 + t92 * t86;
t108 = t111 + 0.2e1 * t79;
t107 = qJ(3) * t109 - t129 * t116 + t15;
t106 = t96 * t120 + t124;
t105 = qJ(3) * t116 + t129 * t109 + t117;
t12 = t86 * pkin(4) - pkin(7) * t127 + t117;
t3 = -t92 * t12 + t91 * t5;
t4 = t91 * t12 + t92 * t5;
t104 = pkin(4) * t12 - pkin(7) * t5 + qJ(3) * t4 - t129 * t3;
t90 = t96 ^ 2;
t73 = t90 * t127;
t99 = qJD(5) ^ 2;
t66 = -t73 - t99;
t42 = t96 * t66 - t125;
t24 = t91 * t42 + t92 * t51;
t26 = t92 * t42 - t91 * t51;
t102 = -pkin(4) * t51 - pkin(7) * t42 + qJ(3) * t26 + t96 * t12 - t129 * t24;
t89 = t93 ^ 2;
t72 = t89 * t127;
t65 = -t72 - t99;
t44 = -t93 * t65 - t123;
t25 = t106 * t92 + t91 * t44;
t27 = -t106 * t91 + t92 * t44;
t101 = -pkin(4) * t106 - pkin(7) * t44 + qJ(3) * t27 - t93 * t12 - t129 * t25;
t59 = (-t89 - t90) * t86;
t60 = t72 + t73;
t35 = t91 * t59 + t92 * t60;
t36 = t92 * t59 - t91 * t60;
t100 = -pkin(4) * t60 - pkin(7) * t59 + qJ(3) * t36 - t129 * t35 - t5;
t49 = pkin(1) * (-t127 * t94 + t97 * t86);
t43 = -t123 - t93 * (t73 - t99);
t41 = -t96 * (-t72 + t99) - t125;
t40 = t51 * t96;
t39 = t106 * t93;
t32 = t106 * t96 - t93 * t51;
t1 = [0, 0, 0, 0, 0, qJDD(1), t118, t112, 0, 0, 0, 0, 0, 0, 0, t86, t37 + t49, -t38 - t128, 0, pkin(1) * (t97 * t37 + t94 * t38), 0, 0, 0, t86, 0, 0, t108 + t49, 0, t114 + t128, pkin(1) * (t94 * t29 - t97 * t31) + t122, 0, 0, 0, 0, 0, t86, pkin(1) * (t109 * t97 + t116 * t94) + t105, pkin(1) * (t109 * t94 - t116 * t97) + t107, 0, pkin(1) * (-t97 * t8 + t94 * t9) + t113, t39, t32, t41, -t40, t43, 0, pkin(1) * (-t97 * t24 + t94 * t26) + t102, pkin(1) * (-t97 * t25 + t94 * t27) + t101, pkin(1) * (-t97 * t35 + t94 * t36) + t100, pkin(1) * (-t97 * t3 + t94 * t4) + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t37, -t38, 0, 0, 0, 0, 0, t86, 0, 0, t108, 0, t114, t122, 0, 0, 0, 0, 0, t86, t105, t107, 0, t113, t39, t32, t41, -t40, t43, 0, t102, t101, t100, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t127, t31, 0, 0, 0, 0, 0, 0, -t109, t116, 0, t8, 0, 0, 0, 0, 0, 0, t24, t25, t35, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, 0, 0, 0, 0, 0, t96 * t63 + t93 * t66, -t93 * t64 + t96 * t65, 0, -t96 * t10 + t93 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t72 - t73, -t124, t67, -t71, qJDD(5), -t10, -t11, 0, 0;];
tauJ_reg = t1;
