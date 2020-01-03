% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRP3
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:17
% DurationCPUTime: 0.51s
% Computational Cost: add. (1133->108), mult. (1607->125), div. (0->0), fcn. (808->6), ass. (0->79)
t74 = qJD(1) + qJD(2);
t72 = t74 ^ 2;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t64 = t83 * t72 * t80;
t59 = qJDD(3) - t64;
t111 = t83 * t59;
t77 = t80 ^ 2;
t118 = t77 * t72;
t86 = qJD(3) ^ 2;
t60 = t86 + t118;
t38 = -t80 * t60 + t111;
t104 = qJD(3) * t74;
t73 = qJDD(1) + qJDD(2);
t113 = t80 * t73;
t48 = 0.2e1 * t83 * t104 + t113;
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t129 = pkin(1) * (t81 * t38 + t84 * t48);
t128 = pkin(2) * t48 + pkin(6) * t38;
t126 = qJ(4) * t48;
t122 = t83 * g(3);
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t95 = t85 * g(1) + t82 * g(2);
t57 = -qJD(1) ^ 2 * pkin(1) - t95;
t99 = t82 * g(1) - t85 * g(2);
t91 = qJDD(1) * pkin(1) + t99;
t28 = t84 * t57 + t81 * t91;
t23 = -t72 * pkin(2) + t73 * pkin(6) + t28;
t16 = t80 * t23 + t122;
t17 = -t80 * g(3) + t83 * t23;
t6 = t80 * t16 + t83 * t17;
t78 = t83 ^ 2;
t117 = t78 * t72;
t34 = t111 + t80 * (-t86 + t117);
t105 = t80 * qJ(4);
t94 = -t83 * pkin(3) - t105;
t120 = t94 * t72;
t12 = -qJDD(3) * pkin(3) - t86 * qJ(4) + (t23 + t120) * t80 + qJDD(4) + t122;
t125 = 2 * qJD(4);
t90 = qJDD(3) * qJ(4) + (qJD(3) * t125) + t83 * t120 + t17;
t11 = -t86 * pkin(3) + t90;
t2 = t83 * t11 + t80 * t12;
t119 = t74 * t80;
t101 = t80 * t104;
t110 = t83 * t73;
t27 = -t81 * t57 + t84 * t91;
t22 = -t73 * pkin(2) - t72 * pkin(6) - t27;
t88 = t22 - (-t101 + t110) * pkin(3) - t126;
t8 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t119 + t88;
t124 = -pkin(2) * t8 + pkin(6) * t2;
t121 = -pkin(2) * t22 + pkin(6) * t6;
t58 = qJDD(3) + t64;
t116 = t80 * t58;
t112 = t83 * t48;
t62 = -t86 - t117;
t37 = t83 * t62 - t116;
t49 = -0.2e1 * t101 + t110;
t108 = pkin(2) * t49 + pkin(6) * t37;
t106 = t77 + t78;
t52 = t106 * t73;
t55 = t106 * t72;
t107 = pkin(2) * t55 + pkin(6) * t52;
t103 = t80 * t22 - t128;
t102 = -t83 * t22 + t108;
t97 = t80 * (qJ(4) * t55 + t12) + t83 * ((t55 - t86) * pkin(3) + t90) + t107;
t96 = t107 + t6;
t24 = t80 * t49 + t112;
t87 = t119 * t125 - t88;
t92 = pkin(3) * t112 + t80 * (-pkin(3) * t101 + t126 + t87) + t128;
t89 = t49 * t105 + t108 + t83 * ((t49 - t101) * pkin(3) + t87);
t56 = (t77 - t78) * t72;
t36 = t116 + t83 * (t86 - t118);
t30 = t48 * t80;
t29 = t49 * t83;
t26 = pkin(1) * (t81 * t52 + t84 * t55);
t15 = pkin(1) * (t81 * t37 + t84 * t49);
t1 = [0, 0, 0, 0, 0, qJDD(1), t99, t95, 0, 0, 0, 0, 0, 0, 0, t73, pkin(1) * (-t81 * t72 + t84 * t73) + t27, pkin(1) * (-t84 * t72 - t81 * t73) - t28, 0, pkin(1) * (t84 * t27 + t81 * t28), t30, t24, t36, t29, t34, 0, t15 + t102, t103 - t129, t26 + t96, pkin(1) * (-t84 * t22 + t81 * t6) + t121, t30, t36, -t24, 0, -t34, t29, t15 + t89, t26 + t97, t92 + t129, pkin(1) * t81 * t2 + (-pkin(1) * t84 + t94) * t8 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t27, -t28, 0, 0, t30, t24, t36, t29, t34, 0, t102, t103, t96, t121, t30, t36, -t24, 0, -t34, t29, t89, t97, t92, t94 * t8 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t56, t113, t64, t110, qJDD(3), -t16, -t17, 0, 0, -t64, t113, -t56, qJDD(3), -t110, t64, pkin(3) * t58 + qJ(4) * t62 - t12, (-pkin(3) * t80 + qJ(4) * t83) * t73, qJ(4) * t59 + (t60 - t86) * pkin(3) + t90, -pkin(3) * t12 + qJ(4) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t113, -t60, t12;];
tauJ_reg = t1;
