% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:13
% EndTime: 2019-12-31 18:01:16
% DurationCPUTime: 0.61s
% Computational Cost: add. (2568->129), mult. (3841->159), div. (0->0), fcn. (1578->8), ass. (0->89)
t73 = (-qJD(1) + qJD(4));
t71 = t73 ^ 2;
t72 = qJDD(1) - qJDD(4);
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t50 = t81 * t71 + t84 * t72;
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t95 = -t84 * t71 + t81 * t72;
t30 = t78 * t50 + t79 * t95;
t111 = t79 * t50 - t78 * t95;
t108 = pkin(1) + pkin(2);
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t59 = t83 * t71 * t80;
t53 = qJDD(5) + t59;
t107 = t80 * t53;
t106 = t80 * t72;
t54 = qJDD(5) - t59;
t105 = t83 * t54;
t87 = qJD(1) ^ 2;
t74 = qJDD(1) * qJ(2);
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t98 = t85 * g(1) + t82 * g(2);
t94 = (2 * qJD(2) * qJD(1)) - t98;
t92 = t74 + t94;
t42 = -t108 * t87 + t92;
t104 = t82 * g(1) - t85 * g(2);
t101 = -qJDD(2) + t104;
t93 = -t87 * qJ(2) - t101;
t43 = -t108 * qJDD(1) + t93;
t28 = t79 * t42 + t78 * t43;
t22 = -t87 * pkin(3) + t28;
t96 = t78 * t42 - t79 * t43;
t89 = -qJDD(1) * pkin(3) - t96;
t14 = t84 * t22 + t81 * t89;
t103 = qJDD(1) * pkin(1);
t102 = 2 * qJD(5) * t73;
t12 = -t71 * pkin(4) - t72 * pkin(7) + t14;
t77 = g(3) + qJDD(3);
t10 = t83 * t12 + t80 * t77;
t9 = t80 * t12 - t83 * t77;
t5 = t83 * t10 + t80 * t9;
t100 = t81 * t22 - t84 * t89;
t11 = t72 * pkin(4) - t71 * pkin(7) + t100;
t99 = -pkin(4) * t11 + pkin(7) * t5;
t61 = t83 * t72;
t97 = t80 * t102 + t61;
t45 = t83 * t102 - t106;
t76 = t83 ^ 2;
t63 = t76 * t71;
t86 = qJD(5) ^ 2;
t58 = -t63 - t86;
t38 = t83 * t58 - t107;
t91 = -pkin(4) * t97 + pkin(7) * t38 - t83 * t11;
t75 = t80 ^ 2;
t62 = t75 * t71;
t57 = -t62 - t86;
t39 = -t80 * t57 - t105;
t90 = pkin(4) * t45 - pkin(7) * t39 - t80 * t11;
t47 = (-t75 - t76) * t72;
t52 = t62 + t63;
t88 = pkin(4) * t52 + pkin(7) * t47 + t5;
t56 = t79 * qJDD(1) + t78 * t87;
t55 = -t78 * qJDD(1) + t79 * t87;
t44 = -t93 + t103;
t37 = t107 + t83 * (-t62 + t86);
t36 = t80 * (t63 - t86) + t105;
t35 = t45 * t80;
t34 = t97 * t83;
t33 = t84 * t47 - t81 * t52;
t32 = t81 * t47 + t84 * t52;
t29 = t83 * t45 - t80 * t97;
t26 = t84 * t39 + t81 * t45;
t25 = t84 * t38 + t81 * t97;
t24 = t81 * t39 - t84 * t45;
t23 = t81 * t38 - t84 * t97;
t18 = t79 * t32 + t78 * t33;
t17 = t78 * t28 - t79 * t96;
t16 = t79 * t24 + t78 * t26;
t15 = t79 * t23 + t78 * t25;
t7 = t100 * t81 + t84 * t14;
t6 = -t100 * t84 + t81 * t14;
t4 = t81 * t11 + t84 * t5;
t3 = -t84 * t11 + t81 * t5;
t2 = t79 * t6 + t78 * t7;
t1 = t79 * t3 + t78 * t4;
t8 = [0, 0, 0, 0, 0, qJDD(1), t104, t98, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t101 + 0.2e1 * t103, 0, 0.2e1 * t74 + t94, qJ(2) * (-t87 * pkin(1) + t92) + pkin(1) * t44, 0, 0, 0, 0, 0, qJDD(1), -qJ(2) * t55 + t108 * t56 + t96, qJ(2) * t56 + t108 * t55 + t28, 0, qJ(2) * (t79 * t28 + t78 * t96) - t108 * t17, 0, 0, 0, 0, 0, t72, pkin(3) * t50 + qJ(2) * t30 + t108 * t111 + t100, -pkin(3) * t95 + qJ(2) * t111 - t108 * t30 + t14, 0, qJ(2) * (-t78 * t6 + t79 * t7) - pkin(3) * t6 - t108 * t2, -t35, -t29, -t37, t34, -t36, 0, qJ(2) * (-t78 * t23 + t79 * t25) - pkin(3) * t23 - t108 * t15 - t91, qJ(2) * (-t78 * t24 + t79 * t26) - pkin(3) * t24 - t108 * t16 + t90, qJ(2) * (-t78 * t32 + t79 * t33) - pkin(3) * t32 - t108 * t18 - t88, qJ(2) * (-t78 * t3 + t79 * t4) - pkin(3) * t3 - t108 * t1 - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t87, -t44, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, t17, 0, 0, 0, 0, 0, 0, -t111, t30, 0, t2, 0, 0, 0, 0, 0, 0, t15, t16, t18, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, 0, 0, 0, 0, t83 * t53 + t80 * t58, -t80 * t54 + t83 * t57, 0, t80 * t10 - t83 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t100, -t14, 0, 0, t35, t29, t37, -t34, t36, 0, t91, -t90, t88, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t62 - t63, -t106, t59, -t61, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg = t8;
