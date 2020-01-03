% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:25
% EndTime: 2019-12-31 17:32:26
% DurationCPUTime: 0.43s
% Computational Cost: add. (374->108), mult. (811->148), div. (0->0), fcn. (647->10), ass. (0->76)
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t24 = t49 * t53 + t50 * t51;
t99 = t24 * qJD(3);
t100 = qJD(5) * t99;
t91 = t49 ^ 2 + t50 ^ 2;
t19 = t24 * qJD(5);
t54 = cos(qJ(3));
t86 = t54 * qJD(2);
t73 = qJD(4) - t86;
t52 = sin(qJ(3));
t88 = qJD(2) * t52;
t67 = t91 * (qJD(3) * qJ(4) + t88);
t98 = -t67 * t54 - (-qJD(3) * pkin(3) + t73) * t52;
t97 = qJD(5) ^ 2;
t95 = t49 * t51;
t94 = t53 * t50;
t55 = qJD(3) ^ 2;
t93 = t54 * t55;
t92 = pkin(6) + qJ(4);
t90 = cos(pkin(7));
t89 = sin(pkin(7));
t87 = qJDD(3) * pkin(3);
t22 = -t94 + t95;
t85 = qJDD(5) * t22;
t84 = qJDD(5) * t24;
t83 = t49 * qJDD(3);
t82 = t50 * qJDD(1);
t81 = t50 * qJDD(3);
t80 = t52 * qJDD(2);
t79 = t54 * qJDD(2);
t75 = qJD(3) * t94;
t78 = qJD(5) * t75 + t51 * t81 + t53 * t83;
t77 = qJD(3) * t95;
t76 = qJD(5) * t95;
t40 = -pkin(4) * t50 - pkin(3);
t41 = qJD(3) * t88;
t74 = t91 * qJDD(3);
t11 = qJDD(3) * qJ(4) + t80 + (qJD(4) + t86) * qJD(3);
t6 = -qJDD(1) * t49 + t50 * t11;
t21 = -t89 * t52 - t90 * t54;
t23 = t90 * t52 - t89 * t54;
t72 = g(1) * t23 - g(2) * t21;
t71 = g(1) * t21 + g(2) * t23;
t5 = -t11 * t49 - t82;
t70 = -t49 * t5 + t50 * t6;
t69 = t51 * t83 - t53 * t81;
t68 = -g(1) * t89 + g(2) * t90;
t26 = t92 * t49;
t27 = t92 * t50;
t66 = -t26 * t53 - t27 * t51;
t65 = -t26 * t51 + t27 * t53;
t64 = qJDD(4) + t41 - t79;
t63 = t54 * qJDD(3) - t55 * t52;
t14 = t64 - t87;
t62 = -t14 + t72;
t18 = -qJD(5) * t94 + t76;
t61 = -qJD(5) * t18 + t84;
t60 = qJD(5) * t19 + t85;
t59 = qJD(5) * t22;
t57 = t70 + t71;
t56 = t41 + t62 + t87;
t48 = qJDD(1) - g(3);
t47 = pkin(8) + qJ(5);
t43 = cos(t47);
t42 = sin(t47);
t20 = t40 * qJD(3) + t73;
t15 = -t75 + t77;
t9 = t40 * qJDD(3) + t64;
t4 = pkin(6) * t81 + t6;
t3 = -t82 + (-pkin(6) * qJDD(3) - t11) * t49;
t2 = qJD(3) * t19 + t69;
t1 = -qJD(3) * t76 + t78;
t7 = [t48, t48, 0, 0, 0, 0, 0, 0, -t49 * t6 - t5 * t50 - g(3), 0, 0, 0, 0, 0, t60, t61; 0, qJDD(2) + t68, 0, t63, -qJDD(3) * t52 - t93, t63 * t50, -t63 * t49, t52 * t74 + t91 * t93, -t98 * qJD(3) - t14 * t54 + t70 * t52 + t68, 0, 0, 0, 0, 0, (-t2 - t100) * t54 + (qJD(3) * t15 + t22 * t97 - t84) * t52, (qJD(3) * t59 - t1) * t54 + (qJD(3) * t99 + t24 * t97 + t85) * t52; 0, 0, qJDD(3), t72 + t79, -t71 - t80, t56 * t50, -t56 * t49, t73 * qJD(3) * t91 + qJ(4) * t74 + t57, t62 * pkin(3) + t57 * qJ(4) + t98 * qJD(2) + t67 * qJD(4), t1 * t24 - t18 * t99, -t1 * t22 + t15 * t18 - t19 * t99 - t2 * t24, t61, -t60, 0, t40 * t2 + t9 * t22 + t20 * t19 + t66 * qJDD(5) + t72 * t43 + (-t24 * qJD(4) - t65 * qJD(5)) * qJD(5) + (-t52 * t15 + t54 * t19) * qJD(2), t40 * t1 + t9 * t24 - t20 * t18 - t65 * qJDD(5) - t72 * t42 + (t22 * qJD(4) - t66 * qJD(5)) * qJD(5) + (-t52 * t99 - t54 * t59) * qJD(2); 0, 0, 0, 0, 0, -t81, t83, -t91 * t55, -t67 * qJD(3) - t62, 0, 0, 0, 0, 0, t69 + 0.2e1 * t100, (-t15 - t77) * qJD(5) + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t15, -t15 ^ 2 + t99 ^ 2, (t15 - t77) * qJD(5) + t78, -t69, qJDD(5), g(3) * t43 - t20 * t99 + t53 * t3 - t51 * t4 - t71 * t42, -g(3) * t42 + t20 * t15 - t51 * t3 - t53 * t4 - t71 * t43;];
tau_reg = t7;
