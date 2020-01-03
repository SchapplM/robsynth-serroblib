% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:20
% EndTime: 2019-12-31 16:49:21
% DurationCPUTime: 0.45s
% Computational Cost: add. (475->109), mult. (992->169), div. (0->0), fcn. (667->12), ass. (0->83)
t62 = sin(pkin(7));
t46 = t62 * pkin(1) + pkin(5);
t105 = pkin(6) + t46;
t58 = qJ(1) + pkin(7);
t50 = sin(t58);
t51 = cos(t58);
t108 = g(1) * t51 + g(2) * t50;
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t87 = t105 * qJD(1);
t14 = t68 * qJD(2) - t87 * t65;
t57 = qJD(3) + qJD(4);
t15 = t65 * qJD(2) + t87 * t68;
t67 = cos(qJ(4));
t99 = qJD(1) * t68;
t90 = t67 * t99;
t100 = qJD(1) * t65;
t64 = sin(qJ(4));
t91 = t64 * t100;
t19 = -t90 + t91;
t21 = -t67 * t100 - t64 * t99;
t104 = t21 * t19;
t103 = t67 * t15;
t59 = t65 ^ 2;
t102 = -t68 ^ 2 + t59;
t101 = qJD(3) * pkin(3);
t63 = cos(pkin(7));
t47 = -t63 * pkin(1) - pkin(2);
t37 = qJD(1) * t47;
t98 = qJD(4) * t64;
t96 = qJDD(2) - g(3);
t95 = t65 * qJDD(1);
t94 = t68 * qJDD(1);
t93 = qJD(1) * qJD(3);
t92 = t65 * t101;
t89 = t68 * t93;
t88 = qJD(3) * t105;
t34 = t46 * qJDD(1);
t86 = pkin(6) * qJDD(1) + t34;
t85 = g(1) * t50 - g(2) * t51;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t84 = g(1) * t66 - g(2) * t69;
t83 = t64 * t95 - t67 * t94;
t25 = t64 * t65 - t67 * t68;
t10 = t57 * t25;
t26 = t64 * t68 + t67 * t65;
t56 = qJDD(3) + qJDD(4);
t82 = -t10 * t57 + t26 * t56;
t13 = t14 + t101;
t81 = -t64 * t13 - t103;
t23 = t105 * t65;
t24 = t105 * t68;
t80 = -t67 * t23 - t64 * t24;
t79 = -t64 * t23 + t67 * t24;
t31 = -t68 * pkin(3) + t47;
t77 = -t37 * qJD(1) + t108 - t34;
t76 = 0.2e1 * qJD(3) * t37 - qJDD(3) * t46;
t7 = qJD(4) * t90 - t57 * t91 + t64 * t94 + (t89 + t95) * t67;
t70 = qJD(3) ^ 2;
t75 = -0.2e1 * qJDD(1) * t47 - t46 * t70 + t85;
t22 = t31 * qJD(1);
t52 = t68 * qJDD(2);
t4 = qJDD(3) * pkin(3) - t15 * qJD(3) - t86 * t65 + t52;
t61 = qJ(3) + qJ(4);
t54 = sin(t61);
t55 = cos(t61);
t74 = t22 * t19 + t15 * t98 + g(3) * t54 + (-t15 * t57 - t4) * t64 + t108 * t55;
t11 = t57 * t26;
t5 = t14 * qJD(3) + t65 * qJDD(2) + t86 * t68;
t73 = -g(3) * t55 + t81 * qJD(4) + t108 * t54 + t22 * t21 + t67 * t4 - t64 * t5;
t71 = qJD(1) ^ 2;
t33 = qJDD(3) * t68 - t70 * t65;
t32 = qJDD(3) * t65 + t70 * t68;
t18 = t68 * t88;
t17 = t65 * t88;
t16 = qJD(1) * t92 + t31 * qJDD(1);
t9 = -t19 ^ 2 + t21 ^ 2;
t8 = t11 * qJD(1) + t83;
t6 = -t11 * t57 - t25 * t56;
t2 = -t83 + (-qJD(1) * t26 - t21) * t57;
t1 = t19 * t57 + t7;
t3 = [qJDD(1), t84, g(1) * t69 + g(2) * t66, (t84 + (t62 ^ 2 + t63 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t59 * qJDD(1) + 0.2e1 * t65 * t89, -0.2e1 * t102 * t93 + 0.2e1 * t65 * t94, t32, t33, 0, t65 * t76 + t68 * t75, -t65 * t75 + t68 * t76, t21 * t10 + t7 * t26, t10 * t19 + t21 * t11 - t7 * t25 - t26 * t8, t82, t6, 0, t19 * t92 + t31 * t8 + t16 * t25 + t22 * t11 + (-qJD(4) * t79 + t64 * t17 - t67 * t18) * t57 + t80 * t56 + t85 * t55, -t21 * t92 + t31 * t7 + t16 * t26 - t22 * t10 - (qJD(4) * t80 - t67 * t17 - t64 * t18) * t57 - t79 * t56 - t85 * t54; 0, 0, 0, t96, 0, 0, 0, 0, 0, t33, -t32, 0, 0, 0, 0, 0, t6, -t82; 0, 0, 0, 0, -t65 * t71 * t68, t102 * t71, t95, t94, qJDD(3), -g(3) * t68 + t65 * t77 + t52, -t65 * t96 + t68 * t77, -t104, t9, t1, t2, t56, -(-t64 * t14 - t103) * t57 + (-t100 * t19 + t67 * t56 - t57 * t98) * pkin(3) + t73, (-qJD(4) * t13 + t14 * t57 - t5) * t67 + (-qJD(4) * t67 * t57 + t100 * t21 - t64 * t56) * pkin(3) + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t9, t1, t2, t56, -t57 * t81 + t73, (-t5 + (-qJD(4) + t57) * t13) * t67 + t74;];
tau_reg = t3;
