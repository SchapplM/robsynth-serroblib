% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:42
% EndTime: 2019-12-31 16:51:43
% DurationCPUTime: 0.53s
% Computational Cost: add. (776->128), mult. (1205->149), div. (0->0), fcn. (621->6), ass. (0->83)
t84 = qJD(1) - qJD(3);
t112 = t84 ^ 2;
t51 = sin(qJ(3));
t56 = qJD(4) ^ 2;
t46 = qJDD(1) - qJDD(3);
t53 = cos(qJ(3));
t95 = t53 * t46;
t113 = (t56 + t112) * t51 + t95;
t104 = t84 * pkin(3);
t55 = -pkin(1) - pkin(2);
t28 = t55 * qJD(1) + qJD(2);
t87 = qJ(2) * qJD(1);
t12 = t53 * t28 - t51 * t87;
t6 = -t12 + t104;
t111 = t84 * t6;
t13 = t51 * t28 + t53 * t87;
t99 = t13 * t84;
t100 = t84 * t12;
t110 = -qJD(3) * t87 + t55 * qJDD(1) + qJDD(2);
t88 = qJD(4) * t84;
t85 = qJD(1) * qJD(2);
t86 = qJ(2) * qJDD(1);
t108 = qJD(3) * t28 + t85 + t86;
t7 = -pkin(6) * t84 + t13;
t50 = sin(qJ(4));
t48 = t50 ^ 2;
t52 = cos(qJ(4));
t49 = t52 ^ 2;
t91 = t48 + t49;
t107 = t7 * t91;
t24 = t53 * qJ(2) + t51 * t55;
t106 = t46 * pkin(3);
t105 = t46 * pkin(6);
t103 = sin(qJ(1));
t23 = -t51 * qJ(2) + t53 * t55;
t10 = t53 * qJD(2) + t23 * qJD(3);
t102 = t10 * t84;
t11 = t51 * qJD(2) + qJD(3) * t24;
t101 = t11 * t84;
t98 = t50 * t52;
t97 = t51 * t46;
t96 = t52 * t46;
t54 = cos(qJ(1));
t94 = t54 * pkin(1) + t103 * qJ(2);
t93 = g(1) * t103 - g(2) * t54;
t92 = t48 - t49;
t90 = pkin(1) * qJDD(1);
t83 = t112 * t98;
t82 = t54 * pkin(2) + t94;
t81 = 0.2e1 * t85;
t3 = t108 * t53 + t110 * t51;
t1 = t3 - t105;
t80 = t91 * t1;
t77 = qJDD(2) - t90;
t76 = -0.2e1 * t88 * t98;
t75 = -t103 * pkin(1) + t54 * qJ(2);
t74 = t108 * t51 - t110 * t53;
t16 = -t103 * t51 - t54 * t53;
t17 = -t103 * t53 + t54 * t51;
t73 = g(1) * t17 - g(2) * t16;
t72 = g(1) * t16 + g(2) * t17;
t2 = t74 + t106;
t71 = -t2 + t73;
t68 = g(1) * t54 + g(2) * t103;
t67 = -t1 - t72 + t111;
t66 = -t103 * pkin(2) + t75;
t65 = t84 * t91;
t64 = -pkin(6) * qJDD(4) + (t12 + t6 + t104) * qJD(4);
t63 = t73 - t74;
t18 = pkin(3) - t23;
t19 = -pkin(6) + t24;
t62 = -qJDD(4) * t19 + (-t18 * t84 - t10 - t6) * qJD(4);
t61 = -qJDD(4) * t51 + 0.2e1 * t53 * t88;
t60 = pkin(6) * t56 + t106 - t71 + t99;
t59 = -t18 * t46 + t19 * t56 - t101 + t71;
t58 = -t3 - t72;
t57 = qJD(1) ^ 2;
t26 = qJDD(4) * t52 - t56 * t50;
t25 = qJDD(4) * t50 + t56 * t52;
t9 = t49 * t46 + t76;
t8 = -t48 * t46 + t76;
t5 = -t50 * t96 + t92 * t88;
t4 = [0, 0, 0, 0, 0, qJDD(1), t93, t68, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t90 + t93, 0, -t68 + t81 + 0.2e1 * t86, -t77 * pkin(1) - g(1) * t75 - g(2) * t94 + (t81 + t86) * qJ(2), 0, 0, 0, 0, 0, t46, -t23 * t46 + t101 - t63, t24 * t46 + t102 - t58, 0, -g(1) * t66 - g(2) * t82 + t13 * t10 - t12 * t11 - t23 * t74 + t3 * t24, -t8, -0.2e1 * t5, -t25, t9, -t26, 0, t62 * t50 - t59 * t52, t59 * t50 + t62 * t52, -t72 + t91 * (-t19 * t46 - t1 - t102), t2 * t18 + t6 * t11 - g(1) * (t17 * pkin(3) + t16 * pkin(6) + t66) - g(2) * (-t16 * pkin(3) + t17 * pkin(6) + t82) + t10 * t107 + t19 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t57, -t57 * qJ(2) + t77 - t93, 0, 0, 0, 0, 0, 0, -t112 * t51 - t95, -t112 * t53 + t97, 0, (-t74 - t99) * t53 + (t3 + t100) * t51 - t93, 0, 0, 0, 0, 0, 0, -t113 * t52 + t61 * t50, t113 * t50 + t61 * t52, t53 * t65 * t84 - t91 * t97, (t80 - t111) * t51 + (-t65 * t7 - t2) * t53 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t63 - t99, t58 - t100, 0, 0, t8, 0.2e1 * t5, t25, -t9, t26, 0, t64 * t50 - t60 * t52, t60 * t50 + t64 * t52, t72 + t91 * (t1 + t100 - t105), -t6 * t13 - t12 * t107 + t71 * pkin(3) + (t80 + t72) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t92 * t112, -t50 * t46, t83, -t96, qJDD(4), g(3) * t52 + t67 * t50, -g(3) * t50 + t67 * t52, 0, 0;];
tau_reg = t4;
