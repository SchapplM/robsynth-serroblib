% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR8
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
% tau_reg [4x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:17
% EndTime: 2019-12-31 16:55:18
% DurationCPUTime: 0.48s
% Computational Cost: add. (445->121), mult. (873->170), div. (0->0), fcn. (565->8), ass. (0->84)
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t105 = -g(1) * t56 + g(2) * t59;
t107 = qJDD(2) + t105;
t55 = sin(qJ(3));
t84 = qJD(1) * qJD(3);
t58 = cos(qJ(3));
t86 = t58 * qJDD(1);
t106 = t55 * t84 - t86;
t94 = pkin(1) * qJDD(1);
t104 = t94 - t107;
t50 = qJD(3) + qJD(4);
t60 = -pkin(1) - pkin(5);
t101 = pkin(6) - t60;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t23 = t54 * t58 + t55 * t57;
t10 = t50 * t23;
t24 = -t54 * t55 + t57 * t58;
t49 = qJDD(3) + qJDD(4);
t100 = -t10 * t50 + t24 * t49;
t93 = qJD(1) * t55;
t82 = t54 * t93;
t92 = qJD(1) * t58;
t18 = -t57 * t92 + t82;
t19 = t23 * qJD(1);
t99 = t18 * t19;
t35 = qJD(1) * t60 + qJD(2);
t15 = -pkin(6) * t93 + t35 * t55;
t98 = t57 * t15;
t52 = t58 ^ 2;
t97 = t55 ^ 2 - t52;
t61 = qJD(3) ^ 2;
t62 = qJD(1) ^ 2;
t96 = -t61 - t62;
t95 = t62 * qJ(2);
t91 = qJD(3) * t55;
t90 = qJD(3) * t58;
t89 = qJD(4) * t54;
t88 = qJDD(3) * t55;
t87 = t55 * qJDD(1);
t85 = qJD(1) * qJD(2);
t83 = qJDD(1) * qJ(2);
t28 = t101 * t58;
t80 = t58 * t84;
t41 = pkin(3) * t55 + qJ(2);
t78 = t50 * t58;
t16 = -pkin(6) * t92 + t35 * t58;
t77 = g(1) * t59 + g(2) * t56;
t11 = -t54 * t91 - t55 * t89 + t57 * t78;
t76 = -t11 * t50 - t23 * t49;
t13 = qJD(3) * pkin(3) + t16;
t75 = -t13 * t54 - t98;
t27 = t101 * t55;
t74 = -t27 * t57 - t28 * t54;
t73 = t27 * t54 - t28 * t57;
t72 = -qJD(4) * t82 - t106 * t54;
t70 = t80 + t87;
t69 = 0.2e1 * qJ(2) * t84 + qJDD(3) * t60;
t68 = -t95 + t105;
t67 = -t77 + 0.2e1 * t85;
t66 = t67 + 0.2e1 * t83;
t65 = -t60 * t61 + t66;
t29 = t41 * qJD(1);
t53 = qJ(3) + qJ(4);
t46 = sin(t53);
t47 = cos(t53);
t34 = qJDD(1) * t60 + qJDD(2);
t25 = t58 * t34;
t7 = qJDD(3) * pkin(3) + pkin(6) * t106 - t35 * t91 + t25;
t64 = t15 * t89 + g(3) * t47 + (-t15 * t50 - t7) * t54 + t29 * t19 - t105 * t46;
t8 = -pkin(6) * t70 + t55 * t34 + t35 * t90;
t63 = g(3) * t46 + t75 * qJD(4) + t105 * t47 + t29 * t18 - t54 * t8 + t57 * t7;
t3 = -qJD(1) * t10 - t54 * t87 + t57 * t86;
t45 = qJDD(3) * t58;
t36 = pkin(3) * t90 + qJD(2);
t22 = qJD(3) * t28;
t21 = t101 * t91;
t14 = pkin(3) * t70 + t83 + t85;
t5 = t18 ^ 2 - t19 ^ 2;
t4 = (qJD(1) * t78 + t87) * t57 + t72;
t2 = -t18 * t50 + (-t50 * t92 - t87) * t57 - t72;
t1 = t19 * t50 + t3;
t6 = [qJDD(1), -t105, t77, -0.2e1 * t94 + t107, t66, t104 * pkin(1) + (t67 + t83) * qJ(2), qJDD(1) * t52 - 0.2e1 * t55 * t80, -0.2e1 * t55 * t86 + 0.2e1 * t84 * t97, -t55 * t61 + t45, -t58 * t61 - t88, 0, t55 * t65 + t58 * t69, -t55 * t69 + t58 * t65, t10 * t18 + t24 * t3, t10 * t19 + t11 * t18 - t23 * t3 - t24 * t4, t100, t76, 0, t36 * t19 + t41 * t4 + t14 * t23 + t29 * t11 + (-qJD(4) * t74 + t57 * t21 + t54 * t22) * t50 + t73 * t49 - t77 * t46, -t36 * t18 + t41 * t3 + t14 * t24 - t29 * t10 - (qJD(4) * t73 + t54 * t21 - t57 * t22) * t50 - t74 * t49 - t77 * t47; 0, 0, 0, qJDD(1), -t62, -t95 - t104, 0, 0, 0, 0, 0, t55 * t96 + t45, t58 * t96 - t88, 0, 0, 0, 0, 0, -qJD(1) * t19 + t100, qJD(1) * t18 + t76; 0, 0, 0, 0, 0, 0, t58 * t62 * t55, -t97 * t62, t86, -t87, qJDD(3), g(3) * t55 + t58 * t68 + t25, g(3) * t58 + (-t34 - t68) * t55, -t99, t5, t1, t2, t49, -(-t16 * t54 - t98) * t50 + (-t19 * t92 + t49 * t57 - t50 * t89) * pkin(3) + t63, (-qJD(4) * t13 + t16 * t50 - t8) * t57 + (-qJD(4) * t50 * t57 + t18 * t92 - t49 * t54) * pkin(3) + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t5, t1, t2, t49, -t50 * t75 + t63, (-t8 + (-qJD(4) + t50) * t13) * t57 + t64;];
tau_reg = t6;
