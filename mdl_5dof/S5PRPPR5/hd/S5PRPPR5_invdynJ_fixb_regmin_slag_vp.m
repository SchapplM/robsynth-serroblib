% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:39
% EndTime: 2019-12-31 17:38:40
% DurationCPUTime: 0.46s
% Computational Cost: add. (447->121), mult. (808->165), div. (0->0), fcn. (582->8), ass. (0->78)
t61 = sin(qJ(2));
t100 = qJD(1) * t61;
t56 = sin(pkin(8));
t58 = cos(pkin(8));
t63 = cos(qJ(2));
t97 = t63 * qJD(1);
t20 = -t58 * t100 + t56 * t97;
t110 = pkin(2) + pkin(3);
t50 = t63 * qJDD(1);
t99 = qJD(2) * t61;
t88 = qJD(1) * t99 + qJDD(3) - t50;
t12 = -t110 * qJDD(2) + t88;
t48 = t61 * qJDD(1);
t91 = qJDD(2) * qJ(3);
t17 = t91 + t48 + (qJD(3) + t97) * qJD(2);
t3 = t12 * t58 - t56 * t17;
t77 = t61 * t56 + t63 * t58;
t107 = t56 * t63;
t78 = -t58 * t61 + t107;
t57 = sin(pkin(7));
t59 = cos(pkin(7));
t81 = g(1) * t59 + g(2) * t57;
t115 = -g(3) * t77 - t81 * t78 - t3 + (qJD(3) * t56 - t20) * qJD(2);
t23 = t77 * qJD(2);
t4 = t56 * t12 + t58 * t17;
t114 = -g(3) * t78 + t81 * t77 - t4;
t30 = -qJ(3) * t56 - t110 * t58;
t27 = pkin(4) - t30;
t31 = t58 * qJ(3) - t110 * t56;
t28 = -pkin(6) + t31;
t65 = qJD(5) ^ 2;
t112 = -t28 * t65 + (pkin(4) + t27) * qJDD(2) + t115;
t111 = g(3) * t61 + t81 * t63 - t48;
t106 = t57 * t61;
t105 = t59 * t61;
t104 = t63 * pkin(2) + t61 * qJ(3);
t60 = sin(qJ(5));
t54 = t60 ^ 2;
t62 = cos(qJ(5));
t103 = -t62 ^ 2 + t54;
t66 = qJD(2) ^ 2;
t102 = t65 + t66;
t101 = qJ(3) * t63;
t98 = qJDD(2) * pkin(2);
t96 = qJDD(5) * t60;
t95 = qJDD(5) * t62;
t94 = t60 * qJDD(2);
t93 = t62 * qJDD(2);
t92 = qJD(2) * qJD(5);
t90 = -g(1) * t105 - g(2) * t106 + g(3) * t63;
t89 = 0.2e1 * t92;
t86 = t50 - t90;
t21 = t77 * qJD(1);
t84 = qJD(3) * t58 - t21;
t83 = t62 * t89;
t82 = qJD(3) - t97;
t26 = -t110 * qJD(2) + t82;
t37 = qJD(2) * qJ(3) + t100;
t7 = t26 * t58 - t37 * t56;
t8 = t56 * t26 + t58 * t37;
t80 = t56 * t7 - t58 * t8;
t79 = (-qJD(2) * pkin(2) + t82) * t61 + t37 * t63;
t22 = -qJD(2) * t107 + t58 * t99;
t75 = -qJD(2) * t22 + qJDD(2) * t77;
t74 = g(1) * t57 - g(2) * t59 + qJDD(4);
t18 = t88 - t98;
t71 = -t65 * t78 - t75;
t70 = -0.2e1 * t23 * qJD(5) + qJDD(5) * t78;
t5 = qJD(2) * pkin(4) - t7;
t69 = qJDD(2) * pkin(6) + qJD(2) * t5 + t114;
t68 = -qJDD(5) * t28 + (-qJD(2) * t27 - t5 - t84) * qJD(5);
t41 = t59 * t101;
t40 = t57 * t101;
t35 = qJDD(2) * t63 - t61 * t66;
t34 = qJDD(2) * t61 + t63 * t66;
t33 = -t60 * t65 + t95;
t32 = -t62 * t65 - t96;
t1 = [qJDD(1) - g(3), 0, t35, -t34, t35, t34, t79 * qJD(2) + t17 * t61 - t18 * t63 - g(3), t75, qJD(2) * t23 - qJDD(2) * t78, t22 * t7 + t23 * t8 - t3 * t77 - t4 * t78 - g(3), 0, 0, 0, 0, 0, t70 * t60 - t71 * t62, t71 * t60 + t70 * t62; 0, qJDD(2), t86, t111, -qJDD(3) + t86 + 0.2e1 * t98, 0.2e1 * qJD(2) * qJD(3) - t111 + 0.2e1 * t91, t17 * qJ(3) + t37 * qJD(3) - t18 * pkin(2) - g(1) * (-pkin(2) * t105 + t41) - g(2) * (-pkin(2) * t106 + t40) - g(3) * t104 - t79 * qJD(1), -qJDD(2) * t30 + t115, t84 * qJD(2) + qJDD(2) * t31 - t114, t4 * t31 + t3 * t30 - t8 * t21 + t7 * t20 - g(1) * t41 - g(2) * t40 - g(3) * (pkin(3) * t63 + t104) - t80 * qJD(3) + t81 * t61 * t110, qJDD(2) * t54 + t60 * t83, -0.2e1 * t103 * t92 + 0.2e1 * t60 * t93, t32, -t33, 0, t112 * t62 + t68 * t60, -t112 * t60 + t68 * t62; 0, 0, 0, 0, -qJDD(2), -t66, -qJD(2) * t37 + t18 + t90, -qJDD(2) * t58 - t56 * t66, qJDD(2) * t56 - t58 * t66, t80 * qJD(2) + t3 * t58 + t4 * t56 + t90, 0, 0, 0, 0, 0, (t60 * t89 - t93) * t58 + (-t102 * t62 - t96) * t56, (t83 + t94) * t58 + (t102 * t60 - t95) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t66 * t62, t103 * t66, -t94, -t93, qJDD(5), t69 * t60 + t74 * t62, -t74 * t60 + t69 * t62;];
tau_reg = t1;
