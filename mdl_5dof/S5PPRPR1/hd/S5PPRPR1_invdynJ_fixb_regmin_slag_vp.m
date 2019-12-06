% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:22
% EndTime: 2019-12-05 15:01:26
% DurationCPUTime: 0.66s
% Computational Cost: add. (524->130), mult. (1161->177), div. (0->0), fcn. (968->14), ass. (0->84)
t61 = sin(pkin(8));
t64 = cos(pkin(8));
t67 = sin(qJ(3));
t69 = cos(qJ(3));
t117 = -t61 * t67 + t64 * t69;
t118 = t117 * qJD(1);
t33 = t61 * t69 + t64 * t67;
t29 = t33 * qJD(3);
t60 = sin(pkin(9));
t63 = cos(pkin(9));
t66 = sin(qJ(5));
t68 = cos(qJ(5));
t32 = t60 * t68 + t63 * t66;
t24 = t32 * qJD(3);
t27 = t32 * qJD(5);
t77 = qJD(4) - t118;
t62 = sin(pkin(7));
t65 = cos(pkin(7));
t85 = g(1) * t65 + g(2) * t62;
t115 = t33 * qJDD(1);
t59 = pkin(8) + qJ(3);
t52 = sin(t59);
t54 = cos(t59);
t74 = -g(3) * t54 + t85 * t52;
t110 = g(3) * t52;
t114 = t85 * t54 + t110;
t113 = qJD(5) ^ 2;
t7 = qJDD(3) * qJ(4) + t115 + (qJD(4) + t118) * qJD(3);
t4 = t60 * qJDD(2) + t63 * t7;
t108 = t54 * t62;
t107 = t54 * t65;
t106 = t60 * t66;
t101 = t68 * t63;
t100 = pkin(6) + qJ(4);
t99 = t60 ^ 2 + t63 ^ 2;
t25 = t33 * qJD(1);
t98 = qJD(3) * t25;
t28 = t117 * qJD(3);
t97 = qJD(3) * t28;
t96 = qJDD(3) * pkin(3);
t30 = -t101 + t106;
t95 = qJDD(5) * t30;
t94 = qJDD(5) * t32;
t93 = t60 * qJDD(3);
t92 = t63 * qJDD(3);
t88 = qJD(3) * t101;
t90 = qJD(5) * t88 + t66 * t92 + t68 * t93;
t89 = qJD(3) * t106;
t46 = -pkin(4) * t63 - pkin(3);
t87 = -g(1) * t62 + g(2) * t65;
t86 = t99 * qJDD(3);
t48 = t63 * qJDD(2);
t3 = -t60 * t7 + t48;
t84 = -t3 * t60 + t4 * t63;
t83 = t66 * t93 - t68 * t92;
t20 = qJD(3) * qJ(4) + t25;
t15 = t63 * qJD(2) - t20 * t60;
t16 = t60 * qJD(2) + t63 * t20;
t81 = t15 * t60 - t16 * t63;
t34 = t100 * t60;
t35 = t100 * t63;
t80 = -t34 * t68 - t35 * t66;
t79 = -t34 * t66 + t35 * t68;
t78 = qJD(3) * t29 - qJDD(3) * t117;
t76 = -qJD(1) * t29 + t117 * qJDD(1);
t26 = t30 * qJD(5);
t75 = qJDD(4) - t76;
t72 = -t74 - t96;
t71 = t74 + t98;
t58 = pkin(9) + qJ(5);
t53 = cos(t58);
t51 = sin(t58);
t21 = -t88 + t89;
t19 = -qJD(3) * pkin(3) + t77;
t17 = t46 * qJD(3) + t77;
t14 = -qJD(5) * t27 - t95;
t13 = -qJD(5) * t26 + t94;
t12 = qJD(3) * t27 + t83;
t11 = -qJD(5) * t89 + t90;
t8 = t75 - t96;
t5 = t46 * qJDD(3) + t75;
t2 = pkin(6) * t92 + t4;
t1 = t48 + (-pkin(6) * qJDD(3) - t7) * t60;
t6 = [qJDD(1) - g(3), -g(3) + (t61 ^ 2 + t64 ^ 2) * qJDD(1), 0, -t78, -qJDD(3) * t33 - t97, -t78 * t63, t78 * t60, t33 * t86 + t99 * t97, -t117 * t8 + t19 * t29 - t81 * t28 + t84 * t33 - g(3), 0, 0, 0, 0, 0, -t117 * t12 + t29 * t21 - t28 * t27 + (t30 * t113 - t94) * t33, -t117 * t11 + t29 * t24 + t28 * t26 + (t32 * t113 + t95) * t33; 0, qJDD(2) + t87, 0, 0, 0, 0, 0, 0, t3 * t63 + t4 * t60 + t87, 0, 0, 0, 0, 0, t14, -t13; 0, 0, qJDD(3), t71 + t76, -t115 + t114, (t71 - t8 + t96) * t63, (t72 + t8 - t98) * t60, t77 * qJD(3) * t99 + qJ(4) * t86 - t114 + t84, -t8 * pkin(3) - t19 * t25 - g(3) * (pkin(3) * t54 + qJ(4) * t52) + (t4 * qJ(4) + t77 * t16) * t63 + (-t3 * qJ(4) - t77 * t15) * t60 + t85 * (pkin(3) * t52 - qJ(4) * t54), t11 * t32 - t24 * t26, -t11 * t30 - t12 * t32 + t21 * t26 - t24 * t27, t13, t14, 0, t80 * qJDD(5) + t46 * t12 + t5 * t30 + t17 * t27 - t25 * t21 + t74 * t53 + (-t79 * qJD(5) - t77 * t32) * qJD(5), -t79 * qJDD(5) + t46 * t11 + t5 * t32 - t17 * t26 - t25 * t24 - t74 * t51 + (-t80 * qJD(5) + t77 * t30) * qJD(5); 0, 0, 0, 0, 0, -t92, t93, -t99 * qJD(3) ^ 2, t81 * qJD(3) + t72 + t75, 0, 0, 0, 0, 0, 0.2e1 * qJD(5) * t24 + t83, (-t21 - t89) * qJD(5) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t21, -t21 ^ 2 + t24 ^ 2, (t21 - t89) * qJD(5) + t90, -t83, qJDD(5), -t66 * t2 + t68 * t1 - t17 * t24 - g(1) * (-t51 * t107 + t53 * t62) - g(2) * (-t51 * t108 - t53 * t65) + t51 * t110, -t68 * t2 - t66 * t1 + t17 * t21 - g(1) * (-t53 * t107 - t51 * t62) - g(2) * (-t53 * t108 + t51 * t65) + t53 * t110;];
tau_reg = t6;
