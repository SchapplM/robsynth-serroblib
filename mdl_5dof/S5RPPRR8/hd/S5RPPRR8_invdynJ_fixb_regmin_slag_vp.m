% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:16
% EndTime: 2019-12-31 18:01:17
% DurationCPUTime: 0.53s
% Computational Cost: add. (791->129), mult. (1166->165), div. (0->0), fcn. (714->10), ass. (0->87)
t123 = qJD(1) - qJD(4);
t57 = t123 ^ 2;
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t33 = t65 * t62 - t68 * t63;
t112 = t123 * t33;
t106 = t62 * qJ(2);
t70 = -pkin(1) - pkin(2);
t40 = t63 * t70 - t106;
t37 = -pkin(3) + t40;
t41 = t63 * qJ(2) + t62 * t70;
t111 = t65 * t37 + t68 * t41;
t103 = qJD(4) * t68;
t104 = qJD(4) * t65;
t45 = t70 * qJDD(1) + qJDD(2);
t39 = t63 * t45;
t93 = -pkin(3) - t106;
t98 = qJD(1) * qJD(2);
t95 = t62 * t98;
t14 = t93 * qJDD(1) + t39 - t95;
t99 = qJ(2) * qJDD(1);
t110 = t62 * t45 + t63 * t99;
t94 = t63 * t98;
t19 = t94 + t110;
t49 = t70 * qJD(1) + qJD(2);
t42 = t63 * t49;
t21 = t93 * qJD(1) + t42;
t100 = qJ(2) * qJD(1);
t24 = t63 * t100 + t62 * t49;
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t97 = pkin(8) + qJ(4);
t91 = sin(t97);
t92 = cos(t97);
t25 = -t66 * t91 - t69 * t92;
t26 = -t66 * t92 + t69 * t91;
t75 = g(1) * t25 + g(2) * t26 + t21 * t103 - t24 * t104 + t65 * t14 + t68 * t19;
t76 = -g(1) * t26 + g(2) * t25 + t24 * t103 + t21 * t104 - t68 * t14 + t65 * t19;
t35 = t68 * t62 + t65 * t63;
t71 = qJD(5) ^ 2;
t58 = qJDD(1) - qJDD(4);
t82 = t33 * t58 - t57 * t35;
t122 = -t35 * t71 + t82;
t121 = t58 * pkin(4);
t120 = t123 * pkin(4);
t119 = (t35 * qJD(2) + t111 * qJD(4)) * t123;
t118 = (t65 * t21 + t68 * t24) * t123;
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t115 = t64 * t67;
t114 = t67 * t58;
t109 = t69 * pkin(1) + t66 * qJ(2);
t108 = g(1) * t66 - g(2) * t69;
t60 = t64 ^ 2;
t107 = -t67 ^ 2 + t60;
t105 = pkin(1) * qJDD(1);
t102 = qJD(5) * t123;
t101 = qJDD(3) + g(3);
t96 = 0.2e1 * t98;
t90 = qJDD(2) - t105;
t87 = g(1) * t69 + g(2) * t66;
t7 = t68 * t21 - t65 * t24;
t86 = (-t62 * t100 + t42) * t62 - t24 * t63;
t85 = t68 * t37 - t65 * t41;
t83 = -t121 - t76;
t3 = -t7 + t120;
t80 = t58 * pkin(7) + t123 * t3 - t75;
t79 = -pkin(7) * qJDD(5) + (t3 + t7 + t120) * qJD(5);
t10 = -pkin(7) + t111;
t5 = -t33 * qJD(2) + t85 * qJD(4);
t9 = pkin(4) - t85;
t78 = -qJDD(5) * t10 + (-t123 * t9 - t3 - t5) * qJD(5);
t77 = -0.2e1 * t112 * qJD(5) - qJDD(5) * t35;
t74 = pkin(7) * t71 + t118 + t121 - t83;
t73 = t10 * t71 - t58 * t9 - t119 + t83;
t72 = qJD(1) ^ 2;
t53 = t69 * qJ(2);
t44 = qJDD(5) * t67 - t71 * t64;
t43 = qJDD(5) * t64 + t71 * t67;
t36 = t66 * t62 + t69 * t63;
t34 = t69 * t62 - t66 * t63;
t22 = -0.2e1 * t102 * t115 - t60 * t58;
t18 = t39 + (-t98 - t99) * t62;
t15 = t107 * t102 - t64 * t114;
t1 = [qJDD(1), t108, t87, -qJDD(2) + 0.2e1 * t105 + t108, -t87 + t96 + 0.2e1 * t99, -t90 * pkin(1) - g(1) * (-t66 * pkin(1) + t53) - g(2) * t109 + (t96 + t99) * qJ(2), 0.2e1 * t95 - g(1) * t34 - g(2) * t36 - t39 + (-t40 + t106) * qJDD(1), -g(1) * t36 + g(2) * t34 + t41 * qJDD(1) + t110 + 0.2e1 * t94, t19 * t41 + t18 * t40 - g(1) * (t70 * t66 + t53) - g(2) * (t69 * pkin(2) + t109) - t86 * qJD(2), t58, -t85 * t58 + t119 + t76, t111 * t58 + t123 * t5 + t75, -t22, -0.2e1 * t15, -t43, -t44, 0, t78 * t64 - t73 * t67, t73 * t64 + t78 * t67; 0, 0, 0, -qJDD(1), -t72, -t72 * qJ(2) - t108 + t90, -t63 * qJDD(1) - t62 * t72, t62 * qJDD(1) - t63 * t72, t86 * qJD(1) + t18 * t63 + t19 * t62 - t108, 0, t82, t112 * t123 + t35 * t58, 0, 0, 0, 0, 0, t122 * t67 + t77 * t64, -t122 * t64 + t77 * t67; 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t76 - t118, -t123 * t7 - t75, t22, 0.2e1 * t15, t43, t44, 0, t79 * t64 - t74 * t67, t74 * t64 + t79 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t115, t107 * t57, -t64 * t58, -t114, qJDD(5), t101 * t67 + t80 * t64, -t101 * t64 + t80 * t67;];
tau_reg = t1;
