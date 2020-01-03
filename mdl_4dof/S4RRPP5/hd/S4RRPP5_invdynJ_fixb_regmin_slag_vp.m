% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:38
% EndTime: 2019-12-31 17:00:40
% DurationCPUTime: 0.80s
% Computational Cost: add. (491->179), mult. (1048->204), div. (0->0), fcn. (518->4), ass. (0->107)
t66 = sin(qJ(2));
t53 = t66 * qJ(3);
t68 = cos(qJ(2));
t113 = t68 * pkin(2) + t53;
t132 = -pkin(1) - t113;
t15 = t132 * qJD(1);
t133 = qJDD(1) * t132;
t126 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t106 = qJD(1) * t68;
t47 = pkin(5) * t106;
t20 = pkin(3) * t106 + t47;
t131 = -qJD(4) - t20;
t115 = pkin(2) + qJ(4);
t130 = t115 * qJD(2);
t129 = t115 * qJDD(2);
t51 = t66 * qJDD(1);
t101 = qJD(1) * qJD(2);
t95 = t68 * t101;
t128 = (-t51 - t95) * pkin(3);
t127 = 0.2e1 * t126;
t102 = qJD(2) * qJ(3);
t14 = t102 - t131;
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t88 = g(1) * t69 + g(2) * t67;
t117 = t66 * t69;
t118 = t66 * t67;
t125 = -g(1) * t117 - g(2) * t118 + g(3) * t68;
t124 = pkin(3) + pkin(5);
t123 = g(1) * t67;
t120 = g(2) * t69;
t119 = g(3) * t66;
t71 = qJD(1) ^ 2;
t116 = t66 * t71;
t52 = t68 * qJDD(1);
t114 = qJ(3) * t52 + qJD(3) * t106;
t63 = t66 ^ 2;
t64 = t68 ^ 2;
t112 = t63 - t64;
t111 = qJ(3) * t68;
t110 = qJD(2) * pkin(2);
t109 = t68 * qJ(4);
t108 = pkin(5) * qJDD(2);
t107 = qJD(1) * t66;
t27 = t124 * t68;
t21 = qJD(2) * t27;
t105 = qJDD(2) * pkin(2);
t104 = t66 * qJD(3);
t46 = pkin(5) * t107;
t18 = -pkin(3) * t107 - t46;
t103 = qJD(3) - t18;
t100 = t68 * t116;
t32 = pkin(5) * t95;
t42 = pkin(5) * t51;
t99 = qJDD(3) + t32 + t42;
t44 = pkin(5) * t52;
t98 = pkin(3) * t52 + qJDD(4) + t44;
t97 = t124 * qJD(2);
t96 = t66 * t101;
t90 = t109 + t113;
t16 = -pkin(1) - t90;
t80 = -t115 * t68 - pkin(1) - t53;
t7 = t80 * qJD(1);
t94 = qJD(1) * t16 + t7;
t93 = t67 * pkin(5) - t132 * t69;
t92 = t42 + t125;
t19 = t66 * t97;
t70 = qJD(2) ^ 2;
t89 = pkin(5) * t70 + t120;
t24 = -t47 - t102;
t86 = (qJD(3) + t46 - t110) * t68 + t24 * t66;
t85 = -qJDD(3) - t92;
t84 = qJ(4) * t66 - t111;
t11 = t99 - t105;
t82 = -0.2e1 * pkin(1) * t101 - t108;
t81 = -t68 * t102 - t104;
t79 = 0.2e1 * qJDD(1) * pkin(1) - t89;
t33 = pkin(2) * t96;
t72 = t84 * qJD(2) - t68 * qJD(4) - t104;
t1 = t72 * qJD(1) + t80 * qJDD(1) + t33;
t49 = t66 * t110;
t5 = t49 + t72;
t78 = -qJD(1) * t5 - qJDD(1) * t16 - t1 - t120;
t77 = -0.2e1 * t15 * qJD(2) + t108;
t76 = -t128 + t99 - t129;
t75 = t7 * t106 - t88 * t68 + t98;
t13 = t49 + t81;
t3 = t81 * qJD(1) + t133 + t33;
t74 = qJD(1) * t13 + t133 + t3 + t89;
t8 = pkin(5) * t96 - t126 - t44;
t73 = t86 * qJD(2) + t11 * t66 - t8 * t68;
t58 = t69 * pkin(5);
t50 = pkin(2) * t107;
t38 = t68 * t123;
t37 = g(1) * t118;
t31 = t69 * t111;
t29 = t67 * t111;
t28 = -t63 * t71 - t70;
t26 = t124 * t66;
t25 = qJDD(2) + t100;
t17 = -qJ(3) * t106 + t50;
t12 = t84 * qJD(1) + t50;
t10 = t103 - t130;
t9 = t15 * t107;
t4 = -qJD(1) * t19 + t126 + t98;
t2 = -qJD(2) * qJD(4) + t76;
t6 = [qJDD(1), -t120 + t123, t88, t63 * qJDD(1) + 0.2e1 * t66 * t95, -0.2e1 * t112 * t101 + 0.2e1 * t66 * t52, qJDD(2) * t66 + t70 * t68, qJDD(2) * t68 - t70 * t66, 0, t82 * t66 + t79 * t68 + t38, -t79 * t66 + t82 * t68 - t37, (t63 + t64) * qJDD(1) * pkin(5) + t73 - t88, t77 * t66 + t74 * t68 - t38, -t74 * t66 + t77 * t68 + t37, t73 * pkin(5) - g(1) * t58 - g(2) * t93 + t15 * t13 + (-t123 + t3) * t132, (t10 * qJD(2) + qJDD(1) * t27 + t4 + (qJD(2) * t26 - t19) * qJD(1)) * t68 + (-t14 * qJD(2) + qJDD(1) * t26 + t2) * t66 - t88, t27 * qJDD(2) + t37 + (-t94 * t68 - t19) * qJD(2) + t78 * t66, -t26 * qJDD(2) + t38 + (t94 * t66 - t21) * qJD(2) + t78 * t68, t1 * t16 + t7 * t5 + t2 * t26 + t10 * t21 + t4 * t27 - t14 * t19 - g(1) * (t69 * pkin(3) + t58) - g(2) * (t69 * t109 + t93) + (-g(1) * (t132 - t109) - g(2) * pkin(3)) * t67; 0, 0, 0, -t100, t112 * t71, t51, t52, qJDD(2), pkin(1) * t116 - t92, t119 - t44 + (pkin(1) * t71 + t88) * t68, -pkin(2) * t51 + (-qJD(2) * t113 - t86) * qJD(1) + t114, -t17 * t106 - 0.2e1 * t105 - t85 + t9, t44 + (qJD(1) * t17 - g(3)) * t66 + (qJD(1) * t15 - t88) * t68 + t127, -t8 * qJ(3) - t24 * qJD(3) - t11 * pkin(2) - t15 * t17 - g(1) * (-pkin(2) * t117 + t31) - g(2) * (-pkin(2) * t118 + t29) - g(3) * t113 - t86 * qJD(1) * pkin(5), -t115 * t51 + (-t10 - t18 - t130) * t106 + t114, -t18 * qJD(2) + (-g(3) + (t12 - t97) * qJD(1)) * t66 + t75 + t127, -t32 + (0.2e1 * qJD(4) + t20) * qJD(2) + (t12 * t68 - t66 * t7) * qJD(1) + 0.2e1 * t129 + t85 + t128, -g(1) * t31 - g(2) * t29 - g(3) * t90 + t4 * qJ(3) + t131 * t10 + t103 * t14 - t7 * t12 + (t66 * t88 - t2) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t25, t28, t24 * qJD(2) + t11 + t125 + t9, t51, t28, -t25, t7 * t107 + (-qJD(4) - t14) * qJD(2) + t76 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, qJDD(2) - t100, -t64 * t71 - t70, -t119 + (-t124 * t107 + t10) * qJD(2) + t75 + t126;];
tau_reg = t6;
