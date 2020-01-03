% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR7
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
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:06
% EndTime: 2019-12-31 16:54:08
% DurationCPUTime: 0.77s
% Computational Cost: add. (973->180), mult. (2390->243), div. (0->0), fcn. (1824->10), ass. (0->106)
t67 = cos(pkin(7));
t72 = cos(qJ(3));
t110 = t72 * t67;
t53 = qJD(1) * t110;
t66 = sin(pkin(7));
t69 = sin(qJ(3));
t115 = t66 * t69;
t99 = qJD(1) * t115;
t37 = t53 - t99;
t32 = qJD(4) - t37;
t42 = t66 * t72 + t67 * t69;
t38 = t42 * qJD(1);
t109 = pkin(5) + qJ(2);
t47 = t109 * t66;
t43 = qJD(1) * t47;
t48 = t109 * t67;
t44 = qJD(1) * t48;
t24 = -t43 * t69 + t44 * t72;
t101 = qJD(1) * qJD(2);
t127 = t109 * qJDD(1) + t101;
t30 = t127 * t66;
t31 = t127 * t67;
t87 = t30 * t72 + t31 * t69;
t5 = -qJDD(3) * pkin(3) + qJD(3) * t24 + t87;
t65 = pkin(7) + qJ(3);
t58 = sin(t65);
t59 = cos(t65);
t70 = sin(qJ(1));
t73 = cos(qJ(1));
t90 = g(1) * t73 + g(2) * t70;
t78 = -g(3) * t59 + t90 * t58;
t130 = -(pkin(3) * t38 + t32 * pkin(6)) * t32 - t5 + t78;
t129 = qJ(2) * qJDD(1);
t102 = t67 * qJDD(1);
t103 = t66 * qJDD(1);
t100 = qJD(3) * t53 + t69 * t102 + t72 * t103;
t21 = -qJD(3) * t99 + t100;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t29 = qJD(3) * t68 + t38 * t71;
t8 = qJD(4) * t29 - t71 * qJDD(3) + t68 * t21;
t41 = -t110 + t115;
t85 = -t47 * t72 - t48 * t69;
t10 = -t41 * qJD(2) + t85 * qJD(3);
t40 = t42 * qJD(3);
t88 = -t72 * t102 + t69 * t103;
t22 = qJD(1) * t40 + t88;
t15 = qJDD(4) + t22;
t23 = -t43 * t72 - t44 * t69;
t16 = -qJD(3) * pkin(3) - t23;
t55 = -pkin(2) * t67 - pkin(1);
t20 = pkin(3) * t41 - pkin(6) * t42 + t55;
t26 = -t47 * t69 + t48 * t72;
t39 = t41 * qJD(3);
t86 = -t30 * t69 + t31 * t72;
t46 = t55 * qJD(1) + qJD(2);
t9 = -pkin(3) * t37 - pkin(6) * t38 + t46;
t97 = qJDD(3) * pkin(6) + qJD(3) * t23 + qJD(4) * t9 + t86;
t126 = -t26 * t15 - t16 * t39 - (qJD(4) * t20 + t10) * t32 - t97 * t41 + t5 * t42;
t125 = g(3) * t58;
t104 = t71 * qJD(3);
t107 = qJD(4) * t68;
t7 = qJD(4) * t104 + t68 * qJDD(3) - t38 * t107 + t71 * t21;
t123 = t7 * t68;
t122 = t15 * t68;
t121 = t16 * t42;
t120 = t20 * t15;
t27 = t38 * t68 - t104;
t119 = t27 * t32;
t118 = t27 * t38;
t117 = t29 * t32;
t116 = t29 * t38;
t114 = t68 * t70;
t113 = t68 * t73;
t112 = t70 * t71;
t12 = t71 * t15;
t111 = t71 * t73;
t108 = t66 ^ 2 + t67 ^ 2;
t106 = qJD(4) * t71;
t105 = qJDD(1) * pkin(1);
t98 = t108 * qJD(1) ^ 2;
t17 = qJD(3) * pkin(6) + t24;
t45 = t55 * qJDD(1) + qJDD(2);
t6 = pkin(3) * t22 - pkin(6) * t21 + t45;
t94 = qJD(4) * t17 - t6;
t93 = t32 * t71;
t91 = 0.2e1 * t108;
t89 = g(1) * t70 - g(2) * t73;
t84 = t12 + (t37 * t68 - t107) * t32;
t83 = -t97 + t125;
t82 = -t42 * t107 - t39 * t71;
t80 = -t89 - t105;
t57 = qJDD(2) - t105;
t79 = -t57 - t80;
t77 = -pkin(6) * t15 + (t16 + t23) * t32;
t76 = t91 * t101 - t90;
t36 = t59 * t111 + t114;
t35 = -t59 * t113 + t112;
t34 = -t59 * t112 + t113;
t33 = t59 * t114 + t111;
t19 = pkin(3) * t40 + pkin(6) * t39;
t11 = t42 * qJD(2) + t26 * qJD(3);
t3 = t71 * t6;
t2 = t17 * t71 + t68 * t9;
t1 = -t17 * t68 + t71 * t9;
t4 = [qJDD(1), t89, t90, t79 * t67, -t79 * t66, t91 * t129 + t76, (-t57 + t89) * pkin(1) + (t108 * t129 + t76) * qJ(2), t21 * t42 - t38 * t39, -t21 * t41 - t22 * t42 - t37 * t39 - t38 * t40, -qJD(3) * t39 + qJDD(3) * t42, -qJD(3) * t40 - qJDD(3) * t41, 0, -qJD(3) * t11 + qJDD(3) * t85 + t22 * t55 + t40 * t46 + t41 * t45 + t89 * t59, -qJD(3) * t10 - qJDD(3) * t26 + t21 * t55 - t39 * t46 + t42 * t45 - t89 * t58, t7 * t71 * t42 + t82 * t29, -(-t27 * t71 - t29 * t68) * t39 + (-t123 - t71 * t8 + (t27 * t68 - t29 * t71) * qJD(4)) * t42, t42 * t12 + t29 * t40 + t82 * t32 + t41 * t7, -t42 * t122 - t27 * t40 - t41 * t8 + (-t42 * t106 + t39 * t68) * t32, t15 * t41 + t32 * t40, -g(1) * t34 - g(2) * t36 + t1 * t40 + t11 * t27 - t85 * t8 + t3 * t41 + (t120 + t19 * t32 + (-t17 * t41 - t26 * t32 + t121) * qJD(4)) * t71 + t126 * t68, -g(1) * t33 - g(2) * t35 + t11 * t29 - t2 * t40 - t85 * t7 + (-(-qJD(4) * t26 + t19) * t32 - t120 + t94 * t41 - qJD(4) * t121) * t68 + t126 * t71; 0, 0, 0, -t102, t103, -t98, -qJ(2) * t98 + qJDD(2) + t80, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t38 + t88, (t37 - t99) * qJD(3) + t100, 0, 0, 0, 0, 0, t84 - t118, -t32 ^ 2 * t71 - t116 - t122; 0, 0, 0, 0, 0, 0, 0, -t38 * t37, -t37 ^ 2 + t38 ^ 2, (-t37 - t99) * qJD(3) + t100, -t88, qJDD(3), -t38 * t46 + t78 - t87, -t37 * t46 + t90 * t59 + t125 - t86, t29 * t93 + t123, (t7 - t119) * t71 + (-t8 - t117) * t68, t32 * t93 - t116 + t122, t84 + t118, -t32 * t38, -pkin(3) * t8 - t1 * t38 + t130 * t71 - t24 * t27 + t77 * t68, -pkin(3) * t7 - t130 * t68 + t2 * t38 - t24 * t29 + t77 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t7 + t119, t117 - t8, t15, -g(1) * t35 + g(2) * t33 - t17 * t106 - t16 * t29 + t2 * t32 + t83 * t68 + t3, g(1) * t36 - g(2) * t34 + t1 * t32 + t16 * t27 + t94 * t68 + t83 * t71;];
tau_reg = t4;
