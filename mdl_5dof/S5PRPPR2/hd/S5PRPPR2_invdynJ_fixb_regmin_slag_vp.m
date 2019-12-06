% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:44
% EndTime: 2019-12-05 15:24:47
% DurationCPUTime: 0.75s
% Computational Cost: add. (603->148), mult. (1261->208), div. (0->0), fcn. (1010->14), ass. (0->100)
t104 = qJD(1) * qJD(2);
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t130 = qJDD(1) * t80 + t82 * t104;
t67 = t82 * qJDD(1);
t37 = qJDD(2) * pkin(2) - t80 * t104 + t67;
t74 = sin(pkin(8));
t77 = cos(pkin(8));
t11 = -t130 * t74 + t37 * t77;
t88 = qJDD(4) - t11;
t8 = -qJDD(2) * pkin(3) + t88;
t72 = qJ(2) + pkin(8);
t64 = sin(t72);
t66 = cos(t72);
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t96 = g(1) * t78 + g(2) * t75;
t86 = -g(3) * t66 + t96 * t64;
t132 = t8 - t86;
t114 = qJD(1) * t82;
t115 = qJD(1) * t80;
t54 = t77 * t115;
t27 = t74 * t114 + t54;
t126 = pkin(2) * t77;
t58 = -pkin(3) - t126;
t131 = qJD(2) * t27 - qJDD(2) * t58 - t132;
t73 = sin(pkin(9));
t76 = cos(pkin(9));
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t41 = t73 * t81 + t76 * t79;
t32 = t41 * qJD(2);
t34 = t41 * qJD(5);
t53 = t74 * t115;
t31 = t77 * t114 - t53;
t112 = qJD(4) - t31;
t116 = t73 ^ 2 + t76 ^ 2;
t129 = qJD(2) * t116;
t127 = qJD(5) ^ 2;
t123 = g(3) * t64;
t56 = pkin(2) * t74 + qJ(4);
t121 = pkin(6) + t56;
t12 = t130 * t77 + t74 * t37;
t7 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t12;
t4 = t73 * qJDD(3) + t76 * t7;
t120 = t66 * t75;
t119 = t66 * t78;
t118 = t73 * t79;
t117 = t81 * t76;
t47 = qJD(2) * pkin(2) + t114;
t24 = t74 * t47 + t54;
t111 = qJDD(1) - g(3);
t39 = -t117 + t118;
t108 = qJDD(5) * t39;
t107 = qJDD(5) * t41;
t106 = t73 * qJDD(2);
t105 = t76 * qJDD(2);
t101 = qJD(2) * t117;
t103 = qJD(5) * t101 + t79 * t105 + t81 * t106;
t102 = qJD(2) * t118;
t100 = -pkin(4) * t76 - pkin(3);
t99 = -g(1) * t75 + g(2) * t78;
t23 = t47 * t77 - t53;
t97 = qJDD(2) * t116;
t60 = t76 * qJDD(3);
t3 = -t7 * t73 + t60;
t95 = -t3 * t73 + t4 * t76;
t94 = qJD(4) - t23;
t93 = -t81 * t105 + t79 * t106;
t22 = qJD(2) * qJ(4) + t24;
t17 = t76 * qJD(3) - t22 * t73;
t18 = t73 * qJD(3) + t76 * t22;
t92 = t17 * t73 - t18 * t76;
t35 = t121 * t73;
t36 = t121 * t76;
t91 = -t35 * t81 - t36 * t79;
t90 = -t35 * t79 + t36 * t81;
t40 = t74 * t82 + t77 * t80;
t38 = t74 * t80 - t77 * t82;
t26 = t40 * qJD(2);
t89 = qJD(2) * t26 + qJDD(2) * t38;
t33 = t39 * qJD(5);
t84 = -g(3) * t82 + t96 * t80;
t83 = qJD(2) ^ 2;
t71 = pkin(9) + qJ(5);
t65 = cos(t71);
t63 = sin(t71);
t43 = t100 - t126;
t30 = t38 * qJD(2);
t28 = -t101 + t102;
t21 = -qJD(2) * pkin(3) + t94;
t19 = t100 * qJD(2) + t94;
t16 = -qJD(5) * t34 - t108;
t15 = -qJD(5) * t33 + t107;
t14 = qJD(2) * t34 + t93;
t13 = -qJD(5) * t102 + t103;
t5 = t100 * qJDD(2) + t88;
t2 = pkin(6) * t105 + t4;
t1 = t60 + (-pkin(6) * qJDD(2) - t7) * t73;
t6 = [t111, 0, qJDD(2) * t82 - t80 * t83, -qJDD(2) * t80 - t82 * t83, -t11 * t38 + t12 * t40 - t23 * t26 - t24 * t30 - g(3), -t89 * t76, t89 * t73, -t30 * t129 + t40 * t97, t21 * t26 + t92 * t30 + t38 * t8 + t95 * t40 - g(3), 0, 0, 0, 0, 0, t38 * t14 + t26 * t28 + t30 * t34 + (t39 * t127 - t107) * t40, t38 * t13 + t26 * t32 - t30 * t33 + (t41 * t127 + t108) * t40; 0, qJDD(2), t67 + t84, -t111 * t80 + t96 * t82, t23 * t27 - t24 * t31 + (t11 * t77 + t12 * t74 + t84) * pkin(2), t131 * t76, -t131 * t73, t112 * t129 + t56 * t97 - t96 * t66 - t123 + t95, t8 * t58 - t21 * t27 - g(3) * (pkin(2) * t82 + pkin(3) * t66 + qJ(4) * t64) + (t112 * t18 + t4 * t56) * t76 + (-t112 * t17 - t3 * t56) * t73 + t96 * (pkin(2) * t80 + pkin(3) * t64 - qJ(4) * t66), t13 * t41 - t32 * t33, -t13 * t39 - t14 * t41 + t28 * t33 - t32 * t34, t15, t16, 0, t91 * qJDD(5) + t43 * t14 + t5 * t39 + t19 * t34 - t27 * t28 + t86 * t65 + (-t90 * qJD(5) - t112 * t41) * qJD(5), -t90 * qJDD(5) + t43 * t13 + t5 * t41 - t19 * t33 - t27 * t32 - t86 * t63 + (-t91 * qJD(5) + t112 * t39) * qJD(5); 0, 0, 0, 0, qJDD(3) + t99, 0, 0, 0, t3 * t76 + t4 * t73 + t99, 0, 0, 0, 0, 0, t16, -t15; 0, 0, 0, 0, 0, -t105, t106, -t116 * t83, t92 * qJD(2) + t132, 0, 0, 0, 0, 0, 0.2e1 * qJD(5) * t32 + t93, (-t28 - t102) * qJD(5) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t28, -t28 ^ 2 + t32 ^ 2, (t28 - t102) * qJD(5) + t103, -t93, qJDD(5), -t79 * t2 + t81 * t1 - t19 * t32 - g(1) * (-t63 * t119 + t65 * t75) - g(2) * (-t63 * t120 - t65 * t78) + t63 * t123, -t81 * t2 - t79 * t1 + t19 * t28 - g(1) * (-t65 * t119 - t63 * t75) - g(2) * (-t65 * t120 + t63 * t78) + t65 * t123;];
tau_reg = t6;
