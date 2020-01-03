% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:28
% EndTime: 2019-12-31 17:32:29
% DurationCPUTime: 0.65s
% Computational Cost: add. (857->155), mult. (1848->209), div. (0->0), fcn. (1432->10), ass. (0->96)
t72 = cos(qJ(3));
t102 = t72 * qJD(2);
t91 = qJD(4) - t102;
t66 = sin(pkin(8));
t69 = sin(qJ(5));
t114 = t69 * t66;
t67 = cos(pkin(8));
t71 = cos(qJ(5));
t39 = -t71 * t67 + t114;
t70 = sin(qJ(3));
t26 = t39 * t70;
t41 = t71 * t66 + t69 * t67;
t33 = t41 * qJD(3);
t104 = t67 * qJD(1);
t103 = t70 * qJD(2);
t47 = qJD(3) * qJ(4) + t103;
t28 = -t66 * qJD(1) + t67 * t47;
t82 = (-t66 * t47 - t104) * t66 - t28 * t67;
t120 = t82 * t72 - (-qJD(3) * pkin(3) + t91) * t70;
t119 = t33 ^ 2;
t111 = pkin(6) + qJ(4);
t43 = t111 * t66;
t44 = t111 * t67;
t16 = -t71 * t43 - t69 * t44;
t76 = t39 * t72;
t118 = qJD(2) * t76 - t39 * qJD(4) + t16 * qJD(5);
t17 = -t69 * t43 + t71 * t44;
t117 = -t17 * qJD(5) - t91 * t41;
t107 = qJD(3) * t67;
t94 = t71 * t107;
t95 = qJD(3) * t114;
t31 = -t94 + t95;
t116 = t33 * t31;
t73 = qJD(3) ^ 2;
t112 = t73 * t72;
t62 = t66 ^ 2;
t63 = t67 ^ 2;
t110 = t62 + t63;
t109 = cos(pkin(7));
t108 = sin(pkin(7));
t106 = qJD(3) * t70;
t105 = qJDD(3) * pkin(3);
t101 = t66 * qJDD(3);
t100 = t67 * qJDD(1);
t99 = t67 * qJDD(3);
t98 = t70 * qJDD(2);
t97 = t72 * qJDD(2);
t96 = qJD(5) * t94 + t71 * t101 + t69 * t99;
t57 = t67 * pkin(4) + pkin(3);
t24 = qJDD(3) * qJ(4) + t98 + (qJD(4) + t102) * qJD(3);
t12 = -t100 + (-pkin(6) * qJDD(3) - t24) * t66;
t15 = -t66 * qJDD(1) + t67 * t24;
t13 = pkin(6) * t99 + t15;
t93 = t71 * t12 - t69 * t13;
t58 = qJD(3) * t103;
t92 = t110 * qJDD(3);
t38 = -t108 * t70 - t109 * t72;
t40 = -t108 * t72 + t109 * t70;
t90 = g(1) * t40 - g(2) * t38;
t89 = g(1) * t38 + g(2) * t40;
t88 = t69 * t101 - t71 * t99;
t87 = -g(1) * t108 + g(2) * t109;
t10 = qJD(5) * t95 - t96;
t36 = t41 * qJD(5);
t86 = t39 * t10 - t36 * t33;
t11 = qJD(3) * t36 + t88;
t35 = t39 * qJD(5);
t85 = -t41 * t11 + t35 * t31;
t84 = t69 * t12 + t71 * t13;
t14 = -t66 * t24 - t100;
t83 = -t14 * t66 + t15 * t67;
t18 = -t104 + (-pkin(6) * qJD(3) - t47) * t66;
t19 = pkin(6) * t107 + t28;
t3 = t71 * t18 - t69 * t19;
t4 = t69 * t18 + t71 * t19;
t81 = qJDD(4) + t58 - t97;
t80 = t72 * qJDD(3) - t73 * t70;
t30 = t81 - t105;
t79 = -t30 + t90;
t78 = -t35 * qJD(5) + t41 * qJDD(5);
t77 = t36 * qJD(5) + t39 * qJDD(5);
t25 = t41 * t70;
t75 = t83 + t89;
t20 = -t57 * qJDD(3) + t81;
t74 = t58 + t79 + t105;
t65 = qJDD(1) - g(3);
t64 = pkin(8) + qJ(5);
t60 = cos(t64);
t59 = sin(t64);
t37 = -t57 * qJD(3) + t91;
t29 = t31 ^ 2;
t8 = qJD(5) * t26 - t72 * t33;
t7 = -qJD(3) * t76 - qJD(5) * t25;
t2 = -t4 * qJD(5) + t93;
t1 = t3 * qJD(5) + t84;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t67 - t15 * t66 - g(3), 0, 0, 0, 0, 0, 0, t77, t78, -t85 + t86, -t1 * t41 + t2 * t39 + t3 * t36 + t4 * t35 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t87, 0, 0, 0, 0, 0, 0, t80, -qJDD(3) * t70 - t112, 0, (t70 ^ 2 + t72 ^ 2) * qJDD(2) + t87, 0, 0, 0, 0, 0, 0, t80 * t67, -t80 * t66, t110 * t112 + t70 * t92, -t120 * qJD(3) - t30 * t72 + t83 * t70 + t87, 0, 0, 0, 0, 0, 0, t8 * qJD(5) - t25 * qJDD(5) + t31 * t106 - t72 * t11, -t7 * qJD(5) + t26 * qJDD(5) + t72 * t10 + t33 * t106, -t25 * t10 + t26 * t11 - t7 * t31 - t8 * t33, -t1 * t26 + t37 * t106 - t2 * t25 - t20 * t72 + t3 * t8 + t4 * t7 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t90 + t97, -t89 - t98, 0, 0, t62 * qJDD(3), 0.2e1 * t66 * t99, 0, t63 * qJDD(3), 0, 0, t74 * t67, -t74 * t66, t91 * qJD(3) * t110 + qJ(4) * t92 + t75, t79 * pkin(3) + t75 * qJ(4) + t120 * qJD(2) - t82 * qJD(4), -t10 * t41 - t33 * t35, t85 + t86, t78, t11 * t39 + t31 * t36, -t77, 0, t117 * qJD(5) + t16 * qJDD(5) - t31 * t103 - t57 * t11 + t20 * t39 + t37 * t36 + t90 * t60, -t118 * qJD(5) - t17 * qJDD(5) + t57 * t10 - t33 * t103 + t20 * t41 - t37 * t35 - t90 * t59, -t1 * t39 + t16 * t10 - t17 * t11 - t117 * t33 - t118 * t31 - t2 * t41 + t3 * t35 - t4 * t36 + t89, t1 * t17 + t2 * t16 - t20 * t57 - t37 * t103 - g(1) * (-t111 * t38 - t40 * t57) - g(2) * (-t111 * t40 + t38 * t57) + t118 * t4 + t117 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t101, -t110 * t73, t82 * qJD(3) - t79, 0, 0, 0, 0, 0, 0, 0.2e1 * t33 * qJD(5) + t88, (-t31 - t95) * qJD(5) + t96, -t29 - t119, t3 * t33 + t4 * t31 + t20 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, -t29 + t119, (t31 - t95) * qJD(5) + t96, -t116, -t88, qJDD(5), g(3) * t60 - t37 * t33 - t89 * t59 + t93, -g(3) * t59 + t37 * t31 - t89 * t60 - t84, 0, 0;];
tau_reg = t5;
