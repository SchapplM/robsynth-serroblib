% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:30
% EndTime: 2019-12-31 16:56:32
% DurationCPUTime: 0.70s
% Computational Cost: add. (516->149), mult. (1050->222), div. (0->0), fcn. (653->6), ass. (0->93)
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t59 = g(1) * t38 - g(2) * t41;
t111 = qJDD(2) - t59;
t84 = pkin(1) * qJDD(1);
t110 = t84 - t111;
t42 = -pkin(1) - pkin(5);
t25 = t42 * qJD(1) + qJD(2);
t40 = cos(qJ(3));
t91 = t40 * t25;
t15 = -qJD(3) * pkin(3) - t91;
t70 = qJD(1) * qJD(3);
t67 = t40 * t70;
t37 = sin(qJ(3));
t72 = t37 * qJDD(1);
t16 = qJDD(4) + t67 + t72;
t22 = t37 * pkin(3) - t40 * pkin(6) + qJ(2);
t27 = t37 * qJD(1) + qJD(4);
t24 = t42 * qJDD(1) + qJDD(2);
t79 = qJD(3) * t40;
t9 = t22 * qJD(1);
t66 = -qJDD(3) * pkin(6) - qJD(4) * t9 - t37 * t24 - t25 * t79;
t80 = qJD(3) * t37;
t7 = -qJDD(3) * pkin(3) - t40 * t24 + t25 * t80;
t109 = -(qJD(4) * t22 + t42 * t79) * t27 + t7 * t40 + (-t15 * qJD(3) - t42 * t16 + t66) * t37;
t39 = cos(qJ(4));
t36 = sin(qJ(4));
t75 = t36 * qJD(3);
t83 = qJD(1) * t40;
t20 = t39 * t83 + t75;
t71 = t40 * qJDD(1);
t97 = t36 * t37;
t4 = t20 * qJD(4) - t39 * qJDD(3) + t36 * t71 - t70 * t97;
t106 = g(3) * t37;
t61 = pkin(3) * t40 + pkin(6) * t37;
t108 = (pkin(6) * qJD(4) + t61 * qJD(1)) * t27 + t59 * t40 + t7 - t106;
t105 = g(3) * t40;
t74 = t39 * qJD(3);
t76 = qJD(4) * t40;
t51 = -t36 * t76 - t37 * t74;
t3 = t51 * qJD(1) + qJD(4) * t74 + t36 * qJDD(3) + t39 * t71;
t104 = t3 * t36;
t103 = t40 * t3;
t18 = t36 * t83 - t74;
t102 = t18 * t27;
t101 = t20 * t27;
t100 = t22 * t16;
t99 = t27 * t39;
t98 = t36 * t16;
t96 = t37 * t25;
t95 = t37 * t42;
t94 = t38 * t36;
t93 = t38 * t39;
t92 = t39 * t16;
t90 = t40 * t42;
t89 = t41 * t36;
t88 = t41 * t39;
t35 = t40 ^ 2;
t87 = t37 ^ 2 - t35;
t43 = qJD(3) ^ 2;
t44 = qJD(1) ^ 2;
t86 = -t43 - t44;
t85 = t44 * qJ(2);
t82 = qJD(3) * t18;
t81 = qJD(3) * t20;
t78 = qJD(4) * t37;
t77 = qJD(4) * t39;
t73 = qJDD(3) * t37;
t69 = qJDD(1) * qJ(2);
t14 = qJD(3) * pkin(6) + t96;
t17 = t61 * qJD(3) + qJD(2);
t6 = t17 * qJD(1) + t22 * qJDD(1);
t64 = qJD(4) * t14 - t6;
t62 = qJD(1) + t78;
t60 = g(1) * t41 + g(2) * t38;
t58 = t66 + t105;
t56 = t27 * t77 + t98;
t55 = -qJD(4) * t36 * t27 + t92;
t52 = 0.2e1 * qJ(2) * t70 + qJDD(3) * t42;
t50 = 0.2e1 * qJD(1) * qJD(2) - t60;
t49 = -t24 + t59 + t85;
t48 = t50 + 0.2e1 * t69;
t47 = -pkin(6) * t16 + (t15 + t91) * t27;
t45 = -t42 * t43 + t48;
t31 = qJDD(3) * t40;
t13 = t37 * t88 - t94;
t12 = t37 * t89 + t93;
t11 = t37 * t93 + t89;
t10 = -t37 * t94 + t88;
t5 = t39 * t6;
t2 = t39 * t14 + t36 * t9;
t1 = -t36 * t14 + t39 * t9;
t8 = [qJDD(1), t59, t60, -0.2e1 * t84 + t111, t48, t110 * pkin(1) + (t50 + t69) * qJ(2), t35 * qJDD(1) - 0.2e1 * t37 * t67, -0.2e1 * t37 * t71 + 0.2e1 * t87 * t70, -t43 * t37 + t31, -t43 * t40 - t73, 0, t45 * t37 + t52 * t40, -t52 * t37 + t45 * t40, t39 * t103 + t51 * t20, (t18 * t39 + t20 * t36) * t80 + (-t104 - t39 * t4 + (t18 * t36 - t20 * t39) * qJD(4)) * t40, (-t27 * t74 + t3) * t37 + (t55 + t81) * t40, (t27 * t75 - t4) * t37 + (-t56 - t82) * t40, t16 * t37 + t27 * t79, -t4 * t90 - g(1) * t13 - g(2) * t11 + t5 * t37 + (t1 * t40 + t18 * t95) * qJD(3) + (t100 + t17 * t27 + (t15 * t40 + (-t27 * t42 - t14) * t37) * qJD(4)) * t39 + t109 * t36, -t3 * t90 + g(1) * t12 - g(2) * t10 + (-t2 * t40 + t20 * t95) * qJD(3) + (-(-t42 * t78 + t17) * t27 - t100 + t64 * t37 - t15 * t76) * t36 + t109 * t39; 0, 0, 0, qJDD(1), -t44, -t85 - t110, 0, 0, 0, 0, 0, t86 * t37 + t31, t86 * t40 - t73, 0, 0, 0, 0, 0, -t40 * t4 + (t82 - t98) * t37 + (-t62 * t39 - t40 * t75) * t27, -t103 + (t81 - t92) * t37 + (t62 * t36 - t40 * t74) * t27; 0, 0, 0, 0, 0, 0, t40 * t44 * t37, -t87 * t44, t71, -t72, qJDD(3), -t49 * t40 + t106, t49 * t37 + t105, t20 * t99 + t104, (t3 - t102) * t39 + (-t4 - t101) * t36, (-t20 * t40 + t37 * t99) * qJD(1) + t56, (t18 * t40 - t27 * t97) * qJD(1) + t55, -t27 * t83, -pkin(3) * t4 - t1 * t83 - t108 * t39 - t18 * t96 + t47 * t36, -pkin(3) * t3 + t108 * t36 + t2 * t83 - t20 * t96 + t47 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t18, -t18 ^ 2 + t20 ^ 2, t3 + t102, t101 - t4, t16, -g(1) * t10 - g(2) * t12 - t14 * t77 - t15 * t20 + t2 * t27 + t58 * t36 + t5, g(1) * t11 - g(2) * t13 + t1 * t27 + t15 * t18 + t64 * t36 + t58 * t39;];
tau_reg = t8;
