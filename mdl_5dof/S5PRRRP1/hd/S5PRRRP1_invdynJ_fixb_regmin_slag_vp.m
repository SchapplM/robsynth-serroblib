% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:16
% EndTime: 2019-12-05 16:40:20
% DurationCPUTime: 0.59s
% Computational Cost: add. (719->145), mult. (1016->184), div. (0->0), fcn. (527->8), ass. (0->94)
t112 = qJ(5) + pkin(7);
t50 = qJDD(2) + qJDD(3);
t103 = t50 * pkin(3);
t51 = pkin(8) + qJ(2);
t45 = qJ(3) + t51;
t37 = cos(t45);
t105 = g(2) * t37;
t59 = cos(qJ(3));
t107 = pkin(2) * t59;
t57 = sin(qJ(3));
t91 = qJD(3) * t57;
t83 = pkin(2) * t91;
t98 = -qJD(2) * t83 + qJDD(2) * t107;
t111 = -t98 - t103 + t105;
t36 = sin(t45);
t110 = g(1) * t37 + g(2) * t36;
t33 = g(1) * t36;
t109 = t105 - t33;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t52 = qJD(2) + qJD(3);
t95 = pkin(2) * qJD(2);
t84 = t57 * t95;
t75 = t112 * t52 + t84;
t9 = t58 * qJD(1) - t75 * t56;
t94 = qJD(4) * pkin(4);
t5 = t9 + t94;
t108 = t5 - t9;
t106 = pkin(3) * t52;
t104 = g(3) * t58;
t102 = t56 * t50;
t101 = t58 * t50;
t92 = qJD(2) * t59;
t80 = pkin(2) * t92;
t24 = -t80 - t106;
t89 = qJD(4) * t56;
t100 = t24 * t89 + t58 * t33;
t53 = t56 ^ 2;
t54 = t58 ^ 2;
t97 = -t53 - t54;
t96 = t53 - t54;
t38 = pkin(2) * t57 + pkin(7);
t93 = -qJ(5) - t38;
t90 = qJD(3) * t59;
t88 = qJDD(1) - g(3);
t87 = qJDD(2) * t57;
t14 = t50 * pkin(7) + (qJD(2) * t90 + t87) * pkin(2);
t63 = qJ(5) * t50 + qJD(1) * qJD(4) + qJD(5) * t52 + t14;
t67 = t75 * qJD(4);
t3 = (qJDD(1) - t67) * t56 + t63 * t58;
t86 = t3 * t58 - t110;
t85 = t24 * qJD(4) * t58 + t111 * t56;
t82 = pkin(2) * t90;
t81 = pkin(4) * t89;
t79 = t52 * t91;
t78 = t52 * t89;
t39 = pkin(4) * t58 + pkin(3);
t76 = t112 * t36 + t37 * t39;
t74 = qJD(4) * t112;
t73 = qJD(4) * t93;
t72 = t52 * t84;
t10 = t56 * qJD(1) + t75 * t58;
t70 = t10 * t58 - t5 * t56;
t69 = t112 * t37 - t36 * t39;
t68 = -t98 + t109;
t66 = -t24 * t52 + t110 - t14;
t60 = qJD(4) ^ 2;
t65 = pkin(7) * t60 - t103 - t72;
t4 = pkin(4) * t78 - t39 * t50 + qJDD(5) - t98;
t40 = -pkin(3) - t107;
t64 = pkin(2) * t79 + t38 * t60 + t40 * t50;
t62 = -pkin(7) * qJDD(4) + (t80 - t106) * qJD(4);
t61 = -qJDD(4) * t38 + (t40 * t52 - t82) * qJD(4);
t49 = t52 ^ 2;
t48 = t58 * qJ(5);
t46 = t58 * qJD(5);
t44 = t58 * qJDD(1);
t43 = cos(t51);
t42 = sin(t51);
t30 = pkin(7) * t58 + t48;
t29 = t112 * t56;
t26 = qJDD(4) * t58 - t56 * t60;
t25 = qJDD(4) * t56 + t58 * t60;
t21 = t38 * t58 + t48;
t20 = t93 * t56;
t17 = -t56 * qJD(5) - t58 * t74;
t16 = -t56 * t74 + t46;
t15 = t50 * t53 + 0.2e1 * t58 * t78;
t12 = -t39 * t52 + qJD(5) - t80;
t8 = -0.2e1 * t96 * t52 * qJD(4) + 0.2e1 * t56 * t101;
t7 = (-qJD(5) - t82) * t56 + t58 * t73;
t6 = t56 * t73 + t58 * t82 + t46;
t2 = qJDD(4) * pkin(4) - t63 * t56 - t58 * t67 + t44;
t1 = [t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t70 * qJD(4) + t2 * t58 + t3 * t56 - g(3); 0, qJDD(2), g(1) * t42 - g(2) * t43, g(1) * t43 + g(2) * t42, t50, (t50 * t59 - t79) * pkin(2) - t68, ((-qJDD(2) - t50) * t57 + (-qJD(2) - t52) * t90) * pkin(2) + t110, t15, t8, t25, t26, 0, t61 * t56 + (-t64 - t111) * t58 + t100, t61 * t58 + (t64 - t33) * t56 + t85, (t21 * t50 + t52 * t6 + (-t20 * t52 - t5) * qJD(4)) * t58 + (-t20 * t50 - t52 * t7 - t2 + (-t21 * t52 - t10) * qJD(4)) * t56 + t86, t3 * t21 + t10 * t6 + t2 * t20 + t5 * t7 + t4 * (-t39 - t107) + t12 * (t81 + t83) - g(1) * (-pkin(2) * t42 + t69) - g(2) * (pkin(2) * t43 + t76); 0, 0, 0, 0, t50, -t68 + t72, (-t87 + (-qJD(3) + t52) * t92) * pkin(2) + t110, t15, t8, t25, t26, 0, t62 * t56 + (-t65 - t111) * t58 + t100, t62 * t58 + (t65 - t33) * t56 + t85, (-qJD(4) * t5 + t30 * t50) * t58 + (-qJD(4) * t10 + t29 * t50 - t2) * t56 + (t16 * t58 - t17 * t56 + (t29 * t58 - t30 * t56) * qJD(4) + t97 * t80) * t52 + t86, t3 * t30 + t10 * t16 - t2 * t29 + t5 * t17 - t4 * t39 + t12 * t81 - g(1) * t69 - g(2) * t76 + (-t12 * t57 - t70 * t59) * t95; 0, 0, 0, 0, 0, 0, 0, -t56 * t49 * t58, t96 * t49, t102, t101, qJDD(4), t66 * t56 - t104 + t44, -t88 * t56 + t66 * t58, -pkin(4) * t102 + (-t94 + t108) * t58 * t52, t108 * t10 + (-t104 + t2 + (-t12 * t52 + t110) * t56) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t49, -t70 * t52 + t109 + t4;];
tau_reg = t1;
