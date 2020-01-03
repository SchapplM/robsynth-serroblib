% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP6
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:45
% EndTime: 2019-12-31 16:30:46
% DurationCPUTime: 0.57s
% Computational Cost: add. (422->147), mult. (990->197), div. (0->0), fcn. (601->6), ass. (0->99)
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t67 = t52 * pkin(3) + t50 * qJ(4);
t63 = pkin(2) + t67;
t117 = t63 * qJD(2);
t51 = sin(qJ(2));
t32 = qJD(2) * pkin(5) + t51 * qJD(1);
t53 = cos(qJ(2));
t91 = t53 * qJD(1);
t97 = qJD(2) * pkin(2);
t33 = -t91 - t97;
t46 = t50 ^ 2;
t47 = t52 ^ 2;
t99 = t46 + t47;
t77 = t99 * t53;
t116 = t32 * t77 + t33 * t51;
t48 = sin(pkin(6));
t49 = cos(pkin(6));
t68 = g(1) * t49 + g(2) * t48;
t66 = pkin(3) * t50 - qJ(4) * t52;
t11 = t66 * qJD(3) - t50 * qJD(4);
t115 = qJD(2) * t11 - t63 * qJDD(2);
t83 = qJD(1) * qJD(2);
t113 = (t68 + t83) * t51;
t102 = t52 * t32;
t92 = qJDD(3) * pkin(3);
t112 = qJD(3) * t102 - t92;
t111 = pkin(2) * t51;
t110 = pkin(5) * t53;
t107 = g(3) * t51;
t105 = t50 * t32;
t104 = t50 * t53;
t103 = t51 * t52;
t101 = t52 * t53;
t100 = -t46 + t47;
t54 = qJD(3) ^ 2;
t55 = qJD(2) ^ 2;
t98 = t54 + t55;
t96 = pkin(5) * qJDD(3);
t94 = qJD(2) * t50;
t93 = qJDD(2) * pkin(2);
t90 = qJDD(1) - g(3);
t89 = qJD(3) * qJ(4);
t88 = qJDD(3) * t50;
t87 = t46 * qJDD(2);
t86 = t47 * qJDD(2);
t85 = t52 * qJDD(2);
t84 = t53 * qJDD(2);
t82 = qJD(2) * qJD(3);
t81 = qJDD(3) * qJ(4);
t80 = g(3) * (t53 * pkin(2) + t51 * pkin(5));
t19 = qJDD(2) * pkin(5) + t51 * qJDD(1) + t53 * t83;
t79 = t99 * t19;
t78 = t99 * t51;
t76 = t50 * t82;
t40 = t51 * t83;
t75 = t33 - t97;
t8 = -t91 - t117;
t74 = t8 - t117;
t73 = t50 * qJD(3) * t91 + t68 * t103 + t52 * t40;
t72 = -qJD(3) * pkin(3) + qJD(4);
t71 = -t53 * qJDD(1) + t40;
t70 = t52 * t76;
t69 = pkin(5) * t54 + g(3) * t53;
t12 = t72 + t105;
t17 = t89 + t102;
t65 = t12 * t50 + t17 * t52;
t10 = t52 * t19;
t14 = t48 * t101 - t49 * t50;
t16 = t49 * t101 + t48 * t50;
t64 = g(1) * t16 + g(2) * t14 - t10;
t13 = t48 * t104 + t49 * t52;
t15 = t49 * t104 - t48 * t52;
t9 = t50 * t19;
t62 = g(1) * t15 + g(2) * t13 + t50 * t107 - t9;
t60 = -qJDD(4) + t62;
t18 = t71 - t93;
t59 = -t18 - t69 + t93;
t1 = t115 + t71;
t58 = -t1 - t115 - t69;
t4 = t81 + t10 + (qJD(4) - t105) * qJD(3);
t5 = qJDD(4) + t9 + t112;
t57 = t4 * t52 + t5 * t50 + (t12 * t52 - t17 * t50) * qJD(3);
t56 = (-t99 * t83 - t68) * t53 - t107 + (t87 + t86) * pkin(5);
t43 = t50 * qJDD(2);
t36 = t49 * t110;
t35 = t48 * t110;
t34 = t50 * t55 * t52;
t27 = t100 * t55;
t26 = qJDD(3) * t52 - t54 * t50;
t25 = t54 * t52 + t88;
t23 = t66 * qJD(2);
t21 = -0.2e1 * t70 + t86;
t20 = 0.2e1 * t70 + t87;
t7 = t100 * t82 + t50 * t85;
t6 = qJDD(2) * t78 + t55 * t77;
t3 = (-0.2e1 * t76 + t85) * t53 + (-t98 * t52 - t88) * t51;
t2 = (-qJDD(3) * t51 - 0.2e1 * t53 * t82) * t52 + (t98 * t51 - t84) * t50;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, 0, 0, 0, 0, 0, -t55 * t51 + t84, -qJDD(2) * t51 - t55 * t53, 0, -g(3) + (t51 ^ 2 + t53 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t3, t2, t6, t116 * qJD(2) - t18 * t53 + t19 * t78 - g(3), 0, 0, 0, 0, 0, 0, t3, t6, -t2, -g(3) + (t65 * qJD(2) - t1) * t53 + (qJD(2) * t8 + t57) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t68 * t51 + t90 * t53, -t90 * t51 + t68 * t53, 0, 0, t20, 0.2e1 * t7, t25, t21, t26, 0, (t75 * qJD(3) - t96) * t50 + t59 * t52 + t73, (-t96 + (t75 + t91) * qJD(3)) * t52 + (-t59 - t113) * t50, t79 + t56, -t18 * pkin(2) - g(1) * (-t49 * t111 + t36) - g(2) * (-t48 * t111 + t35) - t80 + pkin(5) * t79 - t116 * qJD(1), t20, t25, -0.2e1 * t7, 0, -t26, t21, (t74 * qJD(3) - t96) * t50 + t58 * t52 + t73, t56 + t57, (t96 + (-t74 - t91) * qJD(3)) * t52 + (t58 + t113) * t50, -t1 * t63 + t8 * t11 - g(1) * t36 - g(2) * t35 - t80 + (-g(3) * t67 - t65 * qJD(1)) * t53 + t57 * pkin(5) + (-t8 * qJD(1) + t68 * t63) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t27, t43, t34, t85, qJDD(3), -t33 * t94 + t62, (-qJD(2) * t33 + t107) * t52 + t64, 0, 0, -t34, t43, t27, qJDD(3), -t85, t34, 0.2e1 * t92 + (t23 * t52 - t50 * t8) * qJD(2) + t60, -t66 * qJDD(2) + ((t17 - t89) * t50 + (-t12 + t72) * t52) * qJD(2), -g(3) * t103 + 0.2e1 * t81 + 0.2e1 * qJD(3) * qJD(4) + (t23 * t50 + t52 * t8) * qJD(2) - t64, t4 * qJ(4) - t5 * pkin(3) - t8 * t23 - t12 * t102 - g(1) * (-t15 * pkin(3) + t16 * qJ(4)) - g(2) * (-t13 * pkin(3) + t14 * qJ(4)) + t66 * t107 + (qJD(4) + t105) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t34, t43, -t46 * t55 - t54, -t17 * qJD(3) + t8 * t94 + t112 - t60;];
tau_reg = t22;
