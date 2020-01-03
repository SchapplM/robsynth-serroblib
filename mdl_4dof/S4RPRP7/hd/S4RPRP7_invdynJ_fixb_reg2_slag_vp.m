% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP7
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
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:20
% EndTime: 2019-12-31 16:47:20
% DurationCPUTime: 0.42s
% Computational Cost: add. (369->118), mult. (677->133), div. (0->0), fcn. (309->4), ass. (0->83)
t95 = qJDD(1) * qJ(2);
t44 = -pkin(1) - pkin(5);
t22 = qJDD(1) * t44 + qJDD(2);
t40 = sin(qJ(3));
t38 = t40 ^ 2;
t42 = cos(qJ(3));
t39 = t42 ^ 2;
t83 = t38 + t39;
t66 = t83 * t22;
t43 = cos(qJ(1));
t35 = g(2) * t43;
t94 = g(3) * t40 + t42 * t35;
t41 = sin(qJ(1));
t36 = g(1) * t41;
t85 = t36 - t35;
t55 = t40 * pkin(3) - t42 * qJ(4);
t17 = qJ(2) + t55;
t56 = pkin(3) * t42 + qJ(4) * t40;
t5 = t56 * qJD(3) - t42 * qJD(4) + qJD(2);
t1 = qJD(1) * t5 + qJDD(1) * t17;
t15 = t40 * t22;
t69 = qJDD(3) * qJ(4);
t23 = t44 * qJD(1) + qJD(2);
t89 = t42 * t23;
t2 = t69 + t15 + (qJD(4) + t89) * qJD(3);
t16 = t42 * t22;
t76 = qJDD(3) * pkin(3);
t91 = t40 * t23;
t3 = qJD(3) * t91 + qJDD(4) - t16 - t76;
t64 = qJD(3) * pkin(3) - qJD(4);
t6 = -t64 - t89;
t75 = qJD(3) * qJ(4);
t8 = t75 + t91;
t57 = t40 * t6 + t42 * t8;
t49 = t57 * qJD(3) + t2 * t40 - t3 * t42;
t93 = 0.2e1 * qJD(3);
t90 = t40 * t42;
t45 = qJD(3) ^ 2;
t88 = t44 * t45;
t67 = 0.2e1 * qJD(1) * qJD(2);
t87 = (t67 + t95) * qJ(2);
t86 = t43 * pkin(1) + t41 * qJ(2);
t84 = t38 - t39;
t46 = qJD(1) ^ 2;
t82 = t45 + t46;
t79 = t46 * qJ(2);
t7 = qJD(1) * t17;
t78 = t7 * qJD(1);
t77 = pkin(1) * qJDD(1);
t73 = qJDD(3) * t40;
t72 = qJDD(3) * t44;
t71 = t40 * qJDD(1);
t28 = t42 * qJDD(1);
t70 = qJD(1) * qJD(3);
t68 = t16 + t94;
t65 = qJ(2) * t70;
t18 = t83 * qJDD(1);
t63 = g(2) * (t43 * pkin(5) + t86);
t62 = qJDD(2) - t77;
t61 = t70 * t90;
t60 = -t79 - t85;
t59 = g(1) * t43 + g(2) * t41;
t54 = -t79 - t36;
t53 = t78 + t36;
t52 = t7 * t93;
t51 = -t44 * t18 + t85;
t50 = -t59 + t67 + 0.2e1 * t95;
t48 = t50 - t88;
t47 = -0.2e1 * t1 + t59 + t88;
t31 = t43 * qJ(2);
t29 = qJDD(3) * t42;
t25 = t46 * t90;
t24 = t42 * t72;
t21 = t84 * t46;
t20 = -t45 * t40 + t29;
t19 = t45 * t42 + t73;
t14 = t56 * qJD(1);
t12 = t39 * qJDD(1) - 0.2e1 * t61;
t11 = t38 * qJDD(1) + 0.2e1 * t61;
t10 = -t82 * t40 + t29;
t9 = t82 * t42 + t73;
t4 = -t40 * t28 + t84 * t70;
t13 = [0, 0, 0, 0, 0, qJDD(1), t85, t59, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t77 - t85, t50, -t62 * pkin(1) - g(1) * (-t41 * pkin(1) + t31) - g(2) * t86 + t87, t12, 0.2e1 * t4, t20, t11, -t19, 0, t48 * t40 + 0.2e1 * t42 * t65 + t24, (-0.2e1 * t65 - t72) * t40 + t48 * t42, t51 - t66, -g(1) * (t44 * t41 + t31) - t63 + t44 * t66 + t87, t12, t20, -0.2e1 * t4, 0, t19, t11, -t47 * t40 + t42 * t52 + t24, -t49 + t51, (t52 + t72) * t40 + t47 * t42, t1 * t17 + t7 * t5 - g(1) * (t55 * t43 + t31) - t63 + (-g(1) * t44 - g(2) * t55) * t41 + t49 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t46, t60 + t62, 0, 0, 0, 0, 0, 0, t10, -t9, -t18, t66 + t60, 0, 0, 0, 0, 0, 0, t10, -t18, t9, t49 - t78 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t21, t28, -t25, -t71, qJDD(3), t54 * t42 + t68, g(3) * t42 - t15 + (-t54 - t35) * t40, 0, 0, t25, t28, t21, qJDD(3), t71, -t25, -t42 * t36 + 0.2e1 * t76 - qJDD(4) + (-t14 * t40 - t42 * t7) * qJD(1) + t68, -t56 * qJDD(1) + ((t8 - t75) * t42 + (t6 + t64) * t40) * qJD(1), 0.2e1 * t69 + qJD(4) * t93 + t15 + (qJD(1) * t14 - g(3)) * t42 + (-t53 + t35) * t40, -t3 * pkin(3) + g(3) * t55 + t2 * qJ(4) + t8 * qJD(4) - t7 * t14 - t57 * t23 - t85 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t25, t28, -t39 * t46 - t45, -t8 * qJD(3) + t53 * t42 + t3 - t94;];
tau_reg = t13;
