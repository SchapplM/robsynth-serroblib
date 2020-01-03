% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:18
% EndTime: 2019-12-31 16:47:19
% DurationCPUTime: 0.36s
% Computational Cost: add. (324->105), mult. (570->127), div. (0->0), fcn. (262->4), ass. (0->70)
t34 = -pkin(1) - pkin(5);
t15 = t34 * qJDD(1) + qJDD(2);
t30 = sin(qJ(3));
t12 = t30 * t15;
t55 = qJDD(3) * qJ(4);
t16 = t34 * qJD(1) + qJD(2);
t32 = cos(qJ(3));
t74 = t16 * t32;
t2 = t55 + t12 + (qJD(4) + t74) * qJD(3);
t13 = t32 * t15;
t62 = qJDD(3) * pkin(3);
t75 = t16 * t30;
t3 = qJD(3) * t75 + qJDD(4) - t13 - t62;
t50 = qJD(3) * pkin(3) - qJD(4);
t5 = -t50 - t74;
t61 = qJD(3) * qJ(4);
t7 = t61 + t75;
t46 = t30 * t5 + t32 * t7;
t39 = t46 * qJD(3) + t2 * t30 - t3 * t32;
t33 = cos(qJ(1));
t26 = g(2) * t33;
t31 = sin(qJ(1));
t27 = g(1) * t31;
t70 = t27 - t26;
t79 = t39 - t70;
t78 = g(3) * t30 + t32 * t26;
t44 = pkin(3) * t30 - qJ(4) * t32;
t14 = qJ(2) + t44;
t45 = pkin(3) * t32 + qJ(4) * t30;
t4 = t45 * qJD(3) - t32 * qJD(4) + qJD(2);
t1 = qJD(1) * t4 + qJDD(1) * t14;
t77 = 0.2e1 * qJD(3);
t73 = t30 * t32;
t35 = qJD(3) ^ 2;
t72 = t34 * t35;
t71 = t33 * pkin(1) + t31 * qJ(2);
t28 = t30 ^ 2;
t29 = t32 ^ 2;
t69 = t28 - t29;
t36 = qJD(1) ^ 2;
t68 = t35 + t36;
t67 = qJ(2) * t36;
t6 = qJD(1) * t14;
t64 = qJD(1) * t6;
t63 = pkin(1) * qJDD(1);
t59 = qJDD(3) * t30;
t58 = qJDD(3) * t34;
t20 = t32 * qJDD(1);
t57 = qJD(1) * qJD(3);
t56 = qJDD(1) * qJ(2);
t54 = t13 + t78;
t53 = 0.2e1 * qJD(1) * qJD(2);
t52 = qJ(2) * t57;
t51 = (-t28 - t29) * qJDD(1);
t49 = qJDD(2) - t63;
t48 = g(1) * t33 + g(2) * t31;
t43 = -t67 - t27;
t42 = t64 + t27;
t41 = t6 * t77;
t40 = -t48 + t53 + 0.2e1 * t56;
t38 = t40 - t72;
t37 = -0.2e1 * t1 + t48 + t72;
t23 = t33 * qJ(2);
t21 = qJDD(3) * t32;
t18 = t36 * t73;
t17 = t32 * t58;
t11 = t45 * qJD(1);
t9 = -t68 * t30 + t21;
t8 = t68 * t32 + t59;
t10 = [qJDD(1), t70, t48, qJDD(2) - 0.2e1 * t63 - t70, t40, -t49 * pkin(1) - g(1) * (-pkin(1) * t31 + t23) - g(2) * t71 + (t53 + t56) * qJ(2), qJDD(1) * t29 - 0.2e1 * t57 * t73, -0.2e1 * t30 * t20 + 0.2e1 * t69 * t57, -t30 * t35 + t21, -t32 * t35 - t59, 0, t38 * t30 + 0.2e1 * t32 * t52 + t17, (-0.2e1 * t52 - t58) * t30 + t38 * t32, -t37 * t30 + t32 * t41 + t17, t34 * t51 - t79, (t41 + t58) * t30 + t37 * t32, t1 * t14 + t6 * t4 - g(1) * (t44 * t33 + t23) - g(2) * (pkin(5) * t33 + t71) + (-g(1) * t34 - g(2) * t44) * t31 + t39 * t34; 0, 0, 0, qJDD(1), -t36, t49 - t67 - t70, 0, 0, 0, 0, 0, t9, -t8, t9, t51, t8, -t64 + t79; 0, 0, 0, 0, 0, 0, t18, -t69 * t36, t20, -t30 * qJDD(1), qJDD(3), t43 * t32 + t54, g(3) * t32 - t12 + (-t43 - t26) * t30, -t32 * t27 + 0.2e1 * t62 - qJDD(4) + (-t11 * t30 - t32 * t6) * qJD(1) + t54, -t45 * qJDD(1) + ((t7 - t61) * t32 + (t5 + t50) * t30) * qJD(1), 0.2e1 * t55 + qJD(4) * t77 + t12 + (qJD(1) * t11 - g(3)) * t32 + (-t42 + t26) * t30, -t3 * pkin(3) + g(3) * t44 + t2 * qJ(4) + t7 * qJD(4) - t6 * t11 - t46 * t16 - t70 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t18, t20, -t29 * t36 - t35, -qJD(3) * t7 + t42 * t32 + t3 - t78;];
tau_reg = t10;
