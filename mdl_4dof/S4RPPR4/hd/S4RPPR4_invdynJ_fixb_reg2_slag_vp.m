% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:56
% EndTime: 2019-12-31 16:38:57
% DurationCPUTime: 0.32s
% Computational Cost: add. (342->81), mult. (530->107), div. (0->0), fcn. (269->8), ass. (0->58)
t27 = sin(pkin(6));
t15 = pkin(1) * t27 + qJ(3);
t11 = qJD(1) * t15;
t23 = qJ(1) + pkin(6);
t20 = sin(t23);
t21 = cos(t23);
t66 = g(1) * t20 - g(2) * t21;
t42 = -t11 * qJD(1) - t66;
t53 = qJD(3) * qJD(1);
t56 = t15 * qJDD(1);
t8 = t53 + t56;
t61 = pkin(1) * qJDD(1);
t28 = cos(pkin(6));
t19 = -pkin(1) * t28 - pkin(2);
t69 = t19 * qJDD(1);
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t12 = -pkin(5) + t19;
t6 = qJDD(1) * t12 + qJDD(3);
t7 = qJD(1) * t12 + qJD(3);
t3 = -qJD(2) * t29 + t31 * t7;
t63 = t3 * qJD(4);
t1 = t31 * qJDD(2) + t29 * t6 + t63;
t5 = t31 * t6;
t4 = qJD(2) * t31 + t29 * t7;
t62 = t4 * qJD(4);
t2 = -t29 * qJDD(2) + t5 - t62;
t37 = -(t29 * t3 - t31 * t4) * qJD(4) + t1 * t29 + t2 * t31;
t68 = 0.2e1 * qJD(4) * t11 + qJDD(4) * t12;
t67 = t29 * t31;
t24 = t29 ^ 2;
t25 = t31 ^ 2;
t65 = t24 - t25;
t33 = qJD(4) ^ 2;
t34 = qJD(1) ^ 2;
t64 = -t33 - t34;
t26 = qJDD(2) - g(3);
t58 = qJDD(4) * t29;
t57 = qJDD(4) * t31;
t55 = t31 * qJDD(1);
t54 = qJD(1) * qJD(4);
t52 = t34 * t67;
t32 = cos(qJ(1));
t51 = t32 * pkin(1) + t21 * pkin(2) + t20 * qJ(3);
t30 = sin(qJ(1));
t50 = -t30 * pkin(1) + t21 * qJ(3);
t48 = (-t24 - t25) * qJDD(1);
t47 = t54 * t67;
t46 = g(1) * t21 + g(2) * t20;
t45 = g(1) * t30 - g(2) * t32;
t41 = -qJD(4) * t7 - t26;
t40 = t11 * qJD(3) + t8 * t15;
t39 = qJDD(3) + t69;
t38 = -qJD(2) * qJD(4) + t42;
t36 = -t12 * t33 - t46 + 0.2e1 * t8;
t10 = -t29 * t33 + t57;
t9 = -t31 * t33 - t58;
t13 = [0, 0, 0, 0, 0, qJDD(1), t45, g(1) * t32 + g(2) * t30, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t28 * t61 + t66, -0.2e1 * t27 * t61 + t46, 0, (t45 + (t27 ^ 2 + t28 ^ 2) * t61) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t66 + 0.2e1 * t69, -t46 + 0.2e1 * t53 + 0.2e1 * t56, t39 * t19 - g(1) * (-pkin(2) * t20 + t50) - g(2) * t51 + t40, qJDD(1) * t25 - 0.2e1 * t47, -0.2e1 * t29 * t55 + 0.2e1 * t54 * t65, t10, qJDD(1) * t24 + 0.2e1 * t47, t9, 0, t36 * t29 + t31 * t68, -t29 * t68 + t36 * t31, t12 * t48 - t37 + t66, -g(1) * ((-pkin(2) - pkin(5)) * t20 + t50) - g(2) * (pkin(5) * t21 + t51) + t37 * t12 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, t9, -t10, 0, t1 * t31 - t2 * t29 - g(3) + (-t29 * t4 - t3 * t31) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t34, t39 + t42, 0, 0, 0, 0, 0, 0, t29 * t64 + t57, t31 * t64 - t58, t48, t37 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t65 * t34, t55, -t52, -t29 * qJDD(1), qJDD(4), t29 * t41 + t31 * t38 + t5 + t62, t63 + t41 * t31 + (-t38 - t6) * t29, 0, 0;];
tau_reg = t13;
