% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:03
% EndTime: 2019-12-31 16:22:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (251->68), mult. (375->86), div. (0->0), fcn. (180->4), ass. (0->49)
t27 = qJD(2) ^ 2;
t18 = pkin(6) + qJ(2);
t16 = sin(t18);
t17 = cos(t18);
t54 = g(1) * t16 - g(2) * t17;
t37 = -t27 * qJ(3) - t54;
t23 = sin(qJ(4));
t24 = cos(qJ(4));
t25 = -pkin(2) - pkin(5);
t9 = t25 * qJD(2) + qJD(3);
t3 = -t23 * qJD(1) + t24 * t9;
t49 = t3 * qJD(4);
t8 = t25 * qJDD(2) + qJDD(3);
t1 = t24 * qJDD(1) + t23 * t8 + t49;
t4 = t24 * qJD(1) + t23 * t9;
t48 = t4 * qJD(4);
t5 = t24 * t8;
t2 = -t23 * qJDD(1) - t48 + t5;
t29 = -(t23 * t3 - t24 * t4) * qJD(4) + t1 * t23 + t2 * t24;
t57 = 0.2e1 * qJ(3);
t56 = t23 * t24;
t55 = t17 * pkin(2) + t16 * qJ(3);
t41 = 2 * qJD(2) * qJD(3);
t53 = (qJ(3) * qJDD(2) + t41) * qJ(3);
t20 = t23 ^ 2;
t21 = t24 ^ 2;
t52 = t20 - t21;
t26 = qJD(4) ^ 2;
t51 = -t26 - t27;
t47 = pkin(2) * qJDD(2);
t22 = qJDD(1) - g(3);
t46 = qJDD(4) * t23;
t45 = qJDD(4) * t24;
t44 = t24 * qJDD(2);
t43 = qJD(2) * qJD(4);
t42 = t27 * t56;
t40 = (-t20 - t21) * qJDD(2);
t39 = qJDD(3) - t47;
t38 = t43 * t56;
t36 = g(1) * t17 + g(2) * t16;
t33 = -qJD(4) * t9 - t22;
t32 = qJDD(4) * t25 + t43 * t57;
t31 = -qJD(1) * qJD(4) + t37;
t30 = qJDD(2) * t57 - t36 + t41;
t28 = -t25 * t26 + t30;
t11 = t17 * qJ(3);
t7 = -t26 * t23 + t45;
t6 = -t26 * t24 - t46;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, t6, -t7, 0, t1 * t24 - t2 * t23 - g(3) + (-t23 * t4 - t24 * t3) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t54, t36, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t47 - t54, t30, -t39 * pkin(2) - g(1) * (-t16 * pkin(2) + t11) - g(2) * t55 + t53, t21 * qJDD(2) - 0.2e1 * t38, -0.2e1 * t23 * t44 + 0.2e1 * t52 * t43, t7, t20 * qJDD(2) + 0.2e1 * t38, t6, 0, t28 * t23 + t32 * t24, -t32 * t23 + t28 * t24, t25 * t40 - t29 + t54, -g(1) * (t25 * t16 + t11) - g(2) * (t17 * pkin(5) + t55) + t29 * t25 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t27, t37 + t39, 0, 0, 0, 0, 0, 0, t51 * t23 + t45, t51 * t24 - t46, t40, t29 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t52 * t27, t44, -t42, -t23 * qJDD(2), qJDD(4), t33 * t23 + t31 * t24 + t48 + t5, t49 + t33 * t24 + (-t31 - t8) * t23, 0, 0;];
tau_reg = t10;
