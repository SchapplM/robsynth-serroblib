% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPR1
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:04
% EndTime: 2019-03-08 18:21:04
% DurationCPUTime: 0.18s
% Computational Cost: add. (240->61), mult. (299->70), div. (0->0), fcn. (150->4), ass. (0->39)
t44 = qJD(2) - qJD(4);
t34 = -pkin(2) - pkin(3);
t18 = t34 * qJD(2) + qJD(3);
t32 = sin(qJ(4));
t33 = cos(qJ(4));
t47 = qJ(3) * qJD(2);
t5 = t18 * t33 - t32 * t47;
t58 = t44 * t5;
t57 = -qJD(4) * t47 + t34 * qJDD(2) + qJDD(3);
t45 = qJD(2) * qJD(3);
t46 = qJ(3) * qJDD(2);
t56 = qJD(4) * t18 + t45 + t46;
t2 = -t56 * t32 + t57 * t33;
t6 = t18 * t32 + t33 * t47;
t55 = -t44 * t6 + t2;
t29 = pkin(6) + qJ(2);
t26 = sin(t29);
t27 = cos(t29);
t7 = -t26 * t32 - t27 * t33;
t8 = -t26 * t33 + t27 * t32;
t54 = g(1) * t8 - g(2) * t7;
t51 = t27 * pkin(2) + t26 * qJ(3);
t50 = g(1) * t26 - g(2) * t27;
t49 = pkin(2) * qJDD(2);
t43 = 0.2e1 * t45;
t41 = t44 ^ 2;
t40 = qJDD(3) - t49;
t39 = g(1) * t27 + g(2) * t26;
t13 = qJ(3) * t33 + t32 * t34;
t12 = -qJ(3) * t32 + t33 * t34;
t1 = t57 * t32 + t56 * t33;
t36 = -g(1) * t7 - g(2) * t8 - t1;
t35 = qJD(2) ^ 2;
t31 = qJDD(1) - g(3);
t28 = qJDD(2) - qJDD(4);
t20 = t27 * qJ(3);
t4 = -qJD(3) * t32 - t13 * qJD(4);
t3 = qJD(3) * t33 + t12 * qJD(4);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t50, t39, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -qJDD(3) + 0.2e1 * t49 + t50, 0, -t39 + t43 + 0.2e1 * t46, -t40 * pkin(2) - g(1) * (-t26 * pkin(2) + t20) - g(2) * t51 + (t43 + t46) * qJ(3), 0, 0, 0, 0, 0, t28, -t12 * t28 - t4 * t44 - t2 - t54, t13 * t28 + t3 * t44 - t36, 0, t1 * t13 + t6 * t3 + t2 * t12 + t5 * t4 - g(1) * (t34 * t26 + t20) - g(2) * (pkin(3) * t27 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t35, -qJ(3) * t35 + t40 - t50, 0, 0, 0, 0, 0, 0, -t28 * t33 - t32 * t41, t28 * t32 - t33 * t41, 0, t55 * t33 + (t1 + t58) * t32 - t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t54 + t55, t36 - t58, 0, 0;];
tau_reg  = t9;
