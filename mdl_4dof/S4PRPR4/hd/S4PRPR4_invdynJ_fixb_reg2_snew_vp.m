% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (230->51), mult. (423->68), div. (0->0), fcn. (254->6), ass. (0->48)
t52 = pkin(5) + pkin(2);
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t15 = t26 * g(1) - t27 * g(2);
t16 = -t27 * g(1) - t26 * g(2);
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t51 = -t29 * t15 - t31 * t16;
t28 = sin(qJ(4));
t23 = t28 ^ 2;
t33 = qJD(2) ^ 2;
t50 = t23 * t33;
t30 = cos(qJ(4));
t24 = t30 ^ 2;
t49 = t24 * t33;
t39 = t30 * t33 * t28;
t17 = qJDD(4) + t39;
t48 = t28 * t17;
t18 = qJDD(4) - t39;
t47 = t30 * t18;
t46 = t23 + t24;
t45 = t33 * qJ(3);
t44 = qJDD(2) * pkin(2);
t43 = t28 * qJDD(2);
t42 = qJD(2) * qJD(4);
t41 = qJDD(2) * qJ(3);
t40 = (2 * qJD(3) * qJD(2)) - t51;
t38 = t31 * t15 - t29 * t16;
t25 = -g(3) + qJDD(1);
t37 = qJDD(3) - t38;
t35 = t37 - t45;
t34 = -t52 * qJDD(2) + t35;
t2 = t28 * t25 - t30 * t34;
t3 = t30 * t25 + t28 * t34;
t1 = -t30 * t2 + t28 * t3;
t21 = t30 * qJDD(2);
t13 = -0.2e1 * t28 * t42 + t21;
t36 = t40 + t41;
t12 = 0.2e1 * t30 * t42 + t43;
t32 = qJD(4) ^ 2;
t20 = -t32 - t49;
t19 = -t32 - t50;
t14 = t46 * qJDD(2);
t8 = t30 * t20 - t48;
t7 = t28 * t19 + t47;
t6 = t35 - t44;
t5 = -t52 * t33 + t36;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, -t28 * t18 + t30 * t19, -t30 * t17 - t28 * t20, 0, t28 * t2 + t30 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t38, t51, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, t37 - 0.2e1 * t44, t40 + 0.2e1 * t41, -pkin(2) * t6 + qJ(3) * (-t33 * pkin(2) + t36), t13 * t30, -t30 * t12 - t28 * t13, t47 - t28 * (t32 - t49), t12 * t28, t30 * (-t32 + t50) - t48, 0, qJ(3) * t12 + t28 * t5 - t52 * t7, qJ(3) * t13 + t30 * t5 - t52 * t8, t52 * t14 - t46 * t45 - t1, qJ(3) * t5 - t52 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t33, t6, 0, 0, 0, 0, 0, 0, t7, t8, -t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, (-t23 + t24) * t33, t21, -t39, -t43, qJDD(4), -t2, -t3, 0, 0;];
tauJ_reg = t4;
