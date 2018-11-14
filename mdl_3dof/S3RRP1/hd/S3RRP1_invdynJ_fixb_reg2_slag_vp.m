% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% tau_reg [3x(3*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S3RRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:09
% EndTime: 2018-11-14 10:15:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->57), mult. (200->69), div. (0->0), fcn. (86->6), ass. (0->39)
t25 = qJ(1) + qJ(2);
t21 = sin(t25);
t22 = cos(t25);
t49 = -g(1) * t22 - g(2) * t21;
t43 = pkin(1) * qJDD(1);
t23 = qJDD(1) + qJDD(2);
t48 = t23 * pkin(2);
t47 = t22 * pkin(2) + t21 * qJ(3);
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t44 = pkin(1) * qJD(1);
t37 = qJD(2) * t44;
t46 = t26 * t43 + t28 * t37;
t45 = t26 * t37 - t28 * t43;
t42 = qJD(2) * t26;
t41 = qJD(2) * t28;
t40 = t26 * t44;
t39 = t28 * t44;
t24 = qJD(1) + qJD(2);
t38 = t24 * t42;
t36 = -t21 * pkin(2) + t22 * qJ(3);
t35 = t46 + t49;
t19 = t23 * qJ(3);
t20 = t24 * qJD(3);
t1 = t19 + t20 + t46;
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t34 = g(1) * t27 - g(2) * t29;
t33 = g(1) * t21 - g(2) * t22 - t45;
t32 = -qJDD(3) + t33;
t31 = t24 * t39 - t35;
t16 = -t28 * pkin(1) - pkin(2);
t11 = t26 * pkin(1) + qJ(3);
t6 = pkin(1) * t41 + qJD(3);
t5 = t24 * t40;
t4 = t24 * qJ(3) + t40;
t3 = -t24 * pkin(2) + qJD(3) - t39;
t2 = qJDD(3) + t45 - t48;
t7 = [0, 0, 0, 0, 0, qJDD(1), t34, g(1) * t29 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, t23 (t23 * t28 - t38) * pkin(1) + t33 (-t23 * t26 - t24 * t41) * pkin(1) - t35, 0 (t34 + (t26 ^ 2 + t28 ^ 2) * t43) * pkin(1), 0, 0, 0, t23, 0, 0, -pkin(1) * t38 + (pkin(2) - t16) * t23 + t32, 0, t11 * t23 + t6 * t24 + t1 + t49, t1 * t11 + t4 * t6 + t2 * t16 + t3 * pkin(1) * t42 - g(1) * (-t27 * pkin(1) + t36) - g(2) * (t29 * pkin(1) + t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t33 + t5, t31, 0, 0, 0, 0, 0, t23, 0, 0, t32 + t5 + 0.2e1 * t48, 0, 0.2e1 * t19 + 0.2e1 * t20 - t31, t1 * qJ(3) + t4 * qJD(3) - t2 * pkin(2) - g(1) * t36 - g(2) * t47 + (-t26 * t3 - t28 * t4) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24 ^ 2, -t4 * t24 - t32 - t48;];
tau_reg  = t7;
