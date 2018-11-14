% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRPP4
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
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:31
% EndTime: 2018-11-14 14:09:31
% DurationCPUTime: 0.19s
% Computational Cost: add. (195->73), mult. (358->86), div. (0->0), fcn. (265->6), ass. (0->41)
t37 = sin(qJ(2));
t38 = cos(qJ(2));
t46 = qJD(1) * qJD(2);
t53 = qJDD(1) * t37 + t38 * t46;
t32 = qJ(2) + pkin(5);
t29 = sin(t32);
t30 = cos(t32);
t52 = -g(1) * t30 + g(2) * t29;
t50 = qJD(1) * t37;
t49 = qJDD(2) * pkin(2);
t48 = qJDD(1) - g(1);
t31 = t38 * qJDD(1);
t16 = -t37 * t46 + t31 + t49;
t35 = sin(pkin(5));
t36 = cos(pkin(5));
t45 = -t36 * t16 + t53 * t35;
t4 = t35 * t16 + t53 * t36;
t43 = -g(1) * t38 + g(2) * t37;
t22 = qJD(2) * pkin(2) + qJD(1) * t38;
t9 = t35 * t22 + t36 * t50;
t18 = t35 * t38 + t36 * t37;
t17 = t35 * t37 - t36 * t38;
t12 = t18 * qJD(2);
t42 = -qJD(2) * t12 - qJDD(2) * t17;
t8 = t22 * t36 - t35 * t50;
t2 = -qJDD(2) * pkin(3) + qJDD(4) + t45;
t41 = -g(1) * t29 - g(2) * t30 + t4;
t13 = t18 * qJD(1);
t40 = qJD(2) * t13 - t45 + t52;
t39 = qJD(2) ^ 2;
t34 = qJDD(3) + g(3);
t33 = qJDD(2) * qJ(4);
t28 = -pkin(2) * t36 - pkin(3);
t25 = pkin(2) * t35 + qJ(4);
t15 = t17 * qJD(1);
t14 = t17 * qJD(2);
t7 = qJD(2) * qJ(4) + t9;
t6 = -qJD(2) * pkin(3) + qJD(4) - t8;
t5 = -qJD(2) * t14 + qJDD(2) * t18;
t1 = qJD(4) * qJD(2) + t33 + t4;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, qJDD(2) * t38 - t37 * t39, -qJDD(2) * t37 - t38 * t39, 0, -g(1) + (t37 ^ 2 + t38 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t42, -t5, 0, -t12 * t8 - t14 * t9 + t17 * t45 + t18 * t4 - g(1), 0, 0, 0, 0, 0, 0, t42, 0, t5, t1 * t18 + t12 * t6 - t14 * t7 + t17 * t2 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t31 + t43, g(2) * t38 - t48 * t37, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t36 * t49 + t40, -qJD(2) * t15 - t35 * t49 - t41, 0, t13 * t8 + t15 * t9 + (t35 * t4 - t36 * t45 + t43) * pkin(2), 0, 0, 0, qJDD(2), 0, 0, -qJDD(4) + (pkin(3) - t28) * qJDD(2) + t40, 0, qJDD(2) * t25 + t33 + (0.2e1 * qJD(4) + t15) * qJD(2) + t41, t1 * t25 + t2 * t28 - t6 * t13 - g(1) * (pkin(2) * t38 + pkin(3) * t30 + qJ(4) * t29) - g(2) * (-pkin(2) * t37 - pkin(3) * t29 + qJ(4) * t30) + (qJD(4) + t15) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t39, -qJD(2) * t7 + t2 - t52;];
tau_reg  = t3;
