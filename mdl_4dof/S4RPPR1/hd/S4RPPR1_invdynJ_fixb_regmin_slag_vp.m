% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPPR1
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
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4RPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:46
% EndTime: 2018-11-14 13:46:47
% DurationCPUTime: 0.18s
% Computational Cost: add. (184->62), mult. (229->75), div. (0->0), fcn. (127->8), ass. (0->36)
t21 = sin(pkin(6));
t11 = t21 * pkin(1) + qJ(3);
t39 = t11 * qJDD(1);
t18 = qJ(1) + pkin(6);
t14 = sin(t18);
t15 = cos(t18);
t45 = -g(1) * t14 + g(2) * t15;
t37 = qJD(1) - qJD(4);
t44 = qJD(4) + t37;
t22 = cos(pkin(6));
t12 = -t22 * pkin(1) - pkin(2);
t43 = qJDD(1) * t12;
t38 = qJD(3) * qJD(1);
t10 = -pkin(3) + t12;
t16 = qJDD(1) - qJDD(4);
t3 = qJDD(1) * t10 + qJDD(3);
t36 = t10 * t16 + t3;
t35 = t37 ^ 2;
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t1 = -t14 * t23 - t15 * t25;
t2 = -t14 * t25 + t15 * t23;
t34 = -g(1) * t2 + g(2) * t1;
t33 = g(1) * t1 + g(2) * t2;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t32 = g(1) * t24 - g(2) * t26;
t4 = qJD(1) * t10 + qJD(3);
t6 = t11 * qJD(1);
t31 = t23 * t6 - t25 * t4;
t30 = -t23 * t4 - t25 * t6;
t5 = t38 + t39;
t29 = qJD(3) * t37 + t11 * t16 + t5;
t28 = qJDD(3) + t43;
t20 = qJDD(2) - g(3);
t7 = [qJDD(1), t32, g(1) * t26 + g(2) * t24 (t32 + (t21 ^ 2 + t22 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), -qJDD(3) - 0.2e1 * t43 - t45, -g(1) * t15 - g(2) * t14 + 0.2e1 * t38 + 0.2e1 * t39, t5 * t11 + t6 * qJD(3) + t28 * t12 - g(1) * (-t24 * pkin(1) - t14 * pkin(2) + t15 * qJ(3)) - g(2) * (t26 * pkin(1) + t15 * pkin(2) + t14 * qJ(3)) t16, -t36 * t25 + t29 * t23 + (-(-t23 * t10 - t25 * t11) * t37 - t30) * qJD(4) + t34, t36 * t23 + t29 * t25 + ((t25 * t10 - t23 * t11) * t37 - t31) * qJD(4) + t33; 0, 0, 0, t20, 0, 0, t20, 0, 0, 0; 0, 0, 0, 0, -qJDD(1), -qJD(1) ^ 2, -t6 * qJD(1) + t28 + t45, 0, -t25 * t16 - t23 * t35, t23 * t16 - t25 * t35; 0, 0, 0, 0, 0, 0, 0, -t16, -t23 * t5 + t25 * t3 + t44 * t30 - t34, -t23 * t3 - t25 * t5 + t44 * t31 - t33;];
tau_reg  = t7;
