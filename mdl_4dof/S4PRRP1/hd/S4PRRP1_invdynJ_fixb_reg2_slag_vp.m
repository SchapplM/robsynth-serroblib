% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRP1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:40
% EndTime: 2018-11-14 13:43:40
% DurationCPUTime: 0.17s
% Computational Cost: add. (201->59), mult. (200->69), div. (0->0), fcn. (86->6), ass. (0->41)
t27 = pkin(6) + qJ(2);
t25 = qJ(3) + t27;
t15 = sin(t25);
t16 = cos(t25);
t51 = -g(1) * t16 - g(2) * t15;
t45 = pkin(2) * qJDD(2);
t50 = t16 * pkin(3) + t15 * qJ(4);
t26 = qJDD(2) + qJDD(3);
t49 = pkin(3) * t26;
t30 = sin(qJ(3));
t31 = cos(qJ(3));
t46 = pkin(2) * qJD(2);
t38 = qJD(3) * t46;
t48 = t30 * t45 + t31 * t38;
t47 = -t30 * t38 + t31 * t45;
t44 = qJD(3) * t30;
t43 = qJD(3) * t31;
t42 = t30 * t46;
t41 = t31 * t46;
t28 = qJD(2) + qJD(3);
t40 = t28 * t44;
t39 = -pkin(3) * t15 + t16 * qJ(4);
t37 = t48 + t51;
t21 = t26 * qJ(4);
t22 = t28 * qJD(4);
t1 = t21 + t22 + t48;
t23 = sin(t27);
t24 = cos(t27);
t36 = g(1) * t23 - g(2) * t24;
t35 = g(1) * t15 - g(2) * t16 + t47;
t34 = -qJDD(4) + t35;
t33 = t28 * t41 - t37;
t29 = qJDD(1) - g(3);
t18 = -pkin(2) * t31 - pkin(3);
t17 = pkin(2) * t30 + qJ(4);
t8 = pkin(2) * t43 + qJD(4);
t5 = t28 * t42;
t4 = qJ(4) * t28 + t42;
t3 = -pkin(3) * t28 + qJD(4) - t41;
t2 = qJDD(4) - t47 - t49;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t36, g(1) * t24 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, t26 (t26 * t31 - t40) * pkin(2) + t35 (-t26 * t30 - t28 * t43) * pkin(2) - t37, 0 (t36 + (t30 ^ 2 + t31 ^ 2) * t45) * pkin(2), 0, 0, 0, t26, 0, 0, -pkin(2) * t40 + (pkin(3) - t18) * t26 + t34, 0, t17 * t26 + t28 * t8 + t1 + t51, t1 * t17 + t4 * t8 + t2 * t18 + t3 * pkin(2) * t44 - g(1) * (-pkin(2) * t23 + t39) - g(2) * (pkin(2) * t24 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t35 + t5, t33, 0, 0, 0, 0, 0, t26, 0, 0, t34 + t5 + 0.2e1 * t49, 0, 0.2e1 * t21 + 0.2e1 * t22 - t33, t1 * qJ(4) + t4 * qJD(4) - t2 * pkin(3) - g(1) * t39 - g(2) * t50 + (-t3 * t30 - t31 * t4) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t28 ^ 2, -t28 * t4 - t34 - t49;];
tau_reg  = t6;
