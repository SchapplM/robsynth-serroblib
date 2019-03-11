% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PPRR1
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:59
% EndTime: 2019-03-08 18:16:00
% DurationCPUTime: 0.21s
% Computational Cost: add. (190->61), mult. (352->82), div. (0->0), fcn. (278->8), ass. (0->46)
t27 = sin(pkin(6));
t28 = cos(pkin(6));
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t13 = -t27 * t30 - t28 * t32;
t38 = t27 * t32 - t28 * t30;
t52 = -g(1) * t38 - g(2) * t13;
t26 = qJ(3) + qJ(4);
t21 = sin(t26);
t47 = cos(t26);
t8 = -t27 * t21 - t28 * t47;
t9 = t28 * t21 - t27 * t47;
t51 = g(1) * t9 - g(2) * t8;
t24 = qJD(3) + qJD(4);
t20 = t32 * qJDD(2);
t41 = qJD(2) * qJD(3);
t12 = qJDD(3) * pkin(3) - t30 * t41 + t20;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t18 = qJD(3) * pkin(3) + t32 * qJD(2);
t42 = t30 * qJDD(2);
t36 = -t32 * t41 - t42;
t35 = qJD(4) * t18 - t36;
t50 = t29 * t12 + t31 * t35;
t23 = qJDD(3) + qJDD(4);
t46 = pkin(3) * t23;
t43 = qJD(2) * t30;
t40 = -g(1) * t27 + g(2) * t28;
t39 = qJD(4) * t43;
t16 = t29 * t32 + t31 * t30;
t15 = -t29 * t30 + t31 * t32;
t17 = t29 * t39;
t37 = -g(1) * t8 - g(2) * t9 + t17;
t34 = (-pkin(3) * t24 - t18) * qJD(4) + t36;
t7 = t31 * t12;
t2 = -t29 * t35 - t31 * t39 + t7;
t33 = qJD(3) ^ 2;
t25 = qJDD(1) - g(3);
t11 = t15 * qJD(2);
t10 = t16 * qJD(2);
t6 = t29 * t18 + t31 * t43;
t5 = t31 * t18 - t29 * t43;
t4 = t24 * t16;
t3 = t24 * t15;
t1 = -t17 + t50;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t40, 0, 0, 0, 0, 0, 0, t32 * qJDD(3) - t33 * t30, -qJDD(3) * t30 - t33 * t32, 0 (t30 ^ 2 + t32 ^ 2) * qJDD(2) + t40, 0, 0, 0, 0, 0, 0, t15 * t23 - t4 * t24, -t16 * t23 - t3 * t24, 0, t1 * t16 + t2 * t15 + t6 * t3 - t5 * t4 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t20 + t52, -g(1) * t13 + g(2) * t38 - t42, 0, 0, 0, 0, 0, 0, 0, t23, t10 * t24 + t7 + (-t39 + t46) * t31 + t34 * t29 + t51, t11 * t24 + (-t12 - t46) * t29 + t34 * t31 + t37, 0, t5 * t10 - t6 * t11 + (t1 * t29 + t2 * t31 + (-t5 * t29 + t6 * t31) * qJD(4) + t52) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t6 * t24 + t2 + t51, t5 * t24 + t37 - t50, 0, 0;];
tau_reg  = t14;
