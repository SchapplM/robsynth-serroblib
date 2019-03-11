% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP2
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
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:10
% DurationCPUTime: 0.16s
% Computational Cost: add. (191->59), mult. (343->73), div. (0->0), fcn. (243->6), ass. (0->40)
t25 = qJ(2) + qJ(3);
t21 = sin(t25);
t22 = cos(t25);
t45 = g(1) * t22 + g(2) * t21;
t44 = g(1) * t21 - g(2) * t22;
t24 = qJD(2) + qJD(3);
t29 = cos(qJ(2));
t20 = t29 * qJDD(1);
t27 = sin(qJ(2));
t39 = qJD(1) * qJD(2);
t11 = qJDD(2) * pkin(2) - t27 * t39 + t20;
t26 = sin(qJ(3));
t41 = qJD(1) * t27;
t38 = qJD(3) * t41;
t15 = t26 * t38;
t28 = cos(qJ(3));
t16 = qJD(2) * pkin(2) + t29 * qJD(1);
t34 = -t27 * qJDD(1) - t29 * t39;
t33 = qJD(3) * t16 - t34;
t2 = t26 * t11 + t33 * t28 - t15;
t23 = qJDD(2) + qJDD(3);
t43 = pkin(2) * t23;
t40 = qJDD(1) - g(2);
t36 = g(1) * t27 - g(2) * t29;
t6 = t28 * t16 - t26 * t41;
t13 = t26 * t29 + t28 * t27;
t12 = -t26 * t27 + t28 * t29;
t8 = t28 * t11;
t31 = -t33 * t26 - t28 * t38 + t8;
t1 = t23 * pkin(3) + t31;
t35 = (t1 + t44) * pkin(3);
t32 = (-pkin(2) * t24 - t16) * qJD(3) + t34;
t30 = qJD(2) ^ 2;
t10 = t12 * qJD(1);
t9 = t13 * qJD(1);
t7 = t26 * t16 + t28 * t41;
t5 = t24 * pkin(3) + t6;
t4 = t24 * t13;
t3 = t24 * t12;
t14 = [t40, 0, t29 * qJDD(2) - t30 * t27, -qJDD(2) * t27 - t30 * t29, 0, t12 * t23 - t4 * t24, -t13 * t23 - t3 * t24, t1 * t12 + t2 * t13 + t7 * t3 - t5 * t4 - g(2); 0, qJDD(2), t20 + t36, g(1) * t29 - t40 * t27, t23, t9 * t24 + t8 + (-t38 + t43) * t28 + t32 * t26 + t44, t10 * t24 + t15 + (-t11 - t43) * t26 + t32 * t28 + t45, -t7 * t10 + t5 * t9 + t35 + (t1 * t28 + t2 * t26 + (-t26 * t5 + t28 * t7) * qJD(3) + t36) * pkin(2); 0, 0, 0, 0, t23, t7 * t24 + t31 + t44, t6 * t24 - t2 + t45 (t5 - t6) * t7 + t35; 0, 0, 0, 0, 0, 0, 0, qJDD(4) - g(3);];
tau_reg  = t14;
