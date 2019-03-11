% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPP2
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
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:00
% EndTime: 2019-03-08 18:19:01
% DurationCPUTime: 0.20s
% Computational Cost: add. (165->66), mult. (300->80), div. (0->0), fcn. (219->6), ass. (0->40)
t35 = sin(qJ(2));
t36 = cos(qJ(2));
t43 = qJD(1) * qJD(2);
t49 = t35 * qJDD(1) + t36 * t43;
t46 = t36 * qJD(1);
t21 = qJD(2) * pkin(2) + t46;
t34 = cos(pkin(5));
t48 = qJD(1) * t35;
t24 = t34 * t48;
t33 = sin(pkin(5));
t8 = t33 * t21 + t24;
t47 = qJDD(2) * pkin(3);
t45 = qJDD(1) - g(2);
t29 = t36 * qJDD(1);
t15 = qJDD(2) * pkin(2) - t35 * t43 + t29;
t4 = t33 * t15 + t49 * t34;
t3 = t34 * t15 - t49 * t33;
t41 = qJDD(2) * qJ(4) + t4;
t40 = -qJDD(4) + t3;
t39 = g(1) * t35 - g(2) * t36;
t17 = t33 * t36 + t34 * t35;
t16 = t33 * t35 - t34 * t36;
t7 = t34 * t21 - t33 * t48;
t30 = qJ(2) + pkin(5);
t27 = sin(t30);
t28 = cos(t30);
t38 = g(1) * t27 - g(2) * t28 + t40;
t37 = qJD(2) ^ 2;
t32 = qJDD(3) - g(3);
t26 = -t34 * pkin(2) - pkin(3);
t25 = t33 * pkin(2) + qJ(4);
t14 = t16 * qJD(1);
t13 = t16 * qJD(2);
t12 = t33 * t46 + t24;
t11 = t17 * qJD(2);
t6 = qJD(2) * qJ(4) + t8;
t5 = -qJD(2) * pkin(3) + qJD(4) - t7;
t2 = -t40 - t47;
t1 = qJD(4) * qJD(2) + t41;
t9 = [t45, 0, t36 * qJDD(2) - t37 * t35, -qJDD(2) * t35 - t37 * t36, -t7 * t11 - t8 * t13 - t3 * t16 + t4 * t17 - g(2), -t11 * qJD(2) - t16 * qJDD(2), -t13 * qJD(2) + t17 * qJDD(2), t1 * t17 + t5 * t11 - t6 * t13 + t2 * t16 - g(2); 0, qJDD(2), t29 + t39, g(1) * t36 - t45 * t35, t7 * t12 + t8 * t14 + (t3 * t34 + t33 * t4 + t39) * pkin(2), t12 * qJD(2) + (pkin(3) - t26) * qJDD(2) + t38, -g(1) * t28 - g(2) * t27 + t25 * qJDD(2) + (0.2e1 * qJD(4) + t14) * qJD(2) + t41, t1 * t25 + t2 * t26 - t5 * t12 - g(1) * (-t35 * pkin(2) - t27 * pkin(3) + t28 * qJ(4)) - g(2) * (t36 * pkin(2) + t28 * pkin(3) + t27 * qJ(4)) + (qJD(4) + t14) * t6; 0, 0, 0, 0, t32, 0, 0, t32; 0, 0, 0, 0, 0, -qJDD(2), -t37, -t6 * qJD(2) - t38 - t47;];
tau_reg  = t9;
