% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
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
% tau_reg [3x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S3RRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:53
% EndTime: 2019-03-08 18:06:53
% DurationCPUTime: 0.13s
% Computational Cost: add. (158->55), mult. (188->65), div. (0->0), fcn. (82->6), ass. (0->38)
t25 = qJ(1) + qJ(2);
t21 = sin(t25);
t22 = cos(t25);
t47 = -g(1) * t22 - g(2) * t21;
t23 = qJDD(1) + qJDD(2);
t46 = t23 * pkin(2);
t45 = t22 * pkin(2) + t21 * qJ(3);
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t42 = pkin(1) * qJD(1);
t35 = qJD(2) * t42;
t41 = pkin(1) * qJDD(1);
t44 = t26 * t41 + t28 * t35;
t43 = t26 * t35 - t28 * t41;
t40 = qJD(2) * t26;
t39 = qJD(2) * t28;
t38 = t26 * t42;
t37 = t28 * t42;
t24 = qJD(1) + qJD(2);
t36 = t24 * t40;
t34 = -t21 * pkin(2) + t22 * qJ(3);
t33 = t44 + t47;
t19 = t23 * qJ(3);
t20 = t24 * qJD(3);
t1 = t19 + t20 + t44;
t32 = g(1) * t21 - g(2) * t22 - t43;
t31 = -qJDD(3) + t32;
t30 = t24 * t37 - t33;
t29 = cos(qJ(1));
t27 = sin(qJ(1));
t16 = -t28 * pkin(1) - pkin(2);
t11 = t26 * pkin(1) + qJ(3);
t6 = pkin(1) * t39 + qJD(3);
t5 = t24 * t38;
t4 = t24 * qJ(3) + t38;
t3 = -t24 * pkin(2) + qJD(3) - t37;
t2 = qJDD(3) + t43 - t46;
t7 = [qJDD(1), g(1) * t27 - g(2) * t29, g(1) * t29 + g(2) * t27, t23 (t23 * t28 - t36) * pkin(1) + t32 (-t23 * t26 - t24 * t39) * pkin(1) - t33, -pkin(1) * t36 + (pkin(2) - t16) * t23 + t31, t11 * t23 + t6 * t24 + t1 + t47, t1 * t11 + t4 * t6 + t2 * t16 + t3 * pkin(1) * t40 - g(1) * (-t27 * pkin(1) + t34) - g(2) * (t29 * pkin(1) + t45); 0, 0, 0, t23, t32 + t5, t30, t31 + t5 + 0.2e1 * t46, 0.2e1 * t19 + 0.2e1 * t20 - t30, t1 * qJ(3) + t4 * qJD(3) - t2 * pkin(2) - g(1) * t34 - g(2) * t45 + (-t26 * t3 - t28 * t4) * t42; 0, 0, 0, 0, 0, 0, -t23, -t24 ^ 2, -t4 * t24 - t31 - t46;];
tau_reg  = t7;
