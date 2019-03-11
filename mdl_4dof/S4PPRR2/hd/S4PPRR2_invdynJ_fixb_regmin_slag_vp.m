% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:10
% EndTime: 2019-03-08 18:17:10
% DurationCPUTime: 0.16s
% Computational Cost: add. (162->51), mult. (314->65), div. (0->0), fcn. (270->10), ass. (0->35)
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t46 = -t28 * t25 + t30 * t26;
t6 = t46 * qJD(1);
t11 = t30 * t25 + t28 * t26;
t7 = t11 * qJD(1);
t24 = qJD(3) + qJD(4);
t45 = qJD(4) - t24;
t44 = t11 * qJDD(1);
t22 = qJDD(3) + qJDD(4);
t43 = pkin(3) * t22;
t29 = cos(qJ(4));
t42 = t29 * t7;
t23 = pkin(6) + qJ(3);
t21 = qJ(4) + t23;
t17 = sin(t21);
t18 = cos(t21);
t27 = sin(qJ(4));
t39 = qJD(4) * t27 * t7 + g(1) * t18 + g(2) * t17;
t5 = qJD(3) * pkin(3) + t6;
t38 = -pkin(3) * t24 - t5;
t35 = t46 * qJDD(1);
t9 = t11 * qJD(3);
t2 = qJDD(3) * pkin(3) - qJD(1) * t9 + t35;
t37 = -t7 * t24 - t2;
t34 = -t27 * t11 + t29 * t46;
t33 = t29 * t11 + t27 * t46;
t8 = t46 * qJD(3);
t3 = qJD(1) * t8 + t44;
t31 = g(1) * t17 - g(2) * t18 + t29 * t2 - t27 * t3;
t20 = cos(t23);
t19 = sin(t23);
t1 = [qJDD(1) - g(2), -g(2) + (t25 ^ 2 + t26 ^ 2) * qJDD(1), 0, -t9 * qJD(3) + qJDD(3) * t46, -t8 * qJD(3) - t11 * qJDD(3), 0 (-qJD(4) * t33 - t27 * t8 - t29 * t9) * t24 + t34 * t22 -(qJD(4) * t34 - t27 * t9 + t29 * t8) * t24 - t33 * t22; 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0; 0, 0, qJDD(3), g(1) * t19 - g(2) * t20 + t35, g(1) * t20 + g(2) * t19 - t44, t22, t29 * t43 - (-t27 * t6 - t42) * t24 + (t38 * t27 - t42) * qJD(4) + t31 (t37 - t43) * t27 + (t38 * qJD(4) + t6 * t24 - t3) * t29 + t39; 0, 0, 0, 0, 0, t22, t31 + t45 * (-t27 * t5 - t42) t37 * t27 + (-t45 * t5 - t3) * t29 + t39;];
tau_reg  = t1;
