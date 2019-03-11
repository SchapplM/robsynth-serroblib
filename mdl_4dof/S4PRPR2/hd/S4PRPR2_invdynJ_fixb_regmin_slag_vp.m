% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x8]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:54
% EndTime: 2019-03-08 18:21:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (209->65), mult. (392->95), div. (0->0), fcn. (310->8), ass. (0->43)
t28 = qJD(2) + qJD(4);
t49 = qJD(4) - t28;
t29 = sin(pkin(6));
t48 = pkin(2) * t29;
t32 = sin(qJ(2));
t47 = qJD(1) * t32;
t46 = qJDD(1) - g(2);
t45 = qJD(1) * qJD(2);
t26 = qJ(2) + pkin(6) + qJ(4);
t22 = sin(t26);
t23 = cos(t26);
t31 = sin(qJ(4));
t34 = cos(qJ(2));
t18 = qJD(2) * pkin(2) + t34 * qJD(1);
t30 = cos(pkin(6));
t8 = t29 * t18 + t30 * t47;
t44 = qJD(4) * t31 * t8 + g(1) * t23 + g(2) * t22;
t43 = g(1) * t32 - g(2) * t34;
t7 = t30 * t18 - t29 * t47;
t33 = cos(qJ(4));
t6 = qJD(2) * pkin(3) + t7;
t42 = -t31 * t6 - t33 * t8;
t15 = -t29 * t32 + t30 * t34;
t16 = t29 * t34 + t30 * t32;
t41 = t33 * t15 - t31 * t16;
t40 = t31 * t15 + t33 * t16;
t24 = t30 * pkin(2) + pkin(3);
t39 = t31 * t24 + t33 * t48;
t38 = t33 * t24 - t31 * t48;
t25 = t34 * qJDD(1);
t14 = qJDD(2) * pkin(2) - t32 * t45 + t25;
t36 = t32 * qJDD(1) + t34 * t45;
t3 = t30 * t14 - t29 * t36;
t2 = qJDD(2) * pkin(3) + t3;
t4 = t29 * t14 + t30 * t36;
t37 = g(1) * t22 - g(2) * t23 + t33 * t2 - t31 * t4;
t35 = qJD(2) ^ 2;
t27 = qJDD(2) + qJDD(4);
t13 = t15 * qJD(1);
t12 = t15 * qJD(2);
t11 = t16 * qJD(1);
t10 = t16 * qJD(2);
t1 = [t46, 0, t34 * qJDD(2) - t35 * t32, -qJDD(2) * t32 - t35 * t34, -t7 * t10 + t8 * t12 + t3 * t15 + t4 * t16 - g(2), 0 (-t40 * qJD(4) - t33 * t10 - t31 * t12) * t28 + t41 * t27 -(t41 * qJD(4) - t31 * t10 + t33 * t12) * t28 - t40 * t27; 0, qJDD(2), t25 + t43, g(1) * t34 - t46 * t32, t7 * t11 - t8 * t13 + (t29 * t4 + t3 * t30 + t43) * pkin(2), t27, t38 * t27 - (-t33 * t11 - t31 * t13) * t28 + (-t39 * t28 + t42) * qJD(4) + t37, -t39 * t27 - t33 * t4 - t31 * t2 + (-t31 * t11 + t33 * t13) * t28 + (-t38 * t28 - t33 * t6) * qJD(4) + t44; 0, 0, 0, 0, qJDD(3) - g(3), 0, 0, 0; 0, 0, 0, 0, 0, t27, t49 * t42 + t37 (-t8 * t28 - t2) * t31 + (-t49 * t6 - t4) * t33 + t44;];
tau_reg  = t1;
