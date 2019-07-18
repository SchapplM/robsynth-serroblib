% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tauc_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.21s
% Computational Cost: add. (199->61), mult. (307->85), div. (0->0), fcn. (123->4), ass. (0->37)
t24 = -pkin(2) - pkin(3);
t19 = qJD(1) + qJD(2);
t21 = sin(qJ(2));
t42 = pkin(1) * qJD(1);
t32 = t21 * t42;
t11 = t19 * qJ(3) + t32;
t20 = sin(qJ(4));
t43 = t20 * t11;
t23 = cos(qJ(2));
t41 = pkin(1) * qJD(2);
t29 = qJD(1) * t41;
t14 = t23 * t29;
t17 = t19 * qJD(3);
t9 = t14 + t17;
t40 = qJD(2) * t21;
t39 = qJD(4) * t20;
t22 = cos(qJ(4));
t38 = qJD(4) * t22;
t37 = -qJD(4) + t19;
t31 = t23 * t42;
t26 = qJD(3) - t31;
t3 = t24 * t19 + t26;
t36 = t11 * t38 + t20 * t9 + t3 * t39;
t27 = t21 * t29;
t35 = t20 * t27 + t22 * t9 + t3 * t38;
t34 = pkin(1) * t40;
t33 = t23 * t41;
t30 = -t23 * pkin(1) - pkin(2);
t28 = t37 ^ 2;
t25 = t19 * t31 - t14;
t16 = t21 * pkin(1) + qJ(3);
t15 = -pkin(3) + t30;
t13 = qJD(3) + t33;
t10 = -t19 * pkin(2) + t26;
t5 = (-qJD(1) - t19) * t34;
t4 = (-qJD(2) + t19) * t32;
t1 = [0, 0, 0, 0, t5, -t19 * t33 - t14, t5, t13 * t19 + t9, t11 * t13 + t9 * t16 + (t30 * qJD(1) + t10) * t34, 0, -(-qJD(4) * t15 - t13) * t37 * t20 + (-(-qJD(4) * t16 + t34) * t37 - t27) * t22 + t36, (t22 * t13 + t20 * t34) * t37 + ((t15 * t22 - t16 * t20) * t37 - t43) * qJD(4) + t35; 0, 0, 0, 0, t4, t25, t4, 0.2e1 * t17 - t25, t9 * qJ(3) + t11 * qJD(3) + (-t11 * t23 + (-pkin(2) * qJD(2) - t10) * t21) * t42, 0, -(-t20 * qJD(3) + (-qJ(3) * t22 - t20 * t24) * qJD(4)) * t37 + (-t22 * t40 + (-t20 * t23 + t21 * t22) * t37) * t42 + t36, -t43 * qJD(4) + t35 - (-t22 * qJD(3) + (qJ(3) * t20 - t22 * t24) * qJD(4) + (t20 * t21 + t22 * t23) * t42) * t37; 0, 0, 0, 0, 0, 0, 0, -t19 ^ 2, -t11 * t19 + t27, 0, -t20 * t28, -t22 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t27 + (-t22 * t11 - t20 * t3) * t37 - t36, t11 * t39 - (t22 * t3 - t43) * t37 - t35;];
tauc_reg  = t1;
