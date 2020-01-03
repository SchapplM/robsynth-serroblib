% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tauc_reg [4x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (139->53), mult. (278->69), div. (0->0), fcn. (107->4), ass. (0->44)
t15 = qJD(1) + qJD(2);
t14 = t15 ^ 2;
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t39 = qJD(4) * t20;
t21 = cos(qJ(2));
t42 = pkin(1) * qJD(2);
t30 = qJD(1) * t42;
t29 = t21 * t30;
t38 = t15 * qJD(3);
t7 = t29 + t38;
t19 = sin(qJ(2));
t43 = pkin(1) * qJD(1);
t33 = t19 * t43;
t41 = t15 * qJ(3);
t9 = t33 + t41;
t48 = t7 * t18 + t9 * t39;
t23 = qJD(4) ^ 2;
t47 = t23 * t18;
t46 = t23 * t20;
t45 = -t14 - t23;
t44 = t18 ^ 2 - t20 ^ 2;
t40 = qJD(4) * t18;
t37 = qJD(1) + t15;
t36 = qJD(2) - t15;
t35 = t19 * t42;
t34 = t21 * t42;
t32 = t21 * t43;
t31 = -t21 * pkin(1) - pkin(2);
t28 = t36 * t43;
t27 = qJD(3) - t32;
t11 = qJD(3) + t34;
t26 = t11 * t15 - (-pkin(6) + t31) * t23;
t13 = t19 * pkin(1) + qJ(3);
t25 = t13 * t15 + t35;
t24 = -t9 * t15 + t19 * t30;
t22 = -pkin(2) - pkin(6);
t10 = -0.2e1 * t15 * t18 * t39;
t8 = -t15 * pkin(2) + t27;
t5 = t7 * t20;
t3 = t37 * t35;
t2 = t19 * t28;
t1 = 0.2e1 * t44 * t15 * qJD(4);
t4 = [0, 0, 0, 0, -t3, -t37 * t34, t3, t29 + (qJD(3) + t11) * t15, t9 * t11 + t7 * t13 + (t31 * qJD(1) + t8) * t35, t10, t1, -t47, -t46, 0, t26 * t18 + t25 * t39 + t48, t5 + t26 * t20 + (-t25 - t9) * t40; 0, 0, 0, 0, -t2, -t36 * t32, t2, t21 * t28 + 0.2e1 * t38, t7 * qJ(3) + t9 * qJD(3) + (-t21 * t9 + (-pkin(2) * qJD(2) - t8) * t19) * t43, t10, t1, -t47, -t46, 0, -t33 * t39 - t22 * t47 + (qJ(3) * t39 + t27 * t18) * t15 + t48, t5 + (t27 * t15 - t22 * t23) * t20 + (t33 - t9 - t41) * t40; 0, 0, 0, 0, 0, 0, 0, -t14, t24, 0, 0, 0, 0, 0, t45 * t18, t45 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t14 * t18, -t44 * t14, 0, 0, 0, t24 * t20, -t24 * t18;];
tauc_reg = t4;
