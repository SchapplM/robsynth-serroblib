% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% tauc_reg [4x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.13s
% Computational Cost: add. (96->32), mult. (215->56), div. (0->0), fcn. (98->4), ass. (0->29)
t14 = sin(pkin(6));
t19 = qJD(4) ^ 2;
t20 = qJD(1) ^ 2;
t40 = t14 * (t19 + t20);
t39 = 2 * qJD(1);
t18 = (-pkin(1) - pkin(2));
t15 = cos(pkin(6));
t8 = t18 * qJD(1) + qJD(2);
t38 = t15 * t8;
t16 = sin(qJ(4));
t37 = t19 * t16;
t17 = cos(qJ(4));
t36 = t19 * t17;
t35 = t15 * qJ(2) + t14 * t18;
t34 = t16 ^ 2 - t17 ^ 2;
t32 = t14 * qJ(2);
t31 = qJD(2) * t15;
t30 = qJD(1) * qJ(2);
t29 = qJD(2) * t39;
t28 = qJD(4) * t39;
t27 = t17 * t28;
t26 = t14 * t29;
t25 = (-t14 * t30 + t38) * t14 - (t14 * t8 + t15 * t30) * t15;
t24 = t15 * t18 - t32;
t1 = -t38 + (pkin(3) + t32) * qJD(1);
t23 = (t1 - t31) * qJD(1);
t22 = -t19 * (-pkin(5) + t35) + t26;
t21 = qJD(4) * (-qJD(1) * (pkin(3) - t24) - t1 - t31);
t2 = [0, 0, 0, 0, t29, qJ(2) * t29, t26, t15 * t29, ((-t14 * t24 + t15 * t35) * qJD(1) - t25) * qJD(2), t16 * t27, -t34 * t28, -t36, t37, 0, t16 * t21 + t22 * t17, -t22 * t16 + t17 * t21; 0, 0, 0, 0, -t20, -t20 * qJ(2), -t14 * t20, -t15 * t20, t25 * qJD(1), 0, 0, 0, 0, 0, t15 * t16 * t28 - t17 * t40, t15 * t27 + t16 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t20 * t17, t34 * t20, 0, 0, 0, t16 * t23, t17 * t23;];
tauc_reg = t2;
