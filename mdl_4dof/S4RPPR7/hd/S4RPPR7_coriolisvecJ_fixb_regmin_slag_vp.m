% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR7
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
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:41
% EndTime: 2019-12-31 16:41:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (148->47), mult. (382->77), div. (0->0), fcn. (238->4), ass. (0->37)
t25 = sin(pkin(6));
t28 = sin(qJ(4));
t26 = cos(pkin(6));
t29 = cos(qJ(4));
t41 = t29 * t26;
t47 = -t28 * t25 + t41;
t27 = -pkin(1) - qJ(3);
t46 = (t27 * qJD(1));
t40 = t25 ^ 2 + t26 ^ 2;
t45 = t40 * qJD(3);
t44 = 2 * qJD(1);
t43 = -pkin(5) + t27;
t11 = t29 * t25 + t28 * t26;
t39 = qJD(1) * t11;
t38 = qJD(1) * t25;
t10 = t47 * qJD(4);
t37 = t10 * qJD(4);
t24 = qJD(1) * qJ(2);
t21 = qJD(3) + t24;
t36 = t28 * t38;
t35 = qJD(1) * t41;
t34 = qJD(2) * t44;
t9 = t11 * qJD(4);
t32 = t11 * qJD(3);
t31 = t47 * qJD(3);
t30 = qJD(1) ^ 2;
t20 = t25 * pkin(3) + qJ(2);
t17 = qJD(2) + t46;
t16 = qJD(4) * t36;
t15 = pkin(3) * t38 + t21;
t14 = t43 * t26;
t13 = t43 * t25;
t8 = t35 - t36;
t5 = t9 * qJD(4);
t2 = qJD(4) * t35 - t16;
t1 = qJD(1) * t9;
t3 = [0, 0, 0, 0, t34, qJ(2) * t34, t25 * t34, t26 * t34, t44 * t45, (t21 + t24) * qJD(2) + (-t17 - t46) * t45, -t1 * t47 - t8 * t9, t1 * t11 - t8 * t10 - t2 * t47 + t39 * t9, -t5, -t37, 0, t15 * t10 + t20 * t2 + ((-t13 * t29 - t14 * t28) * qJD(4) - t31) * qJD(4) + 0.2e1 * t39 * qJD(2), -t20 * t1 - t15 * t9 + ((t13 * t28 - t14 * t29) * qJD(4) + t32) * qJD(4) + (qJD(1) * t47 + t8) * qJD(2); 0, 0, 0, 0, -t30, -t30 * qJ(2), -t30 * t25, -t30 * t26, 0, (-t21 - t45) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t39 - t5, -qJD(1) * t8 - t37; 0, 0, 0, 0, 0, 0, 0, 0, -t40 * t30, (t40 * t17 + qJD(2)) * qJD(1), 0, 0, 0, 0, 0, -t16 + (t8 + t35) * qJD(4), -0.2e1 * t39 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t39, -t39 ^ 2 + t8 ^ 2, 0, t16 + (t8 - t35) * qJD(4), 0, -qJD(1) * t31 - t15 * t8, qJD(1) * t32 + t15 * t39;];
tauc_reg = t3;
