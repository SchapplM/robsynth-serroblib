% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.24s
% Computational Cost: add. (193->50), mult. (455->80), div. (0->0), fcn. (292->6), ass. (0->39)
t62 = 2 * qJD(3);
t29 = sin(pkin(7)) * pkin(1) + qJ(3);
t61 = qJD(1) * t29;
t35 = sin(pkin(8));
t37 = cos(pkin(8));
t52 = t35 ^ 2 + t37 ^ 2;
t39 = sin(qJ(5));
t40 = cos(qJ(5));
t53 = t40 * t37;
t60 = -t39 * t35 + t53;
t28 = -cos(pkin(7)) * pkin(1) - pkin(2) - qJ(4);
t59 = t28 * qJD(1);
t58 = t52 * qJD(4);
t57 = 0.2e1 * qJD(1);
t55 = -pkin(6) + t28;
t20 = t40 * t35 + t39 * t37;
t51 = qJD(1) * t20;
t49 = qJD(1) * t35;
t18 = t60 * qJD(5);
t8 = t18 * qJD(5);
t48 = t39 * t49;
t47 = qJD(1) * t53;
t46 = qJD(3) * t57;
t22 = qJD(4) + t61;
t45 = t52 * (qJD(3) + t59);
t17 = t20 * qJD(5);
t43 = t20 * qJD(4);
t42 = t60 * qJD(4);
t41 = qJD(1) ^ 2;
t25 = qJD(5) * t48;
t23 = t35 * pkin(4) + t29;
t16 = t47 - t48;
t13 = t55 * t37;
t12 = t55 * t35;
t11 = pkin(4) * t49 + t22;
t7 = t17 * qJD(5);
t6 = qJD(5) * t47 - t25;
t5 = qJD(1) * t17;
t1 = [0, 0, 0, 0, 0, t46, t61 * t62, t35 * t46, t37 * t46, t57 * t58, (t22 + t61) * qJD(3) + (-t52 * t59 - t45) * qJD(4), -t16 * t17 - t5 * t60, -t16 * t18 + t17 * t51 + t5 * t20 - t6 * t60, -t7, -t8, 0, t11 * t18 + t23 * t6 + ((-t12 * t40 - t13 * t39) * qJD(5) - t42) * qJD(5) + t51 * t62, -t11 * t17 - t23 * t5 + ((t12 * t39 - t13 * t40) * qJD(5) + t43) * qJD(5) + (qJD(1) * t60 + t16) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7; 0, 0, 0, 0, 0, -t41, -t61 * qJD(1), -t41 * t35, -t41 * t37, 0, (-t22 - t58) * qJD(1), 0, 0, 0, 0, 0, -qJD(1) * t51 - t7, -qJD(1) * t16 - t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t41, (qJD(3) + t45) * qJD(1), 0, 0, 0, 0, 0, -t25 + (t16 + t47) * qJD(5), -0.2e1 * t51 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t51, t16 ^ 2 - t51 ^ 2, 0, t25 + (t16 - t47) * qJD(5), 0, -qJD(1) * t42 - t11 * t16, qJD(1) * t43 + t11 * t51;];
tauc_reg = t1;
