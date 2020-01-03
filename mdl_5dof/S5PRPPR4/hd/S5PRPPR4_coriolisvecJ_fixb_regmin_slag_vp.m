% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (162->55), mult. (460->91), div. (0->0), fcn. (303->4), ass. (0->47)
t34 = sin(pkin(8));
t32 = t34 ^ 2;
t35 = cos(pkin(8));
t57 = t35 ^ 2 + t32;
t64 = t57 * qJD(2) * qJD(3);
t37 = cos(qJ(5));
t36 = sin(qJ(5));
t59 = t35 * t36;
t18 = t34 * t37 - t59;
t45 = t34 * qJ(4) + pkin(2);
t63 = (t35 * pkin(3) + t45) * qJD(2);
t51 = qJ(3) * qJD(2);
t62 = (t34 * qJD(1) + t35 * t51) * t35;
t61 = t34 * t35;
t58 = -pkin(6) + qJ(3);
t17 = t34 * t36 + t35 * t37;
t56 = qJD(2) * t17;
t55 = qJD(2) * t34;
t54 = qJD(3) * t34;
t53 = qJD(4) * t34;
t10 = t17 * qJD(5);
t52 = t10 * qJD(5);
t49 = qJ(3) * t64 + qJD(3) * t62;
t48 = qJD(2) * t59;
t47 = t37 * t55;
t46 = 0.2e1 * qJD(2) * qJD(4);
t19 = t35 * qJD(1) - t34 * t51;
t43 = 0.2e1 * t56;
t41 = t18 * qJD(3);
t40 = t17 * qJD(3);
t15 = (pkin(3) + pkin(4)) * t35 + t45;
t39 = qJD(2) ^ 2;
t38 = qJD(5) ^ 2;
t24 = qJD(5) * t48;
t23 = t58 * t35;
t22 = t58 * t34;
t21 = t57 * t39;
t16 = qJD(4) - t19;
t14 = t47 - t48;
t11 = t18 * qJD(5);
t9 = 0.2e1 * t64;
t7 = qJD(3) - t63;
t5 = t11 * qJD(5);
t4 = qJD(5) * t47 - t24;
t3 = qJD(2) * t10;
t1 = t15 * qJD(2) - qJD(3);
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t52; 0, 0, 0, 0, 0, 0, t9, -t19 * t54 + t49, t46 * t61, t9, t32 * t46, t16 * t54 + (-t7 + t63) * t53 + t49, -t14 * t10 - t3 * t18, t10 * t56 - t14 * t11 + t3 * t17 - t18 * t4, -t52, -t5, 0, t1 * t11 + t15 * t4 + ((-t22 * t36 - t23 * t37) * qJD(5) + t41) * qJD(5) + t43 * t53, -t1 * t10 - t15 * t3 + ((-t22 * t37 + t23 * t36) * qJD(5) - t40) * qJD(5) + (qJD(2) * t18 + t14) * t53; 0, 0, 0, 0, 0, 0, -t21, (t19 * t34 - t62) * qJD(2), 0, -t21, 0, (-t62 + (-qJD(4) - t16) * t34) * qJD(2), 0, 0, 0, 0, 0, t24 + (-t14 - t47) * qJD(5), t43 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t61, 0, -t32 * t39, (qJD(3) + t7) * t55, 0, 0, 0, 0, 0, -t38 * t36 - t55 * t56, -t14 * t55 - t38 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t56, t14 ^ 2 - t56 ^ 2, 0, t24 + (t14 - t47) * qJD(5), 0, qJD(2) * t41 - t1 * t14, -qJD(2) * t40 + t1 * t56;];
tauc_reg = t2;
