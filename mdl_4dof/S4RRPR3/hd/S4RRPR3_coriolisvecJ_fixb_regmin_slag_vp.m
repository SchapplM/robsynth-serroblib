% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (175->47), mult. (439->78), div. (0->0), fcn. (236->6), ass. (0->48)
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t28 = qJD(1) + qJD(2);
t36 = cos(qJ(2));
t52 = pkin(1) * qJD(1);
t47 = t36 * t52;
t19 = t28 * pkin(2) + t47;
t31 = sin(pkin(7));
t32 = cos(pkin(7));
t34 = sin(qJ(2));
t48 = t34 * t52;
t5 = t32 * t19 - t31 * t48;
t3 = -t28 * pkin(3) - t5;
t50 = qJD(4) * t3;
t51 = pkin(1) * qJD(2);
t57 = t32 * t34;
t14 = (t31 * t36 + t57) * t51;
t8 = qJD(1) * t14;
t59 = t8 * t33 + t35 * t50;
t58 = t31 * t34;
t56 = t33 * t35;
t37 = qJD(4) ^ 2;
t55 = t37 * t33;
t25 = t36 * pkin(1) + pkin(2);
t54 = pkin(1) * t57 + t31 * t25;
t53 = t33 ^ 2 - t35 ^ 2;
t49 = 0.2e1 * qJD(4) * t28;
t38 = pkin(1) * (t32 * t36 - t58);
t16 = qJD(2) * t38;
t9 = qJD(1) * t16;
t46 = -t3 * t28 - t9;
t45 = (-qJD(2) + t28) * t52;
t44 = (-qJD(1) - t28) * t51;
t43 = (pkin(6) + t54) * t37 + t14 * t28;
t21 = t32 * t48;
t13 = t31 * t47 + t21;
t42 = -t13 * t28 + (t31 * pkin(2) + pkin(6)) * t37;
t39 = -pkin(1) * t58 + t32 * t25;
t41 = qJD(4) * ((-pkin(3) - t39) * t28 - t16);
t15 = qJD(1) * t38;
t40 = qJD(4) * ((-t32 * pkin(2) - pkin(3)) * t28 + t15);
t27 = t28 ^ 2;
t26 = t37 * t35;
t18 = t49 * t56;
t10 = t53 * t49;
t6 = t31 * t19 + t21;
t1 = t33 * t50;
t2 = [0, 0, 0, 0, t34 * t44, t36 * t44, -t5 * t14 + t6 * t16 - t8 * t39 + t9 * t54, t18, -t10, t26, -t55, 0, t1 + t33 * t41 + (-t43 - t8) * t35, t43 * t33 + t35 * t41 + t59; 0, 0, 0, 0, t34 * t45, t36 * t45, t5 * t13 - t6 * t15 + (t31 * t9 - t32 * t8) * pkin(2), t18, -t10, t26, -t55, 0, t1 + t33 * t40 + (-t42 - t8) * t35, t42 * t33 + t35 * t40 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t26; 0, 0, 0, 0, 0, 0, 0, -t27 * t56, t53 * t27, 0, 0, 0, t46 * t33, t46 * t35;];
tauc_reg = t2;
