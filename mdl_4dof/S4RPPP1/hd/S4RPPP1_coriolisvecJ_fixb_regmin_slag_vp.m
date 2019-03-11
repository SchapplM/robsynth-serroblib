% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:26
% DurationCPUTime: 0.26s
% Computational Cost: add. (158->80), mult. (607->138), div. (0->0), fcn. (394->4), ass. (0->65)
t32 = cos(pkin(4));
t69 = pkin(1) * t32;
t31 = cos(pkin(6));
t30 = sin(pkin(4));
t56 = qJD(2) * t30;
t16 = t32 * qJD(3) + t31 * t56;
t13 = qJD(1) * t16;
t68 = t13 * t32;
t67 = t16 * t32;
t26 = t30 ^ 2;
t33 = qJD(1) ^ 2;
t66 = t26 * t33;
t29 = sin(pkin(6));
t65 = t29 * t30;
t64 = t30 * t31;
t52 = t30 * qJD(1);
t45 = qJ(2) * t52;
t57 = qJD(1) * t32;
t47 = pkin(1) * t57;
t9 = t29 * t47 + t31 * t45;
t61 = qJ(2) * t30;
t63 = t29 * t69 + t31 * t61;
t25 = t29 ^ 2;
t27 = t31 ^ 2;
t62 = t25 + t27;
t60 = t32 * qJ(3);
t46 = -pkin(1) * t31 - pkin(2);
t40 = -qJ(4) + t46;
t19 = t29 * t45;
t50 = qJD(3) + t19;
t1 = (pkin(3) * t65 + t40 * t32) * qJD(1) + t50;
t59 = qJD(3) + t1;
t2 = qJD(4) + (pkin(3) * t64 + t60) * qJD(1) + t9;
t58 = -qJD(4) - t2;
t55 = qJD(3) * t26;
t54 = t29 * qJD(2);
t53 = t29 * qJD(3);
t51 = t32 * qJD(4);
t49 = t30 * t32 * t33;
t48 = 0.2e1 * qJD(2) * t26;
t44 = qJD(2) * t52;
t43 = -qJ(3) * t29 - pkin(1);
t42 = t29 * t31 * t66;
t41 = t31 * t49;
t39 = t46 * t32;
t38 = (t31 * t47 - t19) * t29 - t9 * t31;
t37 = -0.2e1 * t32 * t44;
t18 = t29 * t44;
t12 = -qJD(1) * t51 + t18;
t36 = t12 * t29 + t13 * t31;
t11 = (-qJD(4) * t31 - t53) * t30;
t35 = (-pkin(2) * t31 + t43) * t30;
t34 = (-pkin(2) - qJ(4)) * t31 + t43;
t28 = t32 ^ 2;
t22 = t29 * t61;
t17 = t29 * t49;
t15 = t30 * t54 - t51;
t14 = (-t25 * t26 - t28) * t33;
t10 = t62 * t66;
t7 = qJD(1) * t11;
t6 = -qJ(3) * t57 - t9;
t5 = qJD(1) * t35 + qJD(2);
t4 = qJD(1) * t39 + t50;
t3 = t34 * t52 + qJD(2);
t8 = [0, 0, 0, t29 * t37, t31 * t37, t62 * qJD(1) * t48 ((t31 * t63 - t29 * (t31 * t69 - t22)) * qJD(1) - t38) * t56, t13 * t64 + (t16 * t64 + t25 * t48) * qJD(1), 0.2e1 * (-t31 * t55 + t32 * t56) * t29 * qJD(1), t68 + (0.2e1 * t25 * t55 + t67) * qJD(1), -t13 * (-t60 - t63) - t6 * t16 + ((t4 * qJD(2) - t5 * qJD(3)) * t29 + ((t22 + t39) * t54 - t35 * t53) * qJD(1)) * t30 ((t15 * t29 + t16 * t31) * qJD(1) + t36) * t30, -t7 * t65 + t68 + (-t11 * t65 + t67) * qJD(1), -t7 * t64 - t12 * t32 + (-t11 * t64 - t15 * t32) * qJD(1), t3 * t11 + t12 * t22 + t1 * t15 + t13 * t63 + t2 * t16 + (t13 * qJ(3) + t12 * t40) * t32 + (t36 * pkin(3) + t7 * t34) * t30; 0, 0, 0, t17, t41, -t10, t38 * t52, -t10, -t17, -t41 (t31 * t6 + (-qJD(3) - t4) * t29) * t52, -t10, -t41, t17 (-t59 * t29 + t58 * t31) * t52; 0, 0, 0, 0, 0, 0, 0, -t41, t42, t14, t18 + (t32 * t6 + t5 * t65) * qJD(1), -t41, t14, -t42, t18 + (t3 * t65 + t58 * t32) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t42 (-t26 * t27 - t28) * t33 (t59 * t32 + (qJD(2) + t3) * t64) * qJD(1);];
tauc_reg  = t8;
