% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:41
% EndTime: 2021-01-15 17:04:43
% DurationCPUTime: 0.28s
% Computational Cost: add. (337->80), mult. (656->114), div. (0->0), fcn. (294->4), ass. (0->60)
t20 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t12 = t20 * qJD(1) + qJD(3);
t31 = sin(qJ(4));
t46 = qJD(2) * qJD(4);
t23 = t31 * t46;
t32 = cos(qJ(4));
t50 = t32 * qJD(4);
t45 = qJ(5) * t50;
t51 = t31 * qJD(5);
t1 = t12 * t50 - t23 + (-t45 - t51) * qJD(1);
t47 = qJD(1) * qJD(4);
t44 = t31 * t47;
t19 = qJ(5) * t44;
t62 = t31 * t12;
t35 = -t32 * qJD(2) - t62;
t49 = t32 * qJD(5);
t2 = -qJD(1) * t49 + t35 * qJD(4) + t19;
t8 = t32 * t12;
t41 = -t31 * qJD(2) + t8;
t48 = qJ(5) * qJD(1);
t4 = -t32 * t48 + t41;
t58 = qJD(4) * pkin(4);
t3 = t4 + t58;
t5 = -t31 * t48 - t35;
t65 = -t1 * t31 - t2 * t32 + (t3 * t31 - t32 * t5) * qJD(4);
t64 = 0.2e1 * qJD(3);
t63 = t3 - t4;
t33 = qJD(4) ^ 2;
t25 = t33 * t31;
t61 = t33 * t32;
t34 = qJD(1) ^ 2;
t60 = t34 * t31;
t26 = qJD(3) * qJD(1);
t43 = t32 * t47;
t13 = pkin(4) * t43 + t26;
t27 = t31 ^ 2;
t28 = t32 ^ 2;
t59 = t27 - t28;
t22 = sin(pkin(7)) * pkin(1) + qJ(3);
t14 = t31 * pkin(4) + t22;
t54 = qJD(1) * t14;
t9 = qJD(5) + t54;
t57 = t9 * qJD(1);
t56 = qJ(5) - t20;
t55 = qJD(5) + t9;
t17 = qJD(1) * t22;
t53 = t17 * qJD(1);
t52 = t31 * qJD(4);
t11 = t56 * t32;
t42 = t9 + t54;
t18 = pkin(4) * t50 + qJD(3);
t40 = qJD(1) * t18 + t13;
t39 = -0.2e1 * t44;
t37 = t3 * t32 + t31 * t5;
t16 = (-t33 - t34) * t32;
t15 = -t25 - t60;
t10 = t56 * t31;
t7 = -qJD(4) * t11 - t51;
t6 = t56 * t52 - t49;
t21 = [0, 0, 0, 0, 0, 0.2e1 * t26, t17 * t64, t32 * t39, 0.2e1 * t59 * t47, -t25, -t61, 0, t17 * t50 - t20 * t25 + (t22 * t50 + t31 * t64) * qJD(1), -t17 * t52 - t20 * t61 + (-t22 * t52 + t32 * t64) * qJD(1), t40 * t31 + (t42 * t32 + t6) * qJD(4), t40 * t32 + (-t42 * t31 - t7) * qJD(4), (-t31 * t7 - t32 * t6 + (t10 * t32 - t11 * t31) * qJD(4)) * qJD(1) + t65, -t1 * t10 - t2 * t11 + t13 * t14 + t9 * t18 + t3 * t6 + t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t25, -t61, t25, 0, -t37 * qJD(4) + t1 * t32 - t2 * t31; 0, 0, 0, 0, 0, -t34, -t53, 0, 0, 0, 0, 0, t15, t16, t15, t16, 0, -t57 - t65; 0, 0, 0, 0, 0, 0, 0, t32 * t60, -t59 * t34, 0, 0, 0, -t32 * t53, t31 * t53 + t23 + (t41 - t8) * qJD(4), t19 + (t5 - t62) * qJD(4) + (-pkin(4) * t60 - t55 * qJD(1) - t46) * t32, -t28 * t34 * pkin(4) + t23 + (t4 - t8) * qJD(4) + (t55 * t31 + t45) * qJD(1), (t58 - t63) * t31 * qJD(1), t63 * t5 + (-t32 * t57 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t43, t39, (-t27 - t28) * t34, t37 * qJD(1) + t13;];
tauc_reg = t21;
