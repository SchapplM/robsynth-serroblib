% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:29
% EndTime: 2021-01-15 14:56:31
% DurationCPUTime: 0.29s
% Computational Cost: add. (278->88), mult. (698->128), div. (0->0), fcn. (374->4), ass. (0->69)
t30 = sin(qJ(3));
t33 = qJD(4) ^ 2;
t34 = qJD(3) ^ 2;
t72 = (t33 + t34) * t30;
t57 = t30 * qJD(2);
t20 = qJD(3) * t57;
t29 = sin(qJ(4));
t52 = qJD(3) * qJD(4);
t46 = t29 * t52;
t11 = pkin(4) * t46 + t20;
t17 = qJD(3) * pkin(6) + t57;
t44 = qJ(5) * qJD(3) + t17;
t31 = cos(qJ(4));
t56 = t31 * qJD(1);
t6 = -t29 * t44 - t56;
t60 = qJD(4) * pkin(4);
t3 = t6 + t60;
t24 = t29 * qJD(1);
t7 = t31 * t44 - t24;
t40 = t29 * t3 - t31 * t7;
t71 = qJD(3) * t40 + t11;
t70 = t3 - t6;
t27 = t29 ^ 2;
t69 = pkin(4) * t27;
t68 = t31 * pkin(4);
t25 = t33 * t29;
t26 = t33 * t31;
t32 = cos(qJ(3));
t67 = t34 * t32;
t66 = qJ(5) + pkin(6);
t53 = qJD(1) * qJD(4);
t58 = t17 * qJD(4);
t65 = t29 * t58 + t31 * t53;
t28 = t31 ^ 2;
t64 = t27 - t28;
t63 = t27 + t28;
t61 = qJD(3) * pkin(3);
t23 = -pkin(3) - t68;
t55 = t32 * qJD(2);
t10 = qJD(3) * t23 + qJD(5) - t55;
t59 = qJD(3) * t10;
t54 = qJ(5) * qJD(4);
t51 = t29 * t34 * t31;
t50 = 0.2e1 * t52;
t49 = t29 * t54;
t48 = t31 * t54;
t47 = qJD(4) * t55;
t45 = qJD(4) * t66;
t43 = -0.2e1 * t32 * t52;
t42 = t31 * t50;
t41 = qJD(5) + t55;
t39 = -t11 + t20;
t38 = -t10 - t41;
t18 = -t55 - t61;
t37 = qJD(4) * (t18 - t61);
t36 = qJD(3) * (-t18 - t55);
t1 = (t31 * t41 - t49) * qJD(3) - t65;
t21 = t29 * t53;
t2 = -t31 * t58 + t21 + (-t29 * t41 - t48) * qJD(3);
t35 = t1 * t31 - t2 * t29 + (-t29 * t7 - t3 * t31) * qJD(4);
t16 = t31 * t47;
t15 = t29 * t47;
t14 = t66 * t31;
t13 = t66 * t29;
t9 = -t29 * qJD(5) - t31 * t45;
t8 = t31 * qJD(5) - t29 * t45;
t5 = t29 * t43 - t31 * t72;
t4 = t29 * t72 + t31 * t43;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, t25, t26, 0, qJD(4) * t40 - t1 * t29 - t2 * t31; 0, 0, 0, -t34 * t30, -t67, 0, 0, 0, 0, 0, t5, t4, t5, t4, t63 * t67, -t71 * t32 + (t35 + t59) * t30; 0, 0, 0, 0, 0, t29 * t42, -t64 * t50, t26, -t25, 0, -pkin(6) * t26 + t29 * t37 + t15, pkin(6) * t25 + t31 * t37 + t16, t15 + t39 * t31 + (t9 + (t10 + (t23 - t68) * qJD(3)) * t29) * qJD(4), t16 - t39 * t29 + (t10 * t31 - t8 + (t23 * t31 + t69) * qJD(3)) * qJD(4), (-t29 * t9 + t31 * t8 + (t13 * t31 - t14 * t29) * qJD(4) - t63 * t55) * qJD(3) + t35, t10 * t29 * t60 + t1 * t14 + t11 * t23 - t2 * t13 + t3 * t9 + t7 * t8 + (-t10 * t30 + t32 * t40) * qJD(2); 0, 0, 0, 0, 0, -t51, t64 * t34, 0, 0, 0, -t24 * qJD(4) + t29 * t36 + t21, (-t17 * t29 - t56) * qJD(4) + t31 * t36 + t65, pkin(4) * t51 + t21 + (-t17 * t31 + t7) * qJD(4) + (t29 * t38 - t48) * qJD(3), -t34 * t69 + t6 * qJD(4) + (t31 * t38 + t49) * qJD(3) + t65, (-t60 + t70) * t31 * qJD(3), t70 * t7 + (-t29 * t59 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t46, t42, -t63 * t34, t71;];
tauc_reg = t12;
