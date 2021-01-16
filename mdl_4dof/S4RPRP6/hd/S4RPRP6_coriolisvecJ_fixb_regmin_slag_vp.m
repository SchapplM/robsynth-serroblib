% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:40
% EndTime: 2021-01-15 10:27:42
% DurationCPUTime: 0.23s
% Computational Cost: add. (203->63), mult. (430->103), div. (0->0), fcn. (169->2), ass. (0->56)
t26 = -pkin(1) - pkin(5);
t15 = t26 * qJD(1) + qJD(2);
t24 = sin(qJ(3));
t38 = qJD(1) * qJD(3);
t35 = t24 * t38;
t17 = qJ(4) * t35;
t25 = cos(qJ(3));
t42 = t25 * qJD(4);
t46 = qJD(3) * t24;
t1 = -qJD(1) * t42 - t15 * t46 + t17;
t45 = qJD(3) * t25;
t36 = qJ(4) * t45;
t43 = t24 * qJD(4);
t2 = t15 * t45 + (-t36 - t43) * qJD(1);
t50 = qJD(3) * pkin(3);
t39 = qJ(4) * qJD(1);
t8 = t25 * t15;
t7 = -t25 * t39 + t8;
t3 = t7 + t50;
t56 = t24 * t15;
t6 = -t24 * t39 + t56;
t59 = -t1 * t25 - t2 * t24 + (t24 * t3 - t25 * t6) * qJD(3);
t21 = qJD(1) * qJD(2);
t58 = 0.2e1 * t21;
t57 = t3 - t7;
t27 = qJD(3) ^ 2;
t55 = t27 * t24;
t54 = t27 * t25;
t28 = qJD(1) ^ 2;
t53 = t28 * t24;
t34 = t25 * t38;
t9 = pkin(3) * t34 + t21;
t22 = t24 ^ 2;
t23 = t25 ^ 2;
t52 = t22 - t23;
t51 = -t27 - t28;
t49 = t28 * qJ(2);
t48 = qJ(4) - t26;
t19 = t24 * pkin(3) + qJ(2);
t47 = qJD(1) * t19;
t10 = qJD(4) + t47;
t44 = t10 * qJD(1);
t41 = qJD(4) + t10;
t40 = qJ(2) * qJD(3);
t37 = 0.2e1 * qJD(1);
t12 = t48 * t25;
t16 = pkin(3) * t45 + qJD(2);
t33 = qJD(1) * t16 + t9;
t32 = t10 + t47;
t31 = -0.2e1 * t35;
t14 = t51 * t25;
t13 = t51 * t24;
t11 = t48 * t24;
t5 = -qJD(3) * t12 - t43;
t4 = t48 * t46 - t42;
t18 = [0, 0, 0, 0, t58, qJ(2) * t58, t25 * t31, 0.2e1 * t52 * t38, -t55, -t54, 0, -t26 * t55 + (qJD(2) * t24 + t25 * t40) * t37, -t26 * t54 + (qJD(2) * t25 - t24 * t40) * t37, t33 * t24 + (t32 * t25 + t4) * qJD(3), t33 * t25 + (-t32 * t24 - t5) * qJD(3), (-t24 * t5 - t25 * t4 + (t11 * t25 - t12 * t24) * qJD(3)) * qJD(1) + t59, -t1 * t12 + t10 * t16 - t2 * t11 + t9 * t19 + t3 * t4 + t6 * t5; 0, 0, 0, 0, -t28, -t49, 0, 0, 0, 0, 0, t13, t14, t13, t14, 0, -t44 - t59; 0, 0, 0, 0, 0, 0, t25 * t53, -t52 * t28, 0, 0, 0, -t25 * t49, t24 * t49, t17 + (t6 - t56) * qJD(3) + (-pkin(3) * t53 - t41 * qJD(1)) * t25, -t23 * t28 * pkin(3) + (t7 - t8) * qJD(3) + (t41 * t24 + t36) * qJD(1), (t50 - t57) * t24 * qJD(1), t57 * t6 + (-t25 * t44 + t1) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, t31, (-t22 - t23) * t28, (t24 * t6 + t25 * t3) * qJD(1) + t9;];
tauc_reg = t18;
