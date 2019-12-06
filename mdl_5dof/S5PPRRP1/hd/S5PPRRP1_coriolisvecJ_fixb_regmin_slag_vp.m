% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauc_reg [5x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:15
% DurationCPUTime: 0.38s
% Computational Cost: add. (321->65), mult. (851->109), div. (0->0), fcn. (586->6), ass. (0->54)
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t68 = -t29 * t32 + t30 * t34;
t14 = t68 * qJD(1);
t19 = t29 * t34 + t30 * t32;
t17 = t19 * qJD(3);
t59 = -qJ(5) - pkin(6);
t31 = sin(qJ(4));
t33 = cos(qJ(4));
t15 = t19 * qJD(1);
t43 = -t59 * qJD(3) + t15;
t4 = t33 * qJD(2) - t43 * t31;
t11 = qJD(1) * t17;
t67 = qJD(3) * t15 - t11;
t5 = qJD(2) * t31 + t43 * t33;
t55 = qJD(4) * pkin(4);
t3 = t4 + t55;
t66 = t3 - t4;
t35 = qJD(4) ^ 2;
t61 = t35 * t31;
t60 = t35 * t33;
t27 = t31 ^ 2;
t28 = t33 ^ 2;
t58 = t27 - t28;
t57 = t27 + t28;
t56 = qJD(3) * pkin(3);
t52 = qJD(4) * t31;
t51 = qJD(4) * t33;
t16 = t68 * qJD(3);
t50 = t16 * qJD(3);
t48 = qJD(3) * qJD(4);
t46 = t31 * t48;
t7 = pkin(4) * t46 + t11;
t47 = -pkin(4) * t33 - pkin(3);
t10 = qJD(1) * t16;
t8 = -t14 - t56;
t45 = -qJD(3) * t8 - t10;
t44 = qJD(4) * t59;
t42 = qJD(3) * qJD(5) + t10;
t41 = t3 * t31 - t33 * t5;
t39 = pkin(6) * t35 - t67;
t38 = qJD(4) * (t14 + t8 - t56);
t1 = t4 * qJD(4) + t42 * t33;
t2 = -qJD(4) * t5 - t42 * t31;
t37 = t1 * t33 - t2 * t31 + (-t3 * t33 - t31 * t5) * qJD(4);
t36 = qJD(3) ^ 2;
t21 = t59 * t33;
t20 = t59 * t31;
t13 = -qJD(5) * t31 + t33 * t44;
t12 = qJD(5) * t33 + t31 * t44;
t6 = t47 * qJD(3) + qJD(5) - t14;
t9 = [0, 0, 0, -t17 * qJD(3), -t50, 0, 0, 0, 0, 0, -t16 * t52 - t19 * t60 + (-t17 * t33 - t52 * t68) * qJD(3), -t16 * t51 + t19 * t61 + (t17 * t31 - t51 * t68) * qJD(3), t57 * t50, -t41 * t16 + t17 * t6 + t37 * t19 - t68 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t60, 0, -t41 * qJD(4) + t1 * t31 + t2 * t33; 0, 0, 0, t67, 0, 0.2e1 * t33 * t46, -0.2e1 * t58 * t48, t60, -t61, 0, t31 * t38 - t39 * t33, t39 * t31 + t33 * t38, (t12 * t33 - t13 * t31 - t57 * t14 + (-t20 * t33 + t21 * t31) * qJD(4)) * qJD(3) + t37, -t1 * t21 + t5 * t12 + t2 * t20 + t3 * t13 + t7 * t47 + (pkin(4) * t52 - t15) * t6 + t41 * t14; 0, 0, 0, 0, 0, -t31 * t36 * t33, t58 * t36, 0, 0, 0, t45 * t31, t45 * t33, (-t55 + t66) * t33 * qJD(3), t66 * t5 + (-qJD(3) * t31 * t6 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t36, t41 * qJD(3) + t7;];
tauc_reg = t9;
