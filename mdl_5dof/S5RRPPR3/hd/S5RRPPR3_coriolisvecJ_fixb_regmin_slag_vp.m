% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:37
% DurationCPUTime: 0.25s
% Computational Cost: add. (295->72), mult. (657->101), div. (0->0), fcn. (344->6), ass. (0->58)
t35 = sin(pkin(8));
t36 = cos(pkin(8));
t40 = cos(qJ(2));
t61 = pkin(1) * qJD(2);
t54 = qJD(1) * t61;
t49 = t40 * t54;
t38 = sin(qJ(2));
t50 = t38 * t54;
t10 = t35 * t49 + t36 * t50;
t32 = qJD(1) + qJD(2);
t62 = pkin(1) * qJD(1);
t56 = t40 * t62;
t20 = t32 * pkin(2) + t56;
t57 = t38 * t62;
t9 = t35 * t20 + t36 * t57;
t4 = t32 * qJ(4) + t9;
t70 = -t4 * t32 + t10;
t31 = t32 ^ 2;
t37 = sin(qJ(5));
t39 = cos(qJ(5));
t59 = t39 * qJD(5);
t11 = -t35 * t50 + t36 * t49;
t7 = t32 * qJD(4) + t11;
t69 = t7 * t37 + t4 * t59;
t67 = t36 * t38;
t41 = qJD(5) ^ 2;
t66 = t41 * t37;
t65 = t41 * t39;
t64 = -t31 - t41;
t63 = t37 ^ 2 - t39 ^ 2;
t60 = t37 * qJD(5);
t25 = t35 * t57;
t18 = t36 * t56 - t25;
t58 = qJD(4) - t18;
t27 = t35 * t38 * pkin(1);
t55 = -t36 * pkin(2) - pkin(3);
t30 = t40 * pkin(1) + pkin(2);
t43 = pkin(1) * t67 + t35 * t30;
t15 = qJ(4) + t43;
t42 = pkin(1) * (t35 * t40 + t67);
t17 = qJD(2) * t42;
t52 = t15 * t32 + t17;
t8 = t36 * t20 - t25;
t51 = t36 * t30 - t27;
t48 = -pkin(3) - t51;
t47 = (-qJD(2) + t32) * t62;
t46 = (-qJD(1) - t32) * t61;
t44 = t36 * t40 * t61 - qJD(2) * t27;
t13 = qJD(4) + t44;
t45 = t13 * t32 - (-pkin(7) + t48) * t41;
t29 = t35 * pkin(2) + qJ(4);
t28 = -pkin(7) + t55;
t19 = -0.2e1 * t32 * t37 * t59;
t16 = qJD(1) * t42;
t12 = 0.2e1 * t63 * t32 * qJD(5);
t6 = t7 * t39;
t3 = -t32 * pkin(3) + qJD(4) - t8;
t1 = [0, 0, 0, 0, t38 * t46, t40 * t46, -t10 * t51 + t11 * t43 - t8 * t17 + t9 * t44, t17 * t32 + t10, (qJD(4) + t13) * t32 + t11, t10 * t48 + t4 * t13 + t7 * t15 + t3 * t17, t19, t12, -t66, -t65, 0, t45 * t37 + t52 * t59 + t69, t6 + t45 * t39 + (-t4 - t52) * t60; 0, 0, 0, 0, t38 * t47, t40 * t47, t8 * t16 - t9 * t18 + (-t10 * t36 + t11 * t35) * pkin(2), -t16 * t32 + t10, (0.2e1 * qJD(4) - t18) * t32 + t11, t10 * t55 - t3 * t16 + t7 * t29 + t58 * t4, t19, t12, -t66, -t65, 0, -t16 * t59 - t28 * t66 + (t29 * t59 + t58 * t37) * t32 + t69, t6 + (-t28 * t41 + t58 * t32) * t39 + (-t29 * t32 + t16 - t4) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t66; 0, 0, 0, 0, 0, 0, 0, 0, -t31, t70, 0, 0, 0, 0, 0, t64 * t37, t64 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t31 * t37, -t63 * t31, 0, 0, 0, t70 * t39, -t70 * t37;];
tauc_reg = t1;
