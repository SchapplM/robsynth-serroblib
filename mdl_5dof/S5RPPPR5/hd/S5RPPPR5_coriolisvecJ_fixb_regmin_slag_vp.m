% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.31s
% Computational Cost: add. (295->64), mult. (661->118), div. (0->0), fcn. (418->6), ass. (0->55)
t38 = sin(pkin(8));
t40 = cos(pkin(8));
t43 = sin(qJ(5));
t44 = cos(qJ(5));
t48 = t44 * t38 + t43 * t40;
t17 = t48 * qJD(1);
t73 = -0.2e1 * t17;
t67 = t38 ^ 2 + t40 ^ 2;
t41 = cos(pkin(7));
t28 = t41 * qJD(2) - qJD(4);
t24 = qJD(1) * t28;
t55 = t67 * t24;
t45 = -pkin(1) - pkin(2);
t39 = sin(pkin(7));
t68 = t41 * qJ(2) + t39 * t45;
t22 = -qJ(4) + t68;
t72 = pkin(6) - t22;
t46 = qJD(1) ^ 2;
t71 = t39 * t46;
t70 = t41 * t46;
t69 = t43 * t38;
t26 = t45 * qJD(1) + qJD(2);
t60 = qJD(1) * qJ(2);
t14 = t39 * t26 + t41 * t60;
t66 = qJD(1) * t40;
t65 = qJD(5) * t41;
t20 = -t44 * t40 + t69;
t18 = t20 * qJD(5);
t64 = t18 * qJD(5);
t19 = t48 * qJD(5);
t63 = t19 * qJD(5);
t62 = t39 * qJD(2);
t61 = t39 * qJD(5) ^ 2;
t59 = qJD(1) * qJD(2);
t58 = qJD(1) * t69;
t57 = t44 * t66;
t56 = 0.2e1 * t59;
t54 = t39 * t59;
t13 = t41 * t26 - t39 * t60;
t53 = -t39 * qJ(2) + t41 * t45;
t52 = 0.2e1 * t54;
t51 = pkin(3) - t53;
t50 = t67 * (-qJD(1) * qJ(4) + t14);
t49 = t13 * t39 - t14 * t41;
t7 = qJD(1) * pkin(3) + qJD(4) - t13;
t47 = t50 + t62;
t25 = qJD(5) * t58;
t16 = -t57 + t58;
t15 = t40 * pkin(4) + t51;
t12 = qJD(1) * t19;
t11 = -qJD(5) * t57 + t25;
t10 = t72 * t40;
t9 = t72 * t38;
t5 = pkin(4) * t66 + t7;
t1 = [0, 0, 0, 0, t56, qJ(2) * t56, t52, t41 * t56, ((-t39 * t53 + t41 * t68) * qJD(1) - t49) * qJD(2), t40 * t52, -0.2e1 * t38 * t54, -0.2e1 * t55, t50 * t28 + t22 * t55 + (qJD(1) * t51 + t7) * t62, -t11 * t48 - t17 * t18, t11 * t20 - t12 * t48 + t18 * t16 - t17 * t19, t64, t63, 0, -t15 * t12 - t5 * t19 + ((t10 * t44 - t43 * t9) * qJD(5) - t48 * t28) * qJD(5) + (-qJD(1) * t20 - t16) * t62, t15 * t11 + t5 * t18 + ((-t10 * t43 - t44 * t9) * qJD(5) + t20 * t28) * qJD(5) + t62 * t73; 0, 0, 0, 0, -t46, -t46 * qJ(2), -t71, -t70, t49 * qJD(1), -t40 * t71, t38 * t71, t67 * t70, t39 * t55 + (-t7 * t39 - t47 * t41) * qJD(1), 0, 0, 0, 0, 0, t41 * t12 + t20 * t61 + (t39 * t16 + t48 * t65) * qJD(1), -t41 * t11 + t48 * t61 + (t39 * t17 - t20 * t65) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67 * t46, t47 * qJD(1), 0, 0, 0, 0, 0, qJD(5) * t73, t25 + (t16 - t57) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t16, -t16 ^ 2 + t17 ^ 2, t25 + (-t16 - t57) * qJD(5), 0, 0, t5 * t17 - t48 * t24, -t5 * t16 + t20 * t24;];
tauc_reg = t1;
