% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:46
% DurationCPUTime: 0.45s
% Computational Cost: add. (384->97), mult. (1071->167), div. (0->0), fcn. (769->8), ass. (0->73)
t51 = sin(qJ(5));
t45 = sin(pkin(9));
t46 = sin(pkin(8));
t75 = qJD(1) * t46;
t66 = t45 * t75;
t48 = cos(pkin(9));
t52 = cos(qJ(5));
t79 = t48 * t52;
t68 = t46 * t79;
t18 = qJD(1) * t68 - t51 * t66;
t49 = cos(pkin(8));
t41 = t49 * qJD(1) - qJD(5);
t87 = qJD(5) + t41;
t86 = pkin(6) * t46;
t33 = -cos(pkin(7)) * pkin(1) - pkin(2) - t46 * qJ(4) - t49 * pkin(3);
t13 = t33 * qJD(1) + qJD(3);
t42 = sin(pkin(7)) * pkin(1) + qJ(3);
t40 = qJD(1) * t42;
t25 = t46 * qJD(2) + t49 * t40;
t4 = t45 * t13 + t48 * t25;
t57 = t45 * t52 + t48 * t51;
t54 = t57 * qJD(5);
t19 = t46 * t54;
t11 = qJD(1) * t19;
t85 = t11 * t49;
t16 = t57 * t75;
t84 = t16 * t41;
t83 = t18 * t41;
t82 = t19 * t41;
t81 = t42 * t49;
t80 = t45 * t51;
t12 = t18 * qJD(5);
t78 = t49 * t12;
t77 = t45 * t33 + t48 * t81;
t43 = t46 ^ 2;
t44 = t49 ^ 2;
t76 = t43 + t44;
t74 = qJD(3) * t46;
t73 = qJD(3) * t49;
t72 = qJD(4) * t46;
t71 = qJD(1) * qJD(3);
t70 = t48 * t86;
t53 = qJD(1) ^ 2;
t69 = t46 * t49 * t53;
t67 = 0.2e1 * qJD(3) * t43;
t65 = t76 * t53;
t3 = t48 * t13 - t45 * t25;
t31 = -t45 * t73 - t48 * t72;
t22 = qJD(1) * t31;
t32 = -t45 * t72 + t48 * t73;
t23 = qJD(1) * t32;
t64 = t52 * t22 - t51 * t23;
t24 = t49 * qJD(2) - t46 * t40;
t1 = (-t49 * pkin(4) - t70) * qJD(1) + t3;
t2 = -pkin(6) * t66 + t4;
t61 = t52 * t1 - t51 * t2;
t60 = -t51 * t1 - t52 * t2;
t59 = t51 * t22 + t52 * t23;
t58 = t24 * t46 - t25 * t49;
t56 = t79 - t80;
t21 = qJD(4) - t24;
t55 = t56 * t41;
t36 = t43 * t42 * t71;
t30 = (pkin(4) * t45 + t42) * t46;
t29 = t48 * t33;
t27 = t56 * t46;
t26 = t57 * t46;
t20 = (-t46 * t80 + t68) * qJD(5);
t8 = pkin(4) * t66 + t21;
t7 = t20 * t41;
t6 = -t45 * t86 + t77;
t5 = -t70 + t29 + (-t42 * t45 - pkin(4)) * t49;
t9 = [0, 0, 0, 0, 0, 0.2e1 * t76 * t71, t36 + (t44 * t40 - t58) * qJD(3), -t22 * t49 + (-t31 * t49 + t45 * t67) * qJD(1), t23 * t49 + (t32 * t49 + t48 * t67) * qJD(1), t23 * t77 + t4 * t32 + t22 * (-t45 * t81 + t29) + t3 * t31 + t36 + t21 * t74, -t11 * t27 - t18 * t19, t11 * t26 - t27 * t12 + t19 * t16 - t18 * t20, t82 + t85, t7 + t78, 0, -(t52 * t31 - t51 * t32) * t41 - t64 * t49 + t30 * t12 + t8 * t20 + (-(-t5 * t51 - t52 * t6) * t41 - t60 * t49) * qJD(5) + (qJD(1) * t26 + t16) * t74, (t51 * t31 + t52 * t32) * t41 + t59 * t49 - t30 * t11 - t8 * t19 + ((t5 * t52 - t51 * t6) * t41 + t61 * t49) * qJD(5) + (qJD(1) * t27 + t18) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t22 * t45 + t23 * t48 - t49 * t71) * t46, 0, 0, 0, 0, 0, t7 - t78, -t82 + t85; 0, 0, 0, 0, 0, -t65, t58 * qJD(1), -t45 * t65, -t48 * t65, t22 * t48 + t23 * t45 + (-t21 * t46 + (t3 * t45 - t4 * t48) * t49) * qJD(1), 0, 0, 0, 0, 0, t41 * t54 + (-t57 * t41 * t49 - t46 * t16) * qJD(1), qJD(5) * t55 + (-t46 * t18 - t49 * t55) * qJD(1); 0, 0, 0, 0, 0, 0, 0, -t48 * t69, t45 * t69, (t3 * t48 + t4 * t45 + qJD(3)) * t75, 0, 0, 0, 0, 0, t12 - t83, -t11 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t16, -t16 ^ 2 + t18 ^ 2, -t11 - t84, -t12 - t83, 0, -t8 * t18 + t87 * t60 + t64, t8 * t16 - t87 * t61 - t59;];
tauc_reg = t9;
