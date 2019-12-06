% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:09
% EndTime: 2019-12-05 15:22:12
% DurationCPUTime: 0.44s
% Computational Cost: add. (329->98), mult. (1002->173), div. (0->0), fcn. (734->6), ass. (0->73)
t47 = sin(qJ(5));
t43 = sin(pkin(9));
t44 = sin(pkin(8));
t73 = qJD(2) * t44;
t63 = t43 * t73;
t45 = cos(pkin(9));
t48 = cos(qJ(5));
t77 = t45 * t48;
t65 = t44 * t77;
t14 = qJD(2) * t65 - t47 * t63;
t46 = cos(pkin(8));
t39 = t46 * qJD(2) - qJD(5);
t85 = qJD(5) + t39;
t84 = pkin(6) * t44;
t9 = t14 * qJD(5);
t83 = t46 * t9;
t53 = t43 * t48 + t45 * t47;
t50 = t53 * qJD(5);
t15 = t44 * t50;
t8 = qJD(2) * t15;
t82 = t8 * t46;
t12 = t53 * t73;
t81 = t12 * t39;
t80 = t14 * t39;
t79 = t15 * t39;
t78 = t43 * t47;
t34 = -t46 * pkin(3) - t44 * qJ(4) - pkin(2);
t23 = t34 * qJD(2) + qJD(3);
t69 = qJ(3) * qJD(2);
t32 = t44 * qJD(1) + t46 * t69;
t4 = t43 * t23 + t45 * t32;
t74 = qJ(3) * t46;
t76 = t43 * t34 + t45 * t74;
t41 = t44 ^ 2;
t42 = t46 ^ 2;
t75 = t41 + t42;
t72 = qJD(3) * t44;
t71 = qJD(3) * t46;
t70 = qJD(4) * t44;
t68 = qJD(2) * qJD(3);
t67 = t45 * t84;
t49 = qJD(2) ^ 2;
t66 = t44 * t46 * t49;
t64 = 0.2e1 * qJD(3) * t41;
t62 = t75 * t49;
t24 = -t43 * t71 - t45 * t70;
t19 = qJD(2) * t24;
t25 = -t43 * t70 + t45 * t71;
t20 = qJD(2) * t25;
t61 = t48 * t19 - t47 * t20;
t3 = t45 * t23 - t43 * t32;
t31 = t46 * qJD(1) - t44 * t69;
t1 = (-t46 * pkin(4) - t67) * qJD(2) + t3;
t2 = -pkin(6) * t63 + t4;
t58 = t48 * t1 - t47 * t2;
t57 = -t47 * t1 - t48 * t2;
t56 = -t19 * t45 - t20 * t43;
t55 = t47 * t19 + t48 * t20;
t54 = t31 * t44 - t32 * t46;
t52 = t77 - t78;
t29 = qJD(4) - t31;
t51 = t52 * t39;
t38 = t41 * qJ(3) * t68;
t30 = (pkin(4) * t43 + qJ(3)) * t44;
t28 = t45 * t34;
t22 = t52 * t44;
t21 = t53 * t44;
t16 = (-t44 * t78 + t65) * qJD(5);
t11 = pkin(4) * t63 + t29;
t7 = t16 * t39;
t6 = -t43 * t84 + t76;
t5 = -t67 + t28 + (-qJ(3) * t43 - pkin(4)) * t46;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t19 * t43 + t20 * t45 - t46 * t68) * t44, 0, 0, 0, 0, 0, t7 - t83, -t79 + t82; 0, 0, 0, 0, 0, 0, 0.2e1 * t75 * t68, t38 + (t42 * t69 - t54) * qJD(3), -t19 * t46 + (-t24 * t46 + t43 * t64) * qJD(2), t20 * t46 + (t25 * t46 + t45 * t64) * qJD(2), ((-t24 * t45 - t25 * t43) * qJD(2) + t56) * t44, t20 * t76 + t4 * t25 + t19 * (-t43 * t74 + t28) + t3 * t24 + t38 + t29 * t72, -t14 * t15 - t8 * t22, t15 * t12 - t14 * t16 + t8 * t21 - t22 * t9, t79 + t82, t7 + t83, 0, -(t48 * t24 - t47 * t25) * t39 - t61 * t46 + t30 * t9 + t11 * t16 + (-(-t47 * t5 - t48 * t6) * t39 - t57 * t46) * qJD(5) + (qJD(2) * t21 + t12) * t72, (t47 * t24 + t48 * t25) * t39 + t55 * t46 - t30 * t8 - t11 * t15 + ((-t47 * t6 + t48 * t5) * t39 + t58 * t46) * qJD(5) + (qJD(2) * t22 + t14) * t72; 0, 0, 0, 0, 0, 0, -t62, t54 * qJD(2), -t43 * t62, -t45 * t62, 0, (-t29 * t44 + (t3 * t43 - t4 * t45) * t46) * qJD(2) - t56, 0, 0, 0, 0, 0, t39 * t50 + (-t53 * t39 * t46 - t44 * t12) * qJD(2), qJD(5) * t51 + (-t44 * t14 - t46 * t51) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t66, t43 * t66, (-t43 ^ 2 - t45 ^ 2) * t49 * t41, (t3 * t45 + t4 * t43 + qJD(3)) * t73, 0, 0, 0, 0, 0, t9 - t80, -t8 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t12, -t12 ^ 2 + t14 ^ 2, -t8 - t81, -t9 - t80, 0, -t11 * t14 + t85 * t57 + t61, t11 * t12 - t85 * t58 - t55;];
tauc_reg = t10;
