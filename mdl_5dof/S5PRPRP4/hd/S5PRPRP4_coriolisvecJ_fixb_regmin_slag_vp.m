% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:06
% EndTime: 2019-12-05 15:36:08
% DurationCPUTime: 0.35s
% Computational Cost: add. (484->100), mult. (1158->136), div. (0->0), fcn. (774->6), ass. (0->71)
t42 = sin(pkin(8));
t43 = cos(pkin(8));
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t28 = t42 * t45 - t43 * t47;
t25 = t28 * qJD(2);
t21 = qJD(1) * t25;
t83 = qJD(3) * qJD(4) - t21;
t46 = cos(qJ(4));
t65 = t47 * qJD(1);
t32 = qJD(2) * pkin(2) + t65;
t70 = qJD(1) * t45;
t34 = t43 * t70;
t19 = t42 * t32 + t34;
t13 = qJD(2) * pkin(6) + t19;
t44 = sin(qJ(4));
t78 = t44 * t13;
t9 = t46 * qJD(3) - t78;
t82 = qJD(5) - t9;
t6 = -qJD(4) * pkin(4) + t82;
t10 = t44 * qJD(3) + t46 * t13;
t7 = qJD(4) * qJ(5) + t10;
t40 = t44 ^ 2;
t41 = t46 ^ 2;
t81 = qJD(2) * (t40 + t41);
t80 = t43 * pkin(2);
t35 = t42 * pkin(2) + pkin(6);
t48 = qJD(4) ^ 2;
t79 = t35 * t48;
t77 = t44 * t46;
t76 = t48 * t44;
t39 = t48 * t46;
t75 = t83 * t46;
t24 = t42 * t65 + t34;
t33 = t42 * t70;
t26 = t43 * t65 - t33;
t67 = qJD(4) * t44;
t68 = qJD(2) * t46;
t74 = t24 * t68 + t26 * t67;
t53 = pkin(4) * t44 - qJ(5) * t46;
t22 = t53 * qJD(4) - t44 * qJD(5);
t73 = t22 - t24;
t72 = t40 - t41;
t69 = qJD(2) * t44;
t66 = qJD(4) * t46;
t49 = qJD(2) ^ 2;
t62 = t49 * t77;
t4 = t13 * t66 + t83 * t44;
t61 = 0.2e1 * qJD(2) * qJD(4);
t60 = t9 + t78;
t29 = t42 * t47 + t43 * t45;
t5 = (t29 * qJD(1) + t22) * qJD(2);
t59 = -t5 - t79;
t23 = t29 * qJD(2);
t20 = qJD(1) * t23;
t58 = -t20 - t79;
t18 = t43 * t32 - t33;
t51 = -t46 * pkin(4) - t44 * qJ(5) - pkin(3);
t27 = t51 - t80;
t8 = t51 * qJD(2) - t18;
t57 = qJD(2) * t27 + t8;
t12 = -qJD(2) * pkin(3) - t18;
t56 = qJD(2) * (-pkin(3) - t80) + t12;
t54 = t44 * t6 + t46 * t7;
t52 = t10 * qJD(4) - t4;
t3 = (qJD(5) - t78) * qJD(4) + t75;
t50 = t3 * t46 + t4 * t44 + (-t44 * t7 + t46 * t6) * qJD(4);
t30 = t53 * qJD(2);
t2 = t25 * t67 - t29 * t39 + (-t23 * t46 + t28 * t67) * qJD(2);
t1 = t25 * t66 + t29 * t76 + (t23 * t44 + t28 * t66) * qJD(2);
t11 = [0, 0, -t49 * t45, -t49 * t47, -t18 * t23 - t19 * t25 + t20 * t28 - t21 * t29, 0, 0, 0, 0, 0, t2, t1, t2, -t25 * t81, -t1, t8 * t23 - t54 * t25 + t5 * t28 + t50 * t29; 0, 0, 0, 0, t18 * t24 - t19 * t26 + (-t20 * t43 - t21 * t42) * pkin(2), t61 * t77, -t72 * t61, t39, -t76, 0, t58 * t46 + t56 * t67 + t74, (-qJD(2) * t24 - t58) * t44 + (t26 + t56) * t66, t57 * t67 + (-qJD(2) * t22 + t59) * t46 + t74, -t26 * t81 + t50, (-t26 - t57) * t66 + (-t73 * qJD(2) + t59) * t44, -t54 * t26 + t5 * t27 + t50 * t35 + t73 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t39, -t76, 0, t39, t54 * qJD(4) + t3 * t44 - t4 * t46; 0, 0, 0, 0, 0, -t62, t72 * t49, 0, 0, 0, -t12 * t69 + t52, t60 * qJD(4) - t12 * t68 - t75, (t30 * t46 - t44 * t8) * qJD(2) + t52, 0, (t30 * t44 + t46 * t8) * qJD(2) + (0.2e1 * qJD(5) - t60) * qJD(4) + t75, -t4 * pkin(4) + t3 * qJ(5) - t6 * t10 - t8 * t30 + t82 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t40 * t49 - t48, -t7 * qJD(4) + t8 * t69 + t4;];
tauc_reg = t11;
