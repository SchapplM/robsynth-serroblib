% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRP3
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
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:00
% EndTime: 2021-01-15 15:14:02
% DurationCPUTime: 0.40s
% Computational Cost: add. (492->114), mult. (1197->157), div. (0->0), fcn. (795->6), ass. (0->80)
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t46 = cos(qJ(2));
t65 = t46 * qJD(1);
t32 = qJD(2) * pkin(2) + t65;
t42 = cos(pkin(8));
t44 = sin(qJ(2));
t72 = qJD(1) * t44;
t34 = t42 * t72;
t41 = sin(pkin(8));
t19 = t41 * t32 + t34;
t12 = qJD(2) * pkin(6) + t19;
t64 = qJ(5) * qJD(2);
t56 = t12 + t64;
t52 = t56 * t45;
t7 = t43 * qJD(3) + t52;
t38 = t45 * qJD(3);
t6 = -t56 * t43 + t38;
t74 = qJD(4) * pkin(4);
t5 = t6 + t74;
t84 = t5 - t6;
t39 = t43 ^ 2;
t83 = pkin(4) * t39;
t82 = t42 * pkin(2);
t81 = t45 * pkin(4);
t48 = qJD(2) ^ 2;
t80 = t45 * t48;
t47 = qJD(4) ^ 2;
t79 = t47 * t43;
t78 = t47 * t45;
t23 = t41 * t65 + t34;
t33 = t41 * t72;
t25 = t42 * t65 - t33;
t67 = t43 * qJD(4);
t70 = qJD(2) * t45;
t77 = t23 * t70 + t25 * t67;
t40 = t45 ^ 2;
t76 = t39 - t40;
t75 = t39 + t40;
t36 = t41 * pkin(2) + pkin(6);
t73 = qJ(5) + t36;
t71 = qJD(2) * t23;
t18 = t42 * t32 - t33;
t11 = -qJD(2) * pkin(3) - t18;
t69 = t11 * qJD(2);
t66 = t45 * qJD(4);
t63 = qJD(2) * qJD(4);
t62 = qJD(2) * qJD(5);
t61 = 0.2e1 * t63;
t60 = -pkin(3) - t81;
t59 = t43 * t63;
t29 = t41 * t46 + t42 * t44;
t22 = t29 * qJD(2);
t20 = qJD(1) * t22;
t58 = -t36 * t47 - t20;
t57 = qJD(4) * t73;
t55 = t45 * t61;
t28 = t41 * t44 - t42 * t46;
t24 = t28 * qJD(2);
t21 = qJD(1) * t24;
t54 = qJD(4) * qJD(3) - t21;
t53 = t43 * t5 - t45 * t7;
t51 = qJD(4) * (qJD(2) * (-pkin(3) - t82) + t11);
t8 = t60 * qJD(2) + qJD(5) - t18;
t50 = (-qJD(5) - t8) * qJD(2) - t54;
t9 = t12 * t67;
t3 = -qJ(5) * t59 - t9 + (t54 + t62) * t45;
t4 = (t21 - t62) * t43 - t7 * qJD(4);
t49 = t3 * t45 - t4 * t43 + (-t43 * t7 - t45 * t5) * qJD(4);
t35 = pkin(4) * t59;
t31 = t60 - t82;
t27 = t73 * t45;
t26 = t73 * t43;
t16 = t25 * t66;
t14 = -t43 * qJD(5) - t45 * t57;
t13 = t45 * qJD(5) - t43 * t57;
t10 = t35 + t20;
t2 = t24 * t67 - t29 * t78 + (-t22 * t45 + t28 * t67) * qJD(2);
t1 = t24 * t66 + t29 * t79 + (t22 * t43 + t28 * t66) * qJD(2);
t15 = [0, 0, -t48 * t44, -t48 * t46, -t18 * t22 - t19 * t24 + t20 * t28 - t21 * t29, 0, 0, 0, 0, 0, t2, t1, t2, t1, -t75 * t24 * qJD(2), t10 * t28 + t8 * t22 + t53 * t24 + t49 * t29; 0, 0, 0, 0, t18 * t23 - t19 * t25 + (-t20 * t42 - t21 * t41) * pkin(2), t43 * t55, -t76 * t61, t78, -t79, 0, t43 * t51 + t58 * t45 + t77, t16 + t45 * t51 + (-t58 - t71) * t43, -t10 * t45 + (t14 + (t8 + (t31 - t81) * qJD(2)) * t43) * qJD(4) + t77, t16 + (t10 - t71) * t43 + (t45 * t8 - t13 + (t31 * t45 + t83) * qJD(2)) * qJD(4), (t13 * t45 - t14 * t43 - t75 * t25 + (t26 * t45 - t27 * t43) * qJD(4)) * qJD(2) + t49, t10 * t31 + t7 * t13 + t5 * t14 - t4 * t26 + t3 * t27 + (pkin(4) * t67 - t23) * t8 + t53 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78, -t79, -t78, 0, -t53 * qJD(4) + t3 * t43 + t4 * t45; 0, 0, 0, 0, 0, -t43 * t80, t76 * t48, 0, 0, 0, (t21 - t69) * t43, t9 + (-t43 * t12 + t38) * qJD(4) + (-t54 - t69) * t45, (t7 - t52) * qJD(4) + (pkin(4) * t80 + t50) * t43, -t48 * t83 + t9 + (t43 * t64 + t6) * qJD(4) + t50 * t45, (-t74 + t84) * t70, t84 * t7 + (-t8 * t43 * qJD(2) + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t59, t55, -t75 * t48, t35 + (t29 * qJD(1) + t53) * qJD(2);];
tauc_reg = t15;
