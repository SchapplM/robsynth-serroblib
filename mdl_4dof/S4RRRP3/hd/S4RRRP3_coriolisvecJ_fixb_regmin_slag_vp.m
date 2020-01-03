% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (413->96), mult. (782->134), div. (0->0), fcn. (347->4), ass. (0->75)
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t39 = cos(qJ(2));
t69 = pkin(1) * qJD(2);
t53 = qJD(1) * t69;
t50 = t39 * t53;
t26 = t38 * t50;
t33 = qJD(1) + qJD(2);
t37 = sin(qJ(2));
t70 = pkin(1) * qJD(1);
t59 = t37 * t70;
t21 = pkin(6) * t33 + t59;
t77 = t36 * t21;
t4 = t26 + (qJD(4) - t77) * qJD(3);
t25 = t36 * t50;
t67 = qJD(3) * t38;
t7 = t21 * t67 + t25;
t84 = t7 * t36 + t4 * t38;
t34 = t36 ^ 2;
t35 = t38 ^ 2;
t83 = t33 * (t34 + t35);
t40 = qJD(3) ^ 2;
t82 = pkin(6) * t40;
t81 = t39 * pkin(1);
t23 = -pkin(3) * t38 - qJ(4) * t36 - pkin(2);
t80 = t23 * t33;
t79 = t33 * t36;
t78 = t33 * t38;
t76 = t38 * t21;
t75 = t40 * t36;
t31 = t40 * t38;
t58 = t39 * t70;
t22 = -pkin(2) * t33 - t58;
t51 = t37 * t53;
t74 = t22 * t67 + t36 * t51;
t66 = qJD(3) * t39;
t55 = t36 * t66;
t73 = t55 * t70 + t59 * t78;
t72 = t34 - t35;
t68 = qJD(3) * t36;
t62 = qJD(3) * qJ(4);
t12 = t62 + t76;
t65 = t12 * qJD(3);
t64 = t36 * qJD(4);
t63 = -qJD(1) - t33;
t32 = t33 ^ 2;
t61 = t36 * t32 * t38;
t60 = t39 * t69;
t57 = t33 * t68;
t56 = t33 * t67;
t47 = pkin(3) * t36 - qJ(4) * t38;
t1 = t51 + (qJD(3) * t47 - t64) * t33;
t54 = -t1 - t82;
t52 = -qJD(3) * pkin(3) + qJD(4);
t49 = (-qJD(2) + t33) * t70;
t48 = t63 * t69;
t10 = t52 + t77;
t46 = t10 * t36 + t12 * t38;
t28 = pkin(1) * t37 + pkin(6);
t13 = pkin(3) * t68 - t38 * t62 - t64;
t9 = t37 * t69 + t13;
t45 = -t28 * t40 - t33 * t9 - t1;
t18 = t23 - t81;
t44 = t18 * t33 - t60;
t43 = t10 * t67 - t36 * t65 + t84;
t42 = -t37 * t79 + t38 * t66;
t41 = (t10 * t38 - t12 * t36) * qJD(3) + t84;
t29 = -pkin(2) - t81;
t20 = 0.2e1 * t36 * t56;
t16 = t22 * t68;
t14 = t47 * t33;
t11 = -0.2e1 * t72 * t33 * qJD(3);
t6 = -t58 + t80;
t3 = t6 * t68;
t2 = [0, 0, 0, 0, t37 * t48, t39 * t48, t20, t11, t31, -t75, 0, t29 * t57 - t28 * t31 + t16 + (t37 * t38 * t63 - t55) * t69, t28 * t75 + t29 * t56 - t42 * t69 + t74, t38 * t45 + t44 * t68 + t3, t60 * t83 + t43, t45 * t36 + (-t44 - t6) * t67, t1 * t18 + t28 * t41 + t46 * t60 + t6 * t9; 0, 0, 0, 0, t37 * t49, t39 * t49, t20, t11, t31, -t75, 0, -pkin(2) * t57 + t16 + (-t51 - t82) * t38 + t73, -pkin(2) * t56 + pkin(6) * t75 + t42 * t70 + t74, t23 * t57 + t3 + (-t13 * t33 + t54) * t38 + t73, -t58 * t83 + t43, (-t58 - t6 - t80) * t67 + ((-t13 + t59) * t33 + t54) * t36, t1 * t23 + t6 * t13 + (-t37 * t6 - t39 * t46) * t70 + t41 * pkin(6); 0, 0, 0, 0, 0, 0, -t61, t72 * t32, 0, 0, 0, -t22 * t79 - t25, -t22 * t78 - t26, -t25 + (t14 * t38 - t36 * t6) * t33, ((t12 - t62) * t36 + (-t10 + t52) * t38) * t33, 0.2e1 * qJD(3) * qJD(4) + t26 + (t14 * t36 + t38 * t6) * t33, -t10 * t76 - t7 * pkin(3) + t4 * qJ(4) - t6 * t14 + (qJD(4) + t77) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, -t32 * t34 - t40, t6 * t79 - t65 + t7;];
tauc_reg = t2;
