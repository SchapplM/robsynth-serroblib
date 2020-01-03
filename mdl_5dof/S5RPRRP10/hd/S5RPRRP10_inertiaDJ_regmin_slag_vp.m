% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:06
% DurationCPUTime: 0.54s
% Computational Cost: add. (813->111), mult. (1901->206), div. (0->0), fcn. (1753->6), ass. (0->72)
t44 = sin(pkin(8));
t47 = sin(qJ(3));
t45 = cos(pkin(8));
t49 = cos(qJ(3));
t76 = t49 * t45;
t30 = t47 * t44 - t76;
t31 = t49 * t44 + t47 * t45;
t38 = -t45 * pkin(2) - pkin(1);
t19 = t30 * pkin(3) - t31 * pkin(7) + t38;
t75 = pkin(6) + qJ(2);
t32 = t75 * t44;
t33 = t75 * t45;
t22 = -t47 * t32 + t49 * t33;
t48 = cos(qJ(4));
t20 = t48 * t22;
t46 = sin(qJ(4));
t73 = t46 * t19 + t20;
t42 = t46 ^ 2;
t43 = t48 ^ 2;
t72 = t42 - t43;
t59 = t72 * qJD(4);
t26 = t31 * qJD(3);
t25 = t30 * qJD(3);
t53 = qJ(5) * t25 - qJD(5) * t31;
t18 = t26 * pkin(3) + t25 * pkin(7);
t77 = t49 * t32;
t9 = qJD(3) * t77 - qJD(2) * t76 + (qJD(2) * t44 + qJD(3) * t33) * t47;
t63 = t48 * t18 + t46 * t9;
t71 = qJ(5) * t31;
t1 = t26 * pkin(4) + t53 * t48 + (-t20 + (-t19 + t71) * t46) * qJD(4) + t63;
t69 = qJD(4) * t48;
t65 = t31 * t69;
t67 = t46 * t18 + t19 * t69 - t48 * t9;
t2 = -qJ(5) * t65 + (-qJD(4) * t22 + t53) * t46 + t67;
t61 = t48 * t19 - t46 * t22;
t5 = t30 * pkin(4) - t48 * t71 + t61;
t6 = -t46 * t71 + t73;
t85 = -t1 * t48 - t2 * t46 + (t46 * t5 - t48 * t6) * qJD(4);
t84 = 0.2e1 * t38;
t83 = t31 * t25;
t82 = t31 * t46;
t81 = t31 * t48;
t80 = t46 * t26;
t79 = t48 * t25;
t78 = t48 * t26;
t74 = -qJ(5) - pkin(7);
t70 = qJD(4) * t46;
t68 = -0.2e1 * pkin(3) * qJD(4);
t66 = pkin(4) * t70;
t64 = t46 * t69;
t62 = -0.4e1 * t46 * t81;
t21 = t47 * t33 + t77;
t60 = qJD(4) * t74;
t58 = 0.2e1 * (t44 ^ 2 + t45 ^ 2) * qJD(2);
t57 = pkin(3) * t25 - pkin(7) * t26;
t56 = pkin(3) * t31 + pkin(7) * t30;
t52 = -t46 * t25 + t65;
t51 = t31 * t70 + t79;
t50 = t30 * t69 + t80;
t10 = t31 * qJD(2) + t22 * qJD(3);
t39 = -t48 * pkin(4) - pkin(3);
t35 = t74 * t48;
t34 = t74 * t46;
t28 = t31 ^ 2;
t24 = -t46 * qJD(5) + t48 * t60;
t23 = t48 * qJD(5) + t46 * t60;
t17 = -t30 * t70 + t78;
t11 = pkin(4) * t82 + t21;
t7 = t52 * pkin(4) + t10;
t4 = -t73 * qJD(4) + t63;
t3 = t22 * t70 - t67;
t8 = [0, 0, 0, 0, 0, t58, qJ(2) * t58, -0.2e1 * t83, 0.2e1 * t25 * t30 - 0.2e1 * t31 * t26, 0, 0, 0, t26 * t84, -t25 * t84, -0.2e1 * t28 * t64 - 0.2e1 * t43 * t83, -t25 * t62 + 0.2e1 * t28 * t59, -0.2e1 * t51 * t30 + 0.2e1 * t31 * t78, -0.2e1 * t52 * t30 - 0.2e1 * t31 * t80, 0.2e1 * t30 * t26, 0.2e1 * t10 * t82 + 0.2e1 * t52 * t21 + 0.2e1 * t61 * t26 + 0.2e1 * t4 * t30, 0.2e1 * t10 * t81 - 0.2e1 * t51 * t21 - 0.2e1 * t73 * t26 + 0.2e1 * t3 * t30, -0.2e1 * (-t46 * t6 - t48 * t5) * t25 + 0.2e1 * t85 * t31, 0.2e1 * t5 * t1 + 0.2e1 * t11 * t7 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, 0, 0, 0, 0, t17, -t50, -(-t42 - t43) * t25, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t10, t9, -t31 * t59 - t46 * t79, qJD(4) * t62 + t72 * t25, t50, t17, 0, -t10 * t48 + t57 * t46 + (t21 * t46 - t56 * t48) * qJD(4), t10 * t46 + t57 * t48 + (t21 * t48 + t56 * t46) * qJD(4), (-t24 * t31 + t25 * t34 + t2 + (t31 * t35 - t5) * qJD(4)) * t48 + (-t23 * t31 - t25 * t35 - t1 + (t31 * t34 - t6) * qJD(4)) * t46, t1 * t34 + t11 * t66 - t2 * t35 + t6 * t23 + t5 * t24 + t7 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t23 + t48 * t24 + (-t34 * t46 - t35 * t48) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, -0.2e1 * t59, 0, 0, 0, t46 * t68, t48 * t68, 0.2e1 * t23 * t48 - 0.2e1 * t24 * t46 + 0.2e1 * (-t34 * t48 + t35 * t46) * qJD(4), -0.2e1 * t35 * t23 + 0.2e1 * t34 * t24 + 0.2e1 * t39 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t52, t26, t4, t3, t51 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t70, 0, -pkin(7) * t69, pkin(7) * t70, -pkin(4) * t69, t24 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
