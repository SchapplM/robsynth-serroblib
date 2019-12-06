% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:22
% EndTime: 2019-12-05 18:20:23
% DurationCPUTime: 0.26s
% Computational Cost: add. (240->58), mult. (679->110), div. (0->0), fcn. (470->8), ass. (0->61)
t74 = 2 * qJD(5);
t38 = sin(pkin(8));
t40 = cos(pkin(8));
t44 = cos(qJ(2));
t59 = pkin(1) * qJD(2);
t54 = t44 * t59;
t42 = sin(qJ(2));
t55 = t42 * t59;
t17 = -t38 * t55 + t40 * t54;
t13 = qJD(4) + t17;
t37 = sin(pkin(9));
t35 = t37 ^ 2;
t11 = t35 * t13;
t32 = t44 * pkin(1) + pkin(2);
t65 = t40 * t42;
t62 = pkin(1) * t65 + t38 * t32;
t15 = qJ(4) + t62;
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t56 = qJD(5) * t43;
t52 = t35 * t56;
t73 = t41 * t11 + t15 * t52;
t72 = t40 * pkin(2);
t16 = (t38 * t44 + t65) * t59;
t39 = cos(pkin(9));
t66 = t39 * t43;
t67 = t39 * t41;
t45 = -t39 * pkin(4) - t37 * pkin(7) - pkin(3);
t48 = -t38 * t42 * pkin(1) + t40 * t32;
t7 = t45 - t48;
t2 = -t13 * t66 - t41 * t16 + (t15 * t67 - t43 * t7) * qJD(5);
t71 = t43 * t11 - t2 * t39;
t33 = t35 * qJD(4);
t18 = t45 - t72;
t31 = t38 * pkin(2) + qJ(4);
t58 = qJD(4) * t39;
t5 = -t43 * t58 + (-t18 * t43 + t31 * t67) * qJD(5);
t70 = t43 * t33 - t5 * t39;
t69 = t16 * t37;
t68 = t16 * t39;
t36 = t39 ^ 2;
t64 = t36 * t13 + t11;
t63 = t31 * t52 + t41 * t33;
t61 = t36 * qJD(4) + t33;
t60 = t35 + t36;
t57 = qJD(5) * t41;
t53 = t35 * t57;
t51 = t37 * t57;
t50 = t37 * t56;
t49 = t13 * t60;
t47 = qJD(4) * t60;
t46 = t37 * t39 * t74;
t28 = t39 * t56;
t27 = t39 * t57;
t22 = -0.2e1 * t41 * t52;
t21 = t43 * t46;
t20 = t41 * t46;
t14 = (t41 ^ 2 - t43 ^ 2) * t35 * t74;
t6 = -t41 * t58 + (-t18 * t41 - t31 * t66) * qJD(5);
t3 = -t13 * t67 + t43 * t16 + (-t15 * t66 - t41 * t7) * qJD(5);
t1 = [0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t54, -0.2e1 * t48 * t16 + 0.2e1 * t62 * t17, -0.2e1 * t68, 0.2e1 * t69, 0.2e1 * t64, 0.2e1 * (-pkin(3) - t48) * t16 + 0.2e1 * t15 * t49, t22, t14, t20, t21, 0, -0.2e1 * t3 * t39 + 0.2e1 * t73, -0.2e1 * t15 * t53 + 0.2e1 * t71; 0, 0, 0, 0, -t55, -t54, (-t16 * t40 + t17 * t38) * pkin(2), -t68, t69, t61 + t64, t16 * (-pkin(3) - t72) + t31 * t49 + t15 * t47, t22, t14, t20, t21, 0, (-t3 - t6) * t39 + t63 + t73, (-t15 - t31) * t53 + t70 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t61, 0.2e1 * t31 * t47, t22, t14, t20, t21, 0, -0.2e1 * t6 * t39 + 0.2e1 * t63, -0.2e1 * t31 * t53 + 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
