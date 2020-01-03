% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:15
% EndTime: 2019-12-31 19:39:16
% DurationCPUTime: 0.37s
% Computational Cost: add. (415->79), mult. (964->157), div. (0->0), fcn. (816->6), ass. (0->55)
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t74 = -t59 * pkin(2) - t57 * qJ(3);
t54 = sin(pkin(8));
t55 = cos(pkin(8));
t34 = t57 * t54 + t59 * t55;
t35 = -t59 * t54 + t57 * t55;
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t14 = t58 * t34 + t56 * t35;
t73 = 2 * qJD(3);
t60 = -pkin(2) - pkin(3);
t70 = pkin(6) - qJ(4);
t68 = t57 * qJD(2);
t27 = -t59 * qJD(4) - t70 * t68;
t44 = t70 * t59;
t62 = qJD(2) * t44 - t57 * qJD(4);
t10 = t55 * t27 + t54 * t62;
t43 = t70 * t57;
t18 = t54 * t43 + t55 * t44;
t50 = t59 * qJD(2);
t69 = qJ(3) * t50 + t57 * qJD(3);
t67 = -0.2e1 * pkin(1) * qJD(2);
t42 = -pkin(1) + t74;
t66 = pkin(6) * t68;
t65 = pkin(6) * t50;
t9 = t54 * t27 - t55 * t62;
t17 = t55 * t43 - t54 * t44;
t33 = t59 * pkin(3) - t42;
t40 = -t54 * qJ(3) + t55 * t60;
t15 = -t56 * t34 + t58 * t35;
t64 = t54 * t58 + t55 * t56;
t63 = t54 * t56 - t55 * t58;
t19 = t60 * t68 + t69;
t61 = t74 * qJD(2) + t59 * qJD(3);
t41 = t55 * qJ(3) + t54 * t60;
t39 = -pkin(4) + t40;
t32 = t34 * qJD(2);
t31 = t54 * t50 - t55 * t68;
t30 = t64 * qJD(5);
t29 = t63 * qJD(5);
t28 = pkin(2) * t68 - t69;
t16 = t34 * pkin(4) + t33;
t13 = t31 * pkin(4) + t19;
t12 = -t34 * pkin(7) + t18;
t11 = -t35 * pkin(7) + t17;
t8 = (t39 * t56 + t41 * t58) * qJD(5) + t64 * qJD(3);
t7 = (-t39 * t58 + t41 * t56) * qJD(5) + t63 * qJD(3);
t6 = -t31 * pkin(7) + t10;
t5 = -t32 * pkin(7) - t9;
t4 = t15 * qJD(5) + t58 * t31 + t56 * t32;
t3 = t14 * qJD(5) + t56 * t31 - t58 * t32;
t2 = -t58 * t5 + t56 * t6 + (t11 * t56 + t12 * t58) * qJD(5);
t1 = -t56 * t5 - t58 * t6 + (-t11 * t58 + t12 * t56) * qJD(5);
t20 = [0, 0, 0, 0.2e1 * t57 * t50, 0.2e1 * (-t57 ^ 2 + t59 ^ 2) * qJD(2), 0, 0, 0, t57 * t67, t59 * t67, -0.2e1 * t28 * t59 + 0.2e1 * t42 * t68, 0, -0.2e1 * t28 * t57 - 0.2e1 * t42 * t50, 0.2e1 * t42 * t28, 0.2e1 * t19 * t34 + 0.2e1 * t33 * t31, 0.2e1 * t19 * t35 + 0.2e1 * t33 * t32, -0.2e1 * t10 * t34 - 0.2e1 * t17 * t32 - 0.2e1 * t18 * t31 + 0.2e1 * t9 * t35, 0.2e1 * t18 * t10 - 0.2e1 * t17 * t9 + 0.2e1 * t33 * t19, -0.2e1 * t15 * t3, 0.2e1 * t3 * t14 - 0.2e1 * t15 * t4, 0, 0, 0, 0.2e1 * t13 * t14 + 0.2e1 * t16 * t4, 0.2e1 * t13 * t15 - 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, t50, -t68, 0, -t65, t66, -t65, t61, -t66, t61 * pkin(6), t9, t10, -t41 * t31 - t40 * t32 + (-t34 * t55 + t35 * t54) * qJD(3), t10 * t41 - t9 * t40 + (-t17 * t54 + t18 * t55) * qJD(3), 0, 0, t3, t4, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, qJ(3) * t73, t54 * t73, t55 * t73, 0, (-t40 * t54 + t41 * t55) * t73, 0, 0, 0, 0, 0, 0.2e1 * t8, -0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t65, 0, 0, -t54 * t31 - t55 * t32, t10 * t54 - t9 * t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, t19, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
