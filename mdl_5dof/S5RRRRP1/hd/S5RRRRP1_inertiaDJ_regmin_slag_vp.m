% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:02
% DurationCPUTime: 0.55s
% Computational Cost: add. (1257->97), mult. (2868->183), div. (0->0), fcn. (2622->6), ass. (0->66)
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t46 = cos(qJ(3));
t47 = cos(qJ(2));
t32 = t44 * t47 + t46 * t45;
t76 = pkin(6) + pkin(7);
t36 = t76 * t45;
t37 = t76 * t47;
t58 = -t46 * t36 - t44 * t37;
t15 = -t32 * pkin(8) + t58;
t51 = t32 * qJD(2);
t62 = qJD(2) * t76;
t33 = t45 * t62;
t34 = t47 * t62;
t59 = t46 * t33 + t44 * t34;
t11 = -pkin(8) * t51 + t15 * qJD(3) - t59;
t57 = t44 * t36 - t46 * t37;
t14 = t57 * qJD(3) + t44 * t33 - t46 * t34;
t56 = t44 * t45 - t46 * t47;
t25 = (-qJD(2) - qJD(3)) * t56;
t12 = -t25 * pkin(8) + t14;
t43 = sin(qJ(4));
t16 = -t56 * pkin(8) - t57;
t74 = cos(qJ(4));
t54 = t74 * t15 - t43 * t16;
t3 = -t54 * qJD(4) - t74 * t11 - t43 * t12;
t73 = pkin(2) * qJD(3);
t72 = qJD(4) * t43;
t71 = t45 * qJD(2);
t70 = t47 * qJD(2);
t69 = t43 * t44 * pkin(2);
t68 = -0.2e1 * pkin(1) * qJD(2);
t42 = pkin(2) * t71;
t67 = t44 * t73;
t66 = t46 * t73;
t65 = pkin(3) * t72;
t64 = t74 * pkin(3);
t41 = -t47 * pkin(2) - pkin(1);
t63 = t74 * t44;
t61 = t74 * qJD(4);
t60 = pkin(3) * t61;
t53 = -t43 * t15 - t74 * t16;
t52 = t41 * t32;
t50 = t74 * t56;
t40 = t46 * pkin(2) + pkin(3);
t18 = -t40 * t61 - t74 * t66 + (qJD(3) + qJD(4)) * t69;
t26 = t56 * pkin(3) + t41;
t24 = t74 * t32 - t43 * t56;
t4 = t53 * qJD(4) - t43 * t11 + t74 * t12;
t49 = -t32 * qJD(3) - t51;
t48 = (-t44 * t61 + (-t43 * t46 - t63) * qJD(3)) * pkin(2);
t17 = -t49 * pkin(3) + t42;
t39 = t64 + pkin(4);
t29 = pkin(2) * t63 + t43 * t40;
t27 = t74 * t40 + pkin(4) - t69;
t23 = t43 * t32 + t50;
t19 = -t40 * t72 + t48;
t13 = -t58 * qJD(3) + t59;
t9 = qJD(4) * t24 + t43 * t25 - t74 * t49;
t8 = qJD(4) * t50 - t74 * t25 + t32 * t72 - t43 * t49;
t7 = -t23 * qJ(5) - t53;
t6 = -t24 * qJ(5) + t54;
t5 = t9 * pkin(4) + t17;
t2 = t8 * qJ(5) - t24 * qJD(5) + t4;
t1 = -t9 * qJ(5) - t23 * qJD(5) - t3;
t10 = [0, 0, 0, 0.2e1 * t45 * t70, 0.2e1 * (-t45 ^ 2 + t47 ^ 2) * qJD(2), 0, 0, 0, t45 * t68, t47 * t68, 0.2e1 * t32 * t25, -0.2e1 * t25 * t56 + 0.2e1 * t32 * t49, 0, 0, 0, 0.2e1 * qJD(3) * t52 + 0.2e1 * (t45 * pkin(2) * t56 + t52) * qJD(2), 0.2e1 * t41 * t25 + 0.2e1 * t32 * t42, -0.2e1 * t24 * t8, 0.2e1 * t8 * t23 - 0.2e1 * t24 * t9, 0, 0, 0, 0.2e1 * t17 * t23 + 0.2e1 * t26 * t9, 0.2e1 * t17 * t24 - 0.2e1 * t26 * t8, -0.2e1 * t1 * t23 - 0.2e1 * t2 * t24 + 0.2e1 * t6 * t8 - 0.2e1 * t7 * t9, 0.2e1 * t7 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * (t23 * pkin(4) + t26) * t5; 0, 0, 0, 0, 0, t70, -t71, 0, -pkin(6) * t70, pkin(6) * t71, 0, 0, t25, t49, 0, t14, t13, 0, 0, -t8, -t9, 0, t4, t3, t18 * t23 - t19 * t24 + t27 * t8 - t29 * t9, t1 * t29 - t7 * t18 + t6 * t19 + t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t66, 0, 0, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t18, 0, -0.2e1 * t29 * t18 + 0.2e1 * t27 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t49, 0, t14, t13, 0, 0, -t8, -t9, 0, t4, t3, t39 * t8 + (-t43 * t9 + (-t74 * t23 + t24 * t43) * qJD(4)) * pkin(3), t2 * t39 + (t1 * t43 + (-t43 * t6 + t74 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t66, 0, 0, 0, 0, 0, (-pkin(3) - t40) * t72 + t48, -t60 + t18, 0, t19 * t39 + (-t18 * t43 + (-t27 * t43 + t74 * t29) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t65, -0.2e1 * t60, 0, 0.2e1 * (t64 - t39) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, t4, t3, pkin(4) * t8, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18, 0, t19 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t60, 0, -pkin(4) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
