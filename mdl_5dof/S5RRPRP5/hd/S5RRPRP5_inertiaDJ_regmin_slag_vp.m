% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:18
% EndTime: 2021-01-15 20:18:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (1020->96), mult. (2300->171), div. (0->0), fcn. (2098->6), ass. (0->58)
t44 = sin(qJ(2));
t45 = cos(qJ(2));
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t54 = t72 * t44 - t73 * t45;
t33 = t73 * t44 + t72 * t45;
t46 = 2 * qJD(5);
t76 = cos(qJ(4));
t75 = -qJ(3) - pkin(6);
t19 = t54 * t75;
t59 = qJD(2) * t72;
t60 = qJD(2) * t73;
t74 = t44 * t60 + t45 * t59;
t43 = sin(qJ(4));
t71 = qJD(4) * t43;
t70 = t44 * qJD(2);
t69 = t45 * qJD(2);
t68 = -0.2e1 * pkin(1) * qJD(2);
t42 = pkin(2) * t70;
t41 = -t45 * pkin(2) - pkin(1);
t67 = t72 * pkin(2);
t66 = qJD(4) * t76;
t61 = qJD(2) * t75;
t58 = t43 * t67;
t20 = t74 * pkin(3) + t42;
t57 = -t44 * t59 + t45 * t60;
t40 = t73 * pkin(2) + pkin(3);
t24 = qJD(4) * t58 - t40 * t66;
t56 = -t44 * qJD(3) + t45 * t61;
t55 = t45 * qJD(3) + t44 * t61;
t53 = t43 * t40 + t76 * t67;
t52 = t76 * t54;
t18 = t33 * t75;
t23 = t54 * pkin(3) + t41;
t17 = t76 * t33 - t43 * t54;
t51 = -t33 * pkin(7) + t18;
t50 = t43 * t51;
t49 = t76 * t51;
t14 = t73 * t55 + t72 * t56;
t13 = -t72 * t55 + t73 * t56;
t48 = -t74 * pkin(7) + t14;
t47 = t57 * pkin(7) - t13;
t30 = -t76 * t40 - pkin(4) + t58;
t29 = qJ(5) + t53;
t25 = t53 * qJD(4);
t22 = qJD(5) - t24;
t21 = 0.2e1 * t25;
t16 = t43 * t33 + t52;
t15 = -t54 * pkin(7) + t19;
t8 = t17 * qJD(4) + t43 * t57 + t76 * t74;
t7 = qJD(4) * t52 + t33 * t71 + t43 * t74 - t76 * t57;
t6 = t76 * t15 + t50;
t5 = t43 * t15 - t49;
t4 = t16 * pkin(4) - t17 * qJ(5) + t23;
t3 = t8 * pkin(4) + t7 * qJ(5) - t17 * qJD(5) + t20;
t2 = qJD(4) * t50 + t15 * t66 + t43 * t48 + t76 * t47;
t1 = -qJD(4) * t49 + t15 * t71 + t43 * t47 - t76 * t48;
t9 = [0, 0, 0, 0.2e1 * t44 * t69, 0.2e1 * (-t44 ^ 2 + t45 ^ 2) * qJD(2), 0, 0, 0, t44 * t68, t45 * t68, 0.2e1 * t41 * t74 + 0.2e1 * t54 * t42, 0.2e1 * t33 * t42 + 0.2e1 * t41 * t57, -0.2e1 * t13 * t33 - 0.2e1 * t14 * t54 - 0.2e1 * t18 * t57 - 0.2e1 * t19 * t74, 0.2e1 * t18 * t13 + 0.2e1 * t19 * t14 + 0.2e1 * t41 * t42, -0.2e1 * t17 * t7, 0.2e1 * t7 * t16 - 0.2e1 * t17 * t8, 0, 0, 0, 0.2e1 * t20 * t16 + 0.2e1 * t23 * t8, 0.2e1 * t20 * t17 - 0.2e1 * t23 * t7, 0.2e1 * t3 * t16 + 0.2e1 * t4 * t8, 0.2e1 * t1 * t16 + 0.2e1 * t2 * t17 - 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t3 * t17 + 0.2e1 * t4 * t7, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, t69, -t70, 0, -pkin(6) * t69, pkin(6) * t70, t13, -t14, (-t73 * t57 - t72 * t74) * pkin(2), (t73 * t13 + t72 * t14) * pkin(2), 0, 0, -t7, -t8, 0, -t2, t1, -t2, -t22 * t16 + t25 * t17 - t29 * t8 - t30 * t7, -t1, -t1 * t29 + t2 * t30 + t6 * t22 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0.2e1 * t24, -t21, 0, 0.2e1 * t22, 0.2e1 * t29 * t22 + 0.2e1 * t30 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t57, 0, t42, 0, 0, 0, 0, 0, t8, -t7, t8, 0, t7, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t2, t1, -t2, pkin(4) * t7 - t8 * qJ(5) - t16 * qJD(5), -t1, -t2 * pkin(4) - t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, -t25, 0, t46 - t24, -t25 * pkin(4) + t22 * qJ(5) + t29 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, qJ(5) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
