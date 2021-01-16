% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP1
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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:08:52
% EndTime: 2021-01-15 20:08:54
% DurationCPUTime: 0.29s
% Computational Cost: add. (341->79), mult. (781->130), div. (0->0), fcn. (549->6), ass. (0->56)
t42 = sin(pkin(8));
t47 = cos(qJ(2));
t58 = pkin(1) * qJD(2);
t43 = cos(pkin(8));
t45 = sin(qJ(2));
t62 = t43 * t45;
t20 = (t42 * t47 + t62) * t58;
t44 = sin(qJ(4));
t37 = t44 * qJD(4);
t36 = pkin(4) * t37;
t10 = t36 + t20;
t34 = t47 * pkin(1) + pkin(2);
t63 = t42 * t45;
t48 = -pkin(1) * t63 + t43 * t34;
t18 = -pkin(3) - t48;
t46 = cos(qJ(4));
t64 = t46 * pkin(4);
t13 = t18 - t64;
t39 = t46 * qJD(4);
t65 = t10 * t44 + t13 * t39;
t61 = t18 * t39 + t20 * t44;
t60 = pkin(1) * t62 + t42 * t34;
t33 = -t43 * pkin(2) - pkin(3);
t27 = t33 - t64;
t41 = t44 ^ 2;
t59 = qJD(4) * t41 * pkin(4) + t27 * t39;
t19 = pkin(7) + t60;
t57 = qJ(5) + t19;
t32 = t42 * pkin(2) + pkin(7);
t56 = qJ(5) + t32;
t55 = t45 * t58;
t54 = t47 * t58;
t53 = pkin(4) * t39;
t4 = t13 * t37;
t22 = t27 * t37;
t52 = t33 * t37;
t51 = t33 * t39;
t50 = t44 * t39;
t49 = t18 * t37 - t20 * t46;
t40 = t46 * qJ(5);
t38 = t46 * qJD(5);
t29 = 0.2e1 * t50;
t26 = 0.2e1 * (t46 ^ 2 - t41) * qJD(4);
t25 = t46 * t32 + t40;
t24 = t56 * t44;
t21 = (t43 * t47 - t63) * t58;
t16 = t46 * t21;
t15 = -t44 * qJD(5) - t56 * t39;
t14 = t56 * t37 - t38;
t9 = t14 * t46;
t7 = t46 * t19 + t40;
t6 = t57 * t44;
t3 = (-qJD(5) - t21) * t44 - t57 * t39;
t2 = t57 * t37 - t16 - t38;
t1 = t2 * t46;
t5 = [0, 0, 0, 0, -0.2e1 * t55, -0.2e1 * t54, -0.2e1 * t48 * t20 + 0.2e1 * t60 * t21, t29, t26, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t61, -0.2e1 * t10 * t46 + 0.2e1 * t4, 0.2e1 * t65, -0.2e1 * t3 * t44 - 0.2e1 * t1 + 0.2e1 * (-t44 * t7 + t46 * t6) * qJD(4), 0.2e1 * t13 * t10 - 0.2e1 * t7 * t2 - 0.2e1 * t6 * t3; 0, 0, 0, 0, -t55, -t54, (-t20 * t43 + t21 * t42) * pkin(2), t29, t26, 0, 0, 0, t49 + t52, t51 + t61, t22 + t4 + (-t10 - t36) * t46, t59 + t65, -t1 - t9 + (-t15 - t3) * t44 + ((t24 + t6) * t46 + (-t25 - t7) * t44) * qJD(4), pkin(4) * t4 + t10 * t27 - t7 * t14 - t6 * t15 - t2 * t25 - t3 * t24; 0, 0, 0, 0, 0, 0, 0, t29, t26, 0, 0, 0, 0.2e1 * t52, 0.2e1 * t51, -0.2e1 * pkin(4) * t50 + 0.2e1 * t22, 0.2e1 * t59, -0.2e1 * t15 * t44 - 0.2e1 * t9 + 0.2e1 * (t24 * t46 - t25 * t44) * qJD(4), 0.2e1 * pkin(4) * t22 - 0.2e1 * t25 * t14 - 0.2e1 * t24 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t44 + t3 * t46 + (t44 * t6 + t46 * t7) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t44 + t15 * t46 + (t24 * t44 + t25 * t46) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, -t19 * t39 - t44 * t21, t19 * t37 - t16, t3, t2, -t53, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, -t32 * t39, t32 * t37, t15, t14, -t53, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t39, -t37, -t39, 0, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
