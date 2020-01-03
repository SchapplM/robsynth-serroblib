% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:28
% EndTime: 2019-12-31 20:58:29
% DurationCPUTime: 0.39s
% Computational Cost: add. (656->95), mult. (1536->155), div. (0->0), fcn. (1239->4), ass. (0->56)
t69 = qJD(2) + qJD(3);
t48 = 2 * qJD(4);
t47 = -pkin(3) - pkin(4);
t68 = -pkin(7) - pkin(6);
t67 = cos(qJ(3));
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t66 = t44 * t45;
t46 = cos(qJ(2));
t65 = t44 * t46;
t29 = t68 * t46;
t15 = -t67 * t29 + t68 * t66;
t56 = t67 * qJD(3);
t42 = pkin(2) * t56;
t31 = t42 + qJD(4);
t37 = pkin(2) * t44 + qJ(4);
t64 = t31 * qJ(4) + t37 * qJD(4);
t63 = qJD(3) * t44;
t62 = t45 * qJD(2);
t61 = t46 * qJD(2);
t60 = -0.2e1 * pkin(1) * qJD(2);
t59 = pkin(2) * t62;
t41 = pkin(2) * t63;
t40 = -t46 * pkin(2) - pkin(1);
t58 = t67 * t46;
t57 = t68 * qJD(2);
t55 = t68 * t67;
t39 = -pkin(2) * t67 - pkin(3);
t54 = t45 * t55;
t23 = t45 * t67 + t65;
t13 = t69 * t23;
t22 = -t58 + t66;
t53 = -qJ(4) * t13 - qJD(4) * t22;
t52 = qJD(2) * t55;
t51 = t23 * qJ(4) - t40;
t5 = -qJD(3) * t54 - t29 * t63 - t45 * t52 - t57 * t65;
t14 = -t44 * t29 - t54;
t12 = -qJD(2) * t58 - t46 * t56 + t66 * t69;
t50 = -t12 * qJ(4) + t23 * qJD(4) - t59;
t49 = t13 * t37 + t22 * t31 - t23 * t41;
t6 = -t29 * t56 - t46 * t52 + (qJD(3) * t68 + t57) * t66;
t43 = qJ(4) * t48;
t36 = -0.2e1 * t41;
t35 = -pkin(4) + t39;
t32 = t48 + t42;
t27 = 0.2e1 * t31;
t20 = t37 * t31;
t10 = pkin(3) * t22 - t51;
t9 = qJ(5) * t22 + t15;
t8 = -t23 * qJ(5) + t14;
t7 = t22 * t47 + t51;
t4 = pkin(3) * t13 - t50;
t3 = t12 * qJ(5) - t23 * qJD(5) + t6;
t2 = qJ(5) * t13 + qJD(5) * t22 - t5;
t1 = t13 * t47 + t50;
t11 = [0, 0, 0, 0.2e1 * t45 * t61, 0.2e1 * (-t45 ^ 2 + t46 ^ 2) * qJD(2), 0, 0, 0, t45 * t60, t46 * t60, -0.2e1 * t23 * t12, 0.2e1 * t12 * t22 - 0.2e1 * t13 * t23, 0, 0, 0, 0.2e1 * t13 * t40 + 0.2e1 * t22 * t59, -0.2e1 * t12 * t40 + 0.2e1 * t23 * t59, 0.2e1 * t10 * t13 + 0.2e1 * t22 * t4, -0.2e1 * t12 * t14 - 0.2e1 * t13 * t15 + 0.2e1 * t22 * t5 + 0.2e1 * t23 * t6, 0.2e1 * t10 * t12 - 0.2e1 * t23 * t4, 0.2e1 * t10 * t4 + 0.2e1 * t14 * t6 - 0.2e1 * t15 * t5, -0.2e1 * t1 * t22 - 0.2e1 * t13 * t7, 0.2e1 * t1 * t23 - 0.2e1 * t12 * t7, 0.2e1 * t12 * t8 + 0.2e1 * t13 * t9 + 0.2e1 * t2 * t22 - 0.2e1 * t23 * t3, 0.2e1 * t1 * t7 + 0.2e1 * t2 * t9 + 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, t61, -t62, 0, -pkin(6) * t61, pkin(6) * t62, 0, 0, -t12, -t13, 0, -t6, t5, -t6, -t12 * t39 - t49, -t5, t14 * t41 + t15 * t31 - t37 * t5 + t39 * t6, -t3, t2, t12 * t35 + t49, t2 * t37 + t3 * t35 + t31 * t9 + t41 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -0.2e1 * t42, t36, 0, t27, 0.2e1 * t39 * t41 + 0.2e1 * t20, t36, t27, 0, 0.2e1 * t35 * t41 + 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, -t6, t5, -t6, pkin(3) * t12 + t53, -t5, -pkin(3) * t6 - qJ(4) * t5 + qJD(4) * t15, -t3, t2, t12 * t47 - t53, qJ(4) * t2 + qJD(4) * t9 + t3 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, -t41, 0, t32, -pkin(3) * t41 + t64, -t41, t32, 0, t41 * t47 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t43, 0, t48, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, t6, 0, 0, t12, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
