% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:43
% EndTime: 2019-03-09 01:51:44
% DurationCPUTime: 0.47s
% Computational Cost: add. (229->82), mult. (501->149), div. (0->0), fcn. (328->4), ass. (0->68)
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t65 = qJ(5) * t31;
t70 = pkin(4) + pkin(8);
t74 = t70 * t33 + t65;
t64 = t33 * qJ(5);
t73 = -t31 * pkin(4) + t64;
t25 = t31 ^ 2;
t27 = t33 ^ 2;
t46 = (t25 - t27) * qJD(4);
t30 = sin(qJ(6));
t24 = t30 ^ 2;
t32 = cos(qJ(6));
t68 = -t32 ^ 2 + t24;
t45 = t68 * qJD(6);
t56 = t31 * qJD(5);
t28 = -pkin(7) + qJ(2);
t69 = pkin(5) - t28;
t16 = t69 * t31;
t58 = t16 * qJD(6);
t72 = (-t31 * t70 + t64) * qJD(4) + t56 + t58;
t71 = 0.2e1 * qJD(5);
t29 = pkin(1) + qJ(3);
t66 = t25 + t27;
t63 = qJD(4) * t16;
t62 = qJD(6) * t30;
t61 = qJD(6) * t32;
t60 = qJD(6) * t33;
t59 = qJD(6) * t70;
t57 = t31 * qJD(2);
t21 = t31 * qJD(4);
t55 = t33 * qJD(2);
t22 = t33 * qJD(4);
t54 = qJ(2) * qJD(2);
t53 = qJ(5) * qJD(6);
t52 = t30 * t60;
t51 = t32 * t60;
t50 = t32 * t22;
t49 = t30 * t61;
t48 = t31 * t22;
t47 = qJD(4) * t69;
t44 = qJD(6) * t66;
t18 = t66 * qJD(2);
t43 = t30 * t50;
t42 = -t33 * qJD(5) + qJD(3);
t17 = t69 * t33;
t15 = t29 - t73;
t8 = t31 * pkin(8) + t15;
t41 = t32 * t17 - t30 * t8;
t40 = t30 * t17 + t32 * t8;
t39 = pkin(4) * t33 + t65;
t5 = -t33 * t47 + t57;
t36 = t74 * qJD(6) + t5;
t7 = t73 * qJD(4) + t56;
t35 = 0.2e1 * qJD(2);
t19 = -0.2e1 * t48;
t14 = t28 * t22 + t57;
t13 = t28 * t21 - t55;
t12 = -t30 * t21 + t51;
t11 = t30 * t22 + t31 * t61;
t10 = t32 * t21 + t52;
t9 = -t31 * t62 + t50;
t6 = t39 * qJD(4) + t42;
t4 = -t31 * t47 - t55;
t3 = t74 * qJD(4) + t42;
t2 = -t40 * qJD(6) - t30 * t3 + t32 * t4;
t1 = -t41 * qJD(6) - t32 * t3 - t30 * t4;
t20 = [0, 0, 0, 0, t35, 0.2e1 * t54, t35, 0.2e1 * qJD(3), 0.2e1 * t29 * qJD(3) + 0.2e1 * t54, t19, 0.2e1 * t46, 0, 0, 0, 0.2e1 * qJD(3) * t31 + 0.2e1 * t29 * t22, 0.2e1 * qJD(3) * t33 - 0.2e1 * t29 * t21, -0.2e1 * t18, -0.2e1 * t15 * t22 - 0.2e1 * t6 * t31, 0.2e1 * t15 * t21 - 0.2e1 * t6 * t33, 0.2e1 * t15 * t6 + 0.2e1 * t28 * t18, 0.2e1 * t24 * t48 + 0.2e1 * t25 * t49, -0.2e1 * t25 * t45 + 0.4e1 * t31 * t43, -0.2e1 * t30 * t46 + 0.2e1 * t31 * t51, -0.2e1 * t31 * t52 - 0.2e1 * t32 * t46, t19, 0.2e1 * (t32 * t63 + t2) * t33 + 0.2e1 * (-t41 * qJD(4) - t30 * t58 - t5 * t32) * t31, 0.2e1 * (-t30 * t63 + t1) * t33 + 0.2e1 * (t40 * qJD(4) + t5 * t30 - t32 * t58) * t31; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t22, t21, 0, t22, -t21, -t6, 0, 0, 0, 0, 0, t12, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, t30 * t44, t32 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, -t13, -t14, -t7, t13, t14, t39 * qJD(2) + t7 * t28, -t31 * t45 + t43, -t68 * t22 - 0.4e1 * t31 * t49, -t10, -t12, 0, t36 * t30 - t72 * t32, t72 * t30 + t36 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, t21, t22, t7, 0, 0, 0, 0, 0, t11, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, qJ(5) * t71, -0.2e1 * t49, 0.2e1 * t45, 0, 0, 0, 0.2e1 * qJD(5) * t30 + 0.2e1 * t32 * t53, 0.2e1 * qJD(5) * t32 - 0.2e1 * t30 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, t13, 0, 0, 0, 0, 0, -t10, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t9, -t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61, 0, t30 * t59, t32 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t20;
