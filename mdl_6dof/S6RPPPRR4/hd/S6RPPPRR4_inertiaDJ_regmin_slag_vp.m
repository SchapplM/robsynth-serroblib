% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:52
% EndTime: 2019-03-09 01:35:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (240->69), mult. (531->147), div. (0->0), fcn. (406->6), ass. (0->61)
t27 = cos(qJ(6));
t21 = t27 ^ 2;
t25 = sin(qJ(6));
t63 = t25 ^ 2 - t21;
t42 = qJD(6) * t63;
t26 = sin(qJ(5));
t20 = t26 ^ 2;
t28 = cos(qJ(5));
t22 = t28 ^ 2;
t41 = (t20 - t22) * qJD(5);
t23 = sin(pkin(9));
t24 = cos(pkin(9));
t29 = -pkin(1) - pkin(2);
t43 = -t23 * qJ(2) + t24 * t29;
t39 = pkin(3) - t43;
t11 = pkin(7) + t39;
t18 = t26 * qJD(5);
t56 = t23 * qJD(2);
t35 = -t11 * t18 + t28 * t56;
t37 = pkin(5) * t28 + pkin(8) * t26;
t67 = qJD(6) * t37 + t35;
t66 = 0.2e1 * qJD(2);
t65 = t24 * t26;
t64 = t24 * qJ(2) + t23 * t29;
t61 = t20 + t22;
t60 = qJD(6) * t25;
t59 = qJD(6) * t26;
t58 = qJD(6) * t27;
t57 = qJD(6) * t28;
t55 = t24 * qJD(2);
t54 = t28 * qJD(5);
t53 = -0.2e1 * pkin(5) * qJD(6);
t12 = -qJ(4) + t64;
t52 = t25 * t57;
t51 = t27 * t57;
t50 = t22 * t56;
t49 = t25 * t54;
t48 = t25 * t58;
t47 = t27 * t18;
t46 = t27 * t54;
t45 = t26 * t54;
t44 = t24 * t54;
t40 = t27 * t45;
t14 = -qJD(4) + t55;
t38 = t37 * qJD(5) + t11 * t59 - t14;
t36 = -t26 * pkin(5) + t28 * pkin(8);
t34 = -t11 * t54 - t26 * t56;
t8 = t47 + t52;
t10 = -t25 * t18 + t51;
t33 = t22 * t60 + t40;
t32 = -t22 * t58 + t25 * t45;
t5 = t36 + t12;
t31 = -qJD(6) * t5 + t34;
t30 = t36 * qJD(5) - t11 * t57;
t9 = t26 * t58 + t49;
t7 = t25 * t59 - t46;
t4 = t25 * t44 + (-t25 * t23 + t27 * t65) * qJD(6);
t3 = -t27 * t44 + (t27 * t23 + t25 * t65) * qJD(6);
t2 = t31 * t25 - t38 * t27;
t1 = t38 * t25 + t31 * t27;
t6 = [0, 0, 0, 0, t66, qJ(2) * t66, 0.2e1 * t56, 0.2e1 * t55 (-t43 * t23 + t64 * t24) * t66, -0.2e1 * t56, -0.2e1 * t14, 0.2e1 * t12 * t14 + 0.2e1 * t39 * t56, -0.2e1 * t45, 0.2e1 * t41, 0, 0, 0, -0.2e1 * t12 * t54 - 0.2e1 * t14 * t26, 0.2e1 * t12 * t18 - 0.2e1 * t14 * t28, -0.2e1 * t21 * t45 - 0.2e1 * t22 * t48, 0.2e1 * t22 * t42 + 0.4e1 * t25 * t40, -0.2e1 * t26 * t52 - 0.2e1 * t27 * t41, 0.2e1 * t25 * t41 - 0.2e1 * t26 * t51, 0.2e1 * t45, -0.2e1 * t11 * t32 - 0.2e1 * t2 * t26 + 0.2e1 * t25 * t50 - 0.2e1 * t5 * t46, -0.2e1 * t1 * t26 - 0.2e1 * t11 * t33 + 0.2e1 * t27 * t50 + 0.2e1 * t5 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t14 - t55) * t23, 0, 0, 0, 0, 0, -t23 * t54, t23 * t18, 0, 0, 0, 0, 0, -t23 * t46 + t24 * t32 - t4 * t26, t23 * t49 + t24 * t33 + t3 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * t58, -t61 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t54, 0, t35, t34, t25 * t47 + t28 * t42, -t63 * t18 + 0.4e1 * t28 * t48, -t9, t7, 0, t30 * t25 + t67 * t27, -t67 * t25 + t30 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t18, t44, 0, 0, 0, 0, 0, t8 * t24, t10 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t18, 0, 0, 0, 0, 0, t7, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t54, 0, 0, 0, 0, 0, -t8, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t48, -0.2e1 * t42, 0, 0, 0, t25 * t53, t27 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, -t54, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t60, 0, -pkin(8) * t58, pkin(8) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
