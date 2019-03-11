% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:02
% EndTime: 2019-03-09 01:32:03
% DurationCPUTime: 0.42s
% Computational Cost: add. (444->77), mult. (986->148), div. (0->0), fcn. (940->8), ass. (0->59)
t43 = cos(qJ(6));
t38 = t43 ^ 2;
t42 = sin(qJ(6));
t63 = t42 ^ 2 - t38;
t52 = qJD(6) * t63;
t39 = sin(pkin(10));
t41 = cos(pkin(10));
t26 = (t39 ^ 2 + t41 ^ 2) * qJD(4);
t71 = sin(qJ(5));
t72 = cos(qJ(5));
t44 = -t71 * t39 + t72 * t41;
t21 = t44 ^ 2;
t76 = 2 * qJD(3);
t30 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t73 = -pkin(7) + t30;
t17 = t73 * t39;
t18 = t73 * t41;
t7 = t72 * t17 + t71 * t18;
t4 = t44 * qJD(4) + t7 * qJD(5);
t75 = t4 * t42;
t74 = t4 * t43;
t53 = qJD(5) * t71;
t54 = qJD(5) * t72;
t19 = -t39 * t54 - t41 * t53;
t70 = t44 * t19;
t20 = -t39 * t53 + t41 * t54;
t23 = t72 * t39 + t71 * t41;
t69 = t23 * t20;
t68 = t42 * t19;
t67 = t43 * t19;
t66 = t43 * t20;
t65 = t23 * t67 + t44 * t66;
t31 = sin(pkin(9)) * pkin(1) + qJ(3);
t62 = qJD(6) * t42;
t61 = qJD(6) * t43;
t60 = t31 * qJD(3);
t59 = -0.2e1 * pkin(5) * qJD(6);
t25 = t39 * pkin(4) + t31;
t58 = t44 * t62;
t57 = t42 * t61;
t56 = t23 ^ 2 + t21;
t55 = -0.4e1 * t44 * t42 * t43;
t51 = -pkin(5) * t19 - pkin(8) * t20;
t50 = -pkin(5) * t44 - pkin(8) * t23;
t5 = t23 * pkin(5) - pkin(8) * t44 + t25;
t49 = t42 * t7 - t43 * t5;
t48 = t42 * t5 + t43 * t7;
t47 = -t19 * t23 - t20 * t44;
t11 = -t44 * t61 - t68;
t46 = -t58 + t67;
t10 = t42 * t20 + t23 * t61;
t45 = -0.2e1 * t69 - 0.2e1 * t70;
t12 = t20 * pkin(5) - t19 * pkin(8) + qJD(3);
t8 = t23 * t62 - t66;
t6 = t71 * t17 - t72 * t18;
t3 = t23 * qJD(4) + t17 * t53 - t18 * t54;
t2 = -t48 * qJD(6) + t43 * t12 + t42 * t3;
t1 = t49 * qJD(6) - t42 * t12 + t43 * t3;
t9 = [0, 0, 0, 0, 0, t76, 0.2e1 * t60, t39 * t76, t41 * t76, 0.2e1 * t26, -0.2e1 * t30 * t26 + 0.2e1 * t60, 0.2e1 * t70, 0.2e1 * t47, 0, 0, 0, 0.2e1 * qJD(3) * t23 + 0.2e1 * t25 * t20, 0.2e1 * qJD(3) * t44 + 0.2e1 * t25 * t19, -0.2e1 * t21 * t57 + 0.2e1 * t38 * t70, t19 * t55 + 0.2e1 * t21 * t52, -0.2e1 * t23 * t58 + 0.2e1 * t65, -0.2e1 * t10 * t44 - 0.2e1 * t23 * t68, 0.2e1 * t69, 0.2e1 * t2 * t23 - 0.2e1 * t49 * t20 + 0.2e1 * t6 * t68 - 0.2e1 * (-t6 * t61 - t75) * t44, 0.2e1 * t1 * t23 - 0.2e1 * t48 * t20 + 0.2e1 * t6 * t67 - 0.2e1 * (t6 * t62 - t74) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t43 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t45 - t56 * t61, t43 * t45 + t56 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, t20, t19, 0, 0, 0, 0, 0, -t8, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t20, 0, -t4, t3, t42 * t67 - t44 * t52, qJD(6) * t55 - t63 * t19, t10, -t8, 0, -t74 + t51 * t42 + (t42 * t6 + t50 * t43) * qJD(6), t75 + t51 * t43 + (-t50 * t42 + t43 * t6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t19, 0, 0, 0, 0, 0, t8, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t20, 0, 0, 0, 0, 0, t46, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t57, -0.2e1 * t52, 0, 0, 0, t42 * t59, t43 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t11, t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t62, 0, -pkin(8) * t61, pkin(8) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
