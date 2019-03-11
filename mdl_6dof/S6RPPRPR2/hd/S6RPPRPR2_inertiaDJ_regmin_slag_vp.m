% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:37
% EndTime: 2019-03-09 01:42:38
% DurationCPUTime: 0.56s
% Computational Cost: add. (637->97), mult. (1442->173), div. (0->0), fcn. (1363->8), ass. (0->63)
t42 = cos(pkin(10));
t79 = cos(qJ(4));
t61 = qJD(4) * t79;
t41 = sin(pkin(10));
t45 = sin(qJ(4));
t74 = t45 * t41;
t25 = qJD(4) * t74 - t42 * t61;
t31 = t79 * t41 + t45 * t42;
t44 = sin(qJ(6));
t46 = cos(qJ(6));
t69 = qJD(6) * t46;
t16 = -t44 * t25 + t31 * t69;
t39 = t44 ^ 2;
t71 = -t46 ^ 2 + t39;
t60 = t71 * qJD(6);
t35 = sin(pkin(9)) * pkin(1) + qJ(3);
t80 = pkin(7) + t35;
t27 = t80 * t41;
t28 = t80 * t42;
t62 = t79 * t42;
t6 = (qJD(3) * t41 + qJD(4) * t28) * t45 - qJD(3) * t62 + t27 * t61;
t83 = 2 * qJD(5);
t82 = pkin(4) + pkin(8);
t26 = t31 * qJD(4);
t3 = -t26 * pkin(5) - t6;
t30 = -t62 + t74;
t81 = t3 * t30;
t78 = t30 * t26;
t77 = t31 * t25;
t75 = t44 * t26;
t73 = t46 * t26;
t70 = qJD(6) * t44;
t68 = qJD(6) * t82;
t67 = qJ(5) * qJD(6);
t66 = t44 * t73;
t63 = t44 * t69;
t59 = 0.2e1 * (t41 ^ 2 + t42 ^ 2) * qJD(3);
t32 = -cos(pkin(9)) * pkin(1) - pkin(2) - t42 * pkin(3);
t50 = -t31 * qJ(5) + t32;
t8 = t82 * t30 + t50;
t17 = t79 * t27 + t45 * t28;
t9 = t31 * pkin(5) + t17;
t58 = t44 * t9 + t46 * t8;
t57 = t44 * t8 - t46 * t9;
t56 = t25 * qJ(5) - t31 * qJD(5);
t55 = -qJ(5) * t26 - qJD(5) * t30;
t53 = t46 * t25 + t31 * t70;
t15 = t30 * t69 + t75;
t52 = t30 * t70 - t73;
t51 = -t45 * t27 + t79 * t28;
t11 = t26 * pkin(4) + t56;
t49 = t3 + (qJ(5) * t30 + t31 * t82) * qJD(6);
t10 = -t30 * pkin(5) + t51;
t48 = -qJD(6) * t10 - t25 * t82 - t55;
t7 = t31 * qJD(3) + t51 * qJD(4);
t29 = t30 ^ 2;
t22 = -0.2e1 * t77;
t12 = t30 * pkin(4) + t50;
t5 = t82 * t26 + t56;
t4 = -t25 * pkin(5) + t7;
t2 = -t58 * qJD(6) + t46 * t4 - t44 * t5;
t1 = t57 * qJD(6) - t44 * t4 - t46 * t5;
t13 = [0, 0, 0, 0, 0, 0, t59, t35 * t59, t22, 0.2e1 * t25 * t30 - 0.2e1 * t26 * t31, 0, 0, 0, 0.2e1 * t32 * t26, -0.2e1 * t32 * t25, -0.2e1 * t17 * t25 - 0.2e1 * t26 * t51 + 0.2e1 * t6 * t30 + 0.2e1 * t7 * t31, -0.2e1 * t11 * t30 - 0.2e1 * t12 * t26, -0.2e1 * t11 * t31 + 0.2e1 * t12 * t25, 0.2e1 * t12 * t11 + 0.2e1 * t17 * t7 - 0.2e1 * t51 * t6, 0.2e1 * t29 * t63 + 0.2e1 * t39 * t78, -0.2e1 * t29 * t60 + 0.4e1 * t30 * t66, 0.2e1 * t16 * t30 + 0.2e1 * t31 * t75, -0.2e1 * t53 * t30 + 0.2e1 * t31 * t73, t22, 0.2e1 * t52 * t10 + 0.2e1 * t2 * t31 + 0.2e1 * t57 * t25 - 0.2e1 * t46 * t81, 0.2e1 * t1 * t31 + 0.2e1 * t15 * t10 + 0.2e1 * t58 * t25 + 0.2e1 * t44 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t26 - t25 * t51 + t7 * t30 - t6 * t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t77 + 0.2e1 * t78, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, -t26, t25, t11, 0, 0, 0, 0, 0, -t16, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t7, t6, pkin(4) * t25 + t55, t7, -t6, -t7 * pkin(4) - t6 * qJ(5) + qJD(5) * t51, -t30 * t60 + t66, -t71 * t26 - 0.4e1 * t30 * t63, -t53, -t16, 0, t49 * t44 - t48 * t46, t48 * t44 + t49 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, 0, t26, -t25, -t11, 0, 0, 0, 0, 0, t16, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, qJ(5) * t83, -0.2e1 * t63, 0.2e1 * t60, 0, 0, 0, 0.2e1 * qJD(5) * t44 + 0.2e1 * t46 * t67, 0.2e1 * qJD(5) * t46 - 0.2e1 * t44 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, t7, 0, 0, 0, 0, 0, -t53, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t52, -t25, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, t44 * t68, t46 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
