% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:23
% DurationCPUTime: 0.72s
% Computational Cost: add. (375->80), mult. (1078->124), div. (0->0), fcn. (848->6), ass. (0->61)
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t78 = cos(qJ(3));
t79 = cos(qJ(2));
t29 = t51 * t52 - t78 * t79;
t84 = qJD(2) + qJD(3);
t11 = t29 * t84;
t50 = sin(qJ(4));
t48 = t50 ^ 2;
t53 = cos(qJ(4));
t49 = t53 ^ 2;
t86 = t48 + t49;
t2 = t86 * t11;
t62 = t78 * pkin(2);
t59 = qJD(3) * t62;
t89 = t86 * t59;
t69 = pkin(2) * qJD(3);
t87 = t69 * t78;
t85 = -pkin(4) * t50 + qJ(5) * t53;
t83 = 0.2e1 * (-t48 + t49) * qJD(4);
t82 = 2 * qJD(5);
t80 = t2 * pkin(7);
t30 = t51 * t79 + t78 * t52;
t12 = t84 * t30;
t77 = t29 * t12;
t46 = t50 * qJD(4);
t47 = t53 * qJD(4);
t22 = pkin(4) * t46 - qJ(5) * t47 - t50 * qJD(5);
t65 = t51 * t69;
t13 = t22 + t65;
t73 = -t13 - t22;
t43 = t51 * pkin(2) + pkin(7);
t72 = t89 * t43;
t71 = t89 * pkin(7);
t44 = -t62 - pkin(3);
t70 = t44 * t47 + t50 * t65;
t67 = pkin(3) * t46;
t66 = pkin(3) * t47;
t64 = pkin(7) * t46;
t63 = pkin(7) * t47;
t61 = t50 * t47;
t60 = -t2 * t43 + t89 * t30;
t57 = -t53 * pkin(4) - t50 * qJ(5);
t54 = t44 * t46 - t53 * t65;
t33 = -pkin(3) + t57;
t21 = t57 * qJD(4) + t53 * qJD(5);
t39 = -0.2e1 * t61;
t38 = 0.2e1 * t61;
t24 = -t62 + t33;
t23 = t33 * t46;
t20 = t86 * t87;
t19 = t24 * t46;
t18 = t43 * t47 + t50 * t59;
t17 = t43 * t46 - t53 * t59;
t16 = 0.2e1 * t20;
t6 = -t12 * t53 + t29 * t46;
t5 = t12 * t50 + t29 * t47;
t4 = -t50 * t11 + t30 * t47;
t3 = t53 * t11 + t30 * t46;
t1 = -0.2e1 * t30 * t2 + 0.2e1 * t77;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t30 * t11 + 0.2e1 * t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * qJD(2), -t79 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, (-t78 * t12 - t11 * t51 + (t29 * t51 + t78 * t30) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t6, t5, -t2, t12 * t44 + t29 * t65 + t60, 0, 0, 0, 0, 0, 0, t6, -t2, -t5, t12 * t24 + t29 * t13 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t65, -0.2e1 * t59, 0, 0, t38, t83, 0, t39, 0, 0, 0.2e1 * t54, 0.2e1 * t70, t16, 0.2e1 * t44 * t65 + 0.2e1 * t72, t38, 0, -t83, 0, 0, t39, -0.2e1 * t13 * t53 + 0.2e1 * t19, t16, -0.2e1 * t13 * t50 - 0.2e1 * t24 * t47, 0.2e1 * t24 * t13 + 0.2e1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, -t2, -t12 * pkin(3) - t80, 0, 0, 0, 0, 0, 0, t6, -t2, -t5, t12 * t33 + t29 * t22 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t59, 0, 0, t38, t83, 0, t39, 0, 0, t54 - t67, -t66 + t70, t20, -pkin(3) * t65 + t71, t38, 0, -t83, 0, 0, t39, t73 * t53 + t19 + t23, t20, t73 * t50 + (-t24 - t33) * t47, t13 * t33 + t24 * t22 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t83, 0, t39, 0, 0, -0.2e1 * t67, -0.2e1 * t66, 0, 0, t38, 0, -t83, 0, 0, t39, -0.2e1 * t22 * t53 + 0.2e1 * t23, 0, -0.2e1 * t22 * t50 - 0.2e1 * t33 * t47, 0.2e1 * t33 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, -t3, -t85 * t11 + t21 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t46, 0, -t18, t17, 0, 0, 0, t47, 0, 0, t46, 0, -t18, t21, -t17, t21 * t43 + t85 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t46, 0, -t63, t64, 0, 0, 0, t47, 0, 0, t46, 0, -t63, t21, -t64, t21 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, qJ(5) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
