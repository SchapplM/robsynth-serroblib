% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:35
% DurationCPUTime: 0.62s
% Computational Cost: add. (398->116), mult. (1187->217), div. (0->0), fcn. (950->6), ass. (0->79)
t48 = sin(qJ(2));
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t73 = t49 * qJD(3);
t65 = t47 * t73;
t60 = 0.2e1 * t65;
t44 = t47 ^ 2;
t50 = cos(qJ(2));
t72 = t50 * qJD(2);
t66 = t44 * t72;
t95 = t48 * t60 + t66;
t45 = sin(pkin(8));
t42 = t45 ^ 2;
t46 = cos(pkin(8));
t94 = (t46 ^ 2 + t42) * qJD(4);
t75 = t47 * qJD(3);
t64 = t48 * t75;
t93 = t49 * t72 - t64;
t28 = -t46 * pkin(4) - t45 * qJ(5) - pkin(3);
t81 = t47 * qJ(4);
t92 = qJD(3) * (-t28 * t49 + t81);
t41 = pkin(6) * t73;
t56 = pkin(4) * t45 - qJ(5) * t46;
t79 = qJD(5) * t47;
t7 = -t46 * t79 + t56 * t73 + t41;
t91 = t7 * t45;
t90 = t7 * t46;
t84 = t48 * t49;
t68 = t46 * t84;
t23 = -t50 * t45 + t68;
t83 = t50 * t46;
t11 = -t46 * t64 + (t45 * t48 + t49 * t83) * qJD(2);
t6 = t11 * t46;
t89 = t23 * t46 * qJD(4) + qJ(4) * t6;
t88 = t45 * t49;
t21 = -t47 * qJD(4) + (pkin(3) * t47 - qJ(4) * t49) * qJD(3);
t87 = t46 * t21;
t57 = -t49 * pkin(3) - t81;
t29 = -pkin(2) + t57;
t86 = t46 * t29;
t85 = t46 * t49;
t38 = pkin(6) * t85;
t15 = t45 * t29 + t38;
t82 = qJ(4) * t94;
t80 = qJD(4) * t49;
t78 = t42 * qJD(5);
t77 = t45 * qJD(4);
t76 = t45 * qJD(5);
t74 = t48 * qJD(2);
t70 = pkin(6) * t88;
t69 = -0.2e1 * pkin(2) * qJD(3);
t67 = pkin(6) * t75;
t62 = pkin(6) * t45 + pkin(4);
t10 = t93 * t45 - t46 * t74;
t61 = t10 * t45 + t6;
t8 = t45 * t67 + t87;
t19 = t45 * t21;
t9 = -t46 * t67 + t19;
t58 = -t8 * t45 + t9 * t46;
t22 = t45 * t84 + t83;
t54 = qJ(4) * t10 + qJD(4) * t22;
t53 = 0.2e1 * t22 * t10 + 0.2e1 * t23 * t11 + 0.2e1 * (t65 * t48 + t66) * t48;
t24 = t47 * t72 + t48 * t73;
t52 = t10 * t49 - t22 * t75 + t95 * t45;
t51 = (t10 * t46 - t11 * t45) * t47 + (t22 * t46 - t23 * t45) * t73;
t36 = t46 * t73;
t35 = t45 * t73;
t34 = t49 * t77;
t27 = 0.2e1 * t94;
t20 = (pkin(6) + t56) * t47;
t17 = t24 * t46;
t16 = t24 * t45;
t14 = -t70 + t86;
t13 = t62 * t49 - t86;
t12 = -t49 * qJ(5) + t15;
t4 = -t62 * t75 - t87;
t3 = -t49 * qJD(5) + t19 + (-pkin(6) * t46 + qJ(5)) * t75;
t1 = t46 * t66 + t11 * t49 + (-t23 + 0.2e1 * t68) * t75;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, t53; 0, 0, -t74, -t72, 0, 0, 0, 0, 0, -t49 * t74 - t50 * t75, t47 * t74 - t50 * t73, t52, t1, t51, t95 * pkin(6) - t10 * t14 + t11 * t15 - t22 * t8 + t23 * t9, t52, t51, -t1, t47 * t48 * t7 + t10 * t13 + t11 * t12 + t20 * t24 + t22 * t4 + t23 * t3; 0, 0, 0, 0, t60, 0.2e1 * (t49 ^ 2 - t44) * qJD(3), 0, 0, 0, t47 * t69, t49 * t69, -0.2e1 * t8 * t49 + 0.2e1 * (t14 + 0.2e1 * t70) * t75, 0.2e1 * t9 * t49 + 0.2e1 * (-t15 + 0.2e1 * t38) * t75, 0.2e1 * (-t45 * t9 - t46 * t8) * t47 + 0.2e1 * (-t14 * t46 - t15 * t45) * t73, 0.2e1 * pkin(6) ^ 2 * t65 + 0.2e1 * t14 * t8 + 0.2e1 * t15 * t9, 0.2e1 * t47 * t91 + 0.2e1 * t4 * t49 + 0.2e1 * (-t13 * t47 + t20 * t88) * qJD(3), 0.2e1 * (-t3 * t45 + t4 * t46) * t47 + 0.2e1 * (-t12 * t45 + t13 * t46) * t73, -0.2e1 * t47 * t90 - 0.2e1 * t3 * t49 + 0.2e1 * (t12 * t47 - t20 * t85) * qJD(3), 0.2e1 * t12 * t3 + 0.2e1 * t13 * t4 + 0.2e1 * t20 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t93, -t17, t16, t61, -pkin(3) * t24 + t45 * t54 + t89, -t17, t61, -t16, t24 * t28 + (-t48 * t79 + t54) * t45 + t89; 0, 0, 0, 0, 0, 0, t73, -t75, 0, -t41, t67, t34 + (t45 * t57 - t38) * qJD(3), t46 * t80 + (t57 * t46 + t70) * qJD(3), t58, -pkin(3) * t41 + (-t14 * t45 + t15 * t46) * qJD(4) + t58 * qJ(4), -t45 * t92 - t47 * t78 + t34 - t90, t3 * t46 + t4 * t45, -t91 + (t47 * t76 - t80 + t92) * t46, t7 * t28 + (qJ(4) * t3 + qJD(4) * t12) * t46 + (qJ(4) * t4 + qJD(4) * t13 - qJD(5) * t20) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0.2e1 * t82, 0.2e1 * t46 * t76, t27, 0.2e1 * t78, -0.2e1 * t28 * t76 + 0.2e1 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, t41, t35, 0, -t36, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t36, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
