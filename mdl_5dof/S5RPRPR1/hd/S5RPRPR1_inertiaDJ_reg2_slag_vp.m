% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:47
% DurationCPUTime: 0.73s
% Computational Cost: add. (1074->88), mult. (2114->162), div. (0->0), fcn. (1980->6), ass. (0->63)
t50 = sin(pkin(8));
t51 = cos(pkin(8));
t54 = cos(qJ(3));
t76 = t54 * qJD(3);
t53 = sin(qJ(3));
t77 = t53 * qJD(3);
t33 = t50 * t77 - t51 * t76;
t38 = t50 * t54 + t51 * t53;
t34 = t38 * qJD(3);
t52 = sin(qJ(5));
t66 = t50 * t53 - t51 * t54;
t82 = cos(qJ(5));
t61 = t52 * t38 + t66 * t82;
t5 = qJD(5) * t61 + t82 * t33 + t52 * t34;
t62 = t82 * t38 - t52 * t66;
t88 = t5 * t62;
t55 = -pkin(1) - pkin(6);
t79 = qJ(4) - t55;
t41 = t79 * t54;
t85 = t79 * t53;
t91 = t50 * t41 + t51 * t85;
t14 = -t38 * pkin(7) - t91;
t12 = t91 * qJD(3) + t66 * qJD(4);
t56 = t34 * pkin(7) + t12;
t19 = -t51 * t41 + t50 * t85;
t13 = t19 * qJD(3) - t38 * qJD(4);
t57 = t33 * pkin(7) + t13;
t65 = pkin(7) * t66 + t19;
t59 = t82 * t65;
t78 = qJD(5) * t52;
t1 = -qJD(5) * t59 + t14 * t78 - t52 * t56 - t82 * t57;
t60 = t52 * t65;
t70 = qJD(5) * t82;
t2 = -qJD(5) * t60 - t14 * t70 - t52 * t57 + t82 * t56;
t4 = t82 * t14 + t60;
t97 = t1 * t62 + t2 * t61 + t4 * t5;
t72 = t51 * pkin(3) + pkin(4);
t83 = pkin(3) * t50;
t30 = -t52 * t83 + t82 * t72;
t23 = t30 * qJD(5);
t31 = t52 * t72 + t82 * t83;
t24 = t31 * qJD(5);
t96 = -t23 * t62 - t24 * t61 + t5 * t31;
t71 = t52 * t33 - t82 * t34;
t6 = t38 * t70 - t66 * t78 - t71;
t95 = t61 * t6;
t7 = -qJD(5) * t62 + t71;
t93 = -t61 * t7 - t88;
t87 = (t33 * t50 + t34 * t51) * pkin(3);
t84 = 2 * qJD(2);
t81 = t38 * t33;
t80 = t66 * t34;
t46 = t53 * pkin(3) + qJ(2);
t43 = pkin(3) * t76 + qJD(2);
t75 = qJ(2) * qJD(3);
t73 = t53 * t76;
t68 = -t80 + t81;
t58 = -t12 * t66 + t13 * t38 - t19 * t34 + t33 * t91;
t49 = qJ(2) * t84;
t22 = t38 * pkin(4) + t46;
t21 = -t33 * pkin(4) + t43;
t3 = -t52 * t14 + t59;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t49, -0.2e1 * t73, 0.2e1 * (t53 ^ 2 - t54 ^ 2) * qJD(3), 0, 0.2e1 * t73, 0, 0, 0.2e1 * qJD(2) * t53 + 0.2e1 * t54 * t75, 0.2e1 * qJD(2) * t54 - 0.2e1 * t53 * t75, 0, t49, 0.2e1 * t80, -0.2e1 * t33 * t66 + 0.2e1 * t34 * t38, 0, -0.2e1 * t81, 0, 0, -0.2e1 * t46 * t33 + 0.2e1 * t43 * t38, -0.2e1 * t46 * t34 - 0.2e1 * t43 * t66, -0.2e1 * t58, 0.2e1 * t19 * t12 - 0.2e1 * t13 * t91 + 0.2e1 * t46 * t43, 0.2e1 * t95, -0.2e1 * t5 * t61 + 0.2e1 * t6 * t62, 0, -0.2e1 * t88, 0, 0, 0.2e1 * t21 * t62 - 0.2e1 * t22 * t5, -0.2e1 * t21 * t61 - 0.2e1 * t22 * t6, 0.2e1 * t3 * t6 + 0.2e1 * t97, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t22 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t68, t58, 0, 0, 0, 0, 0, 0, 0, 0, t88 - t93 - t95, t3 * t7 - t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, -t76, 0, -t55 * t77, -t55 * t76, 0, 0, 0, 0, -t34, 0, t33, 0, t12, -t13, t87, (t12 * t51 + t13 * t50) * pkin(3), 0, 0, -t6, 0, t5, 0, t2, t1, t30 * t6 + t96, -t1 * t31 + t2 * t30 + t4 * t23 - t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, -t87, 0, 0, 0, 0, 0, 0, t7, t5, 0, t7 * t30 - t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24, -0.2e1 * t23, 0, 0.2e1 * t31 * t23 - 0.2e1 * t30 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, 0, t43, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, t5, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
