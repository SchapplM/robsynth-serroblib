% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR13_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:52
% EndTime: 2019-12-31 18:32:55
% DurationCPUTime: 0.72s
% Computational Cost: add. (1116->104), mult. (2470->197), div. (0->0), fcn. (2347->6), ass. (0->77)
t41 = cos(pkin(8));
t82 = pkin(6) + qJ(2);
t30 = t82 * t41;
t40 = sin(pkin(8));
t87 = cos(qJ(3));
t69 = t87 * qJD(2);
t70 = qJD(3) * t87;
t43 = sin(qJ(3));
t79 = qJD(3) * t43;
t84 = t43 * t41;
t14 = (-t82 * t79 + t69) * t40 + qJD(2) * t84 + t30 * t70;
t23 = t40 * t79 - t41 * t70;
t94 = -t23 * pkin(4) + t14;
t42 = sin(qJ(5));
t38 = t42 ^ 2;
t44 = cos(qJ(5));
t39 = t44 ^ 2;
t68 = (t38 - t39) * qJD(5);
t71 = t82 * t40;
t27 = t87 * t71;
t13 = (qJD(2) * t40 + qJD(3) * t30) * t43 + qJD(3) * t27 - t41 * t69;
t93 = -0.2e1 * t23;
t29 = t87 * t40 + t84;
t24 = t29 * qJD(3);
t92 = 0.2e1 * t24;
t91 = 2 * qJD(4);
t90 = pkin(3) + pkin(7);
t28 = t43 * t40 - t87 * t41;
t8 = -t24 * pkin(4) - t13;
t88 = t28 * t8;
t86 = t38 * t24;
t22 = t39 * t24;
t85 = t42 * t24;
t83 = t44 * t24;
t78 = qJD(5) * t42;
t77 = qJD(5) * t44;
t76 = qJD(5) * t90;
t75 = qJ(4) * qJD(5);
t74 = t28 * t92;
t19 = t29 * t93;
t73 = t42 * t83;
t72 = t42 * t77;
t34 = -t41 * pkin(2) - pkin(1);
t20 = t43 * t30 + t27;
t26 = t28 ^ 2;
t67 = t26 * t72;
t66 = 0.2e1 * (t40 ^ 2 + t41 ^ 2) * qJD(2);
t60 = -t29 * qJ(4) + t34;
t11 = t90 * t28 + t60;
t59 = t29 * pkin(4) + t20;
t54 = t44 * t59;
t4 = -t42 * t11 + t54;
t5 = t44 * t11 + t42 * t59;
t65 = -t4 * t42 + t44 * t5;
t21 = t87 * t30 - t43 * t71;
t64 = -t21 * t13 + t20 * t14;
t63 = t23 * qJ(4) - t29 * qJD(4);
t62 = -qJ(4) * t24 - qJD(4) * t28;
t58 = -t42 * t23 + t29 * t77;
t57 = t44 * t23 + t29 * t78;
t56 = t28 * t77 + t85;
t55 = t28 * t78 - t83;
t53 = 0.2e1 * t23 * t28 - 0.2e1 * t29 * t24;
t52 = t8 + (qJ(4) * t28 + t29 * t90) * qJD(5);
t50 = t90 * t24 + t63;
t15 = -t28 * pkin(4) + t21;
t49 = -qJD(5) * t15 - t23 * t90 - t62;
t2 = -qJD(5) * t54 + t11 * t78 - t42 * t94 - t44 * t50;
t3 = -qJD(5) * t5 - t42 * t50 + t94 * t44;
t1 = t65 * qJD(5) - t2 * t42 + t3 * t44;
t47 = -t2 * t44 - t3 * t42 + (-t4 * t44 - t42 * t5) * qJD(5);
t46 = 0.2e1 * t13 * t28 + 0.2e1 * t14 * t29 - 0.2e1 * t20 * t23 - 0.2e1 * t21 * t24;
t35 = qJ(4) * t91;
t18 = t28 * pkin(3) + t60;
t10 = t24 * pkin(3) + t63;
t9 = -t28 * t68 + t73;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, qJ(2) * t66, t19, t53, 0, t74, 0, 0, t34 * t92, t34 * t93, t46, 0.2e1 * t64, 0, 0, 0, t19, t53, t74, t46, -0.2e1 * t10 * t28 - 0.2e1 * t18 * t24, -0.2e1 * t10 * t29 + 0.2e1 * t18 * t23, 0.2e1 * t18 * t10 + 0.2e1 * t64, 0.2e1 * t28 * t86 + 0.2e1 * t67, -0.2e1 * t26 * t68 + 0.4e1 * t28 * t73, 0.2e1 * t58 * t28 + 0.2e1 * t29 * t85, 0.2e1 * t28 * t22 - 0.2e1 * t67, -0.2e1 * t57 * t28 + 0.2e1 * t29 * t83, t19, 0.2e1 * t55 * t15 - 0.2e1 * t4 * t23 + 0.2e1 * t3 * t29 - 0.2e1 * t44 * t88, 0.2e1 * t56 * t15 + 0.2e1 * t2 * t29 + 0.2e1 * t5 * t23 + 0.2e1 * t42 * t88, 0.2e1 * t65 * t24 + 0.2e1 * t47 * t28, 0.2e1 * t15 * t8 - 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t23, t10, 0, 0, 0, 0, 0, 0, -t58, t57, t22 + t86, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, -t14, t13, 0, 0, 0, t23, t24, 0, 0, 0, pkin(3) * t23 + t62, t14, -t13, -t14 * pkin(3) - t13 * qJ(4) + t21 * qJD(4), t9, -0.4e1 * t28 * t72 + t22 - t86, -t57, -t9, -t58, 0, t52 * t42 - t49 * t44, t49 * t42 + t52 * t44, -t1, t8 * qJ(4) + t15 * qJD(4) - t1 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t35, -0.2e1 * t72, 0.2e1 * t68, 0, 0.2e1 * t72, 0, 0, 0.2e1 * qJD(4) * t42 + 0.2e1 * t44 * t75, 0.2e1 * qJD(4) * t44 - 0.2e1 * t42 * t75, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, t14, 0, 0, 0, 0, 0, 0, -t57, -t58, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t55, -t23, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, -t77, 0, t42 * t76, t44 * t76, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
