% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPR14
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:18
% DurationCPUTime: 0.40s
% Computational Cost: add. (471->83), mult. (1020->164), div. (0->0), fcn. (907->6), ass. (0->64)
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t67 = sin(pkin(8));
t53 = qJD(3) * t67;
t68 = cos(pkin(8));
t54 = qJD(3) * t68;
t15 = -t36 * t54 - t38 * t53;
t57 = t68 * t38;
t18 = -t67 * t36 + t57;
t74 = t18 * t15;
t14 = t36 * t53 - t38 * t54;
t56 = t67 * t38;
t19 = -t68 * t36 - t56;
t75 = t14 * t19;
t43 = 0.2e1 * t74 + 0.2e1 * t75;
t80 = (-t67 * t14 + t68 * t15) * pkin(3);
t17 = t18 ^ 2;
t79 = 2 * qJD(2);
t78 = 2 * qJD(5);
t63 = t38 * qJD(3);
t39 = -pkin(1) - pkin(6);
t69 = qJ(4) - t39;
t13 = -t36 * qJD(4) - t69 * t63;
t64 = t36 * qJD(3);
t42 = -t38 * qJD(4) + t69 * t64;
t3 = t67 * t13 - t68 * t42;
t35 = sin(qJ(5));
t77 = t3 * t35;
t37 = cos(qJ(5));
t76 = t3 * t37;
t73 = t35 * t15;
t72 = t37 * t15;
t34 = t37 ^ 2;
t71 = t35 ^ 2 - t34;
t70 = t36 * pkin(3) + qJ(2);
t66 = qJD(5) * t35;
t65 = qJD(5) * t37;
t22 = pkin(3) * t63 + qJD(2);
t62 = qJ(2) * qJD(3);
t30 = -t68 * pkin(3) - pkin(4);
t61 = t30 * t78;
t60 = t35 * t65;
t59 = t19 ^ 2 + t17;
t58 = -0.4e1 * t18 * t35 * t37;
t55 = t71 * qJD(5);
t21 = t69 * t36;
t10 = -t68 * t21 - t69 * t56;
t8 = -t19 * pkin(4) - t18 * pkin(7) + t70;
t52 = t37 * t10 + t35 * t8;
t51 = t35 * t10 - t37 * t8;
t29 = t67 * pkin(3) + pkin(7);
t49 = t14 * t29 + t15 * t30;
t48 = t18 * t30 + t19 * t29;
t47 = t35 * t14 + t19 * t65;
t46 = -t37 * t14 + t19 * t66;
t45 = -t18 * t65 - t73;
t44 = -t18 * t66 + t72;
t4 = t68 * t13 + t67 * t42;
t9 = -t67 * t21 + t69 * t57;
t40 = t10 * t14 + t9 * t15 + t3 * t18 + t4 * t19;
t5 = -t14 * pkin(4) - t15 * pkin(7) + t22;
t2 = -t52 * qJD(5) - t35 * t4 + t37 * t5;
t1 = t51 * qJD(5) - t35 * t5 - t37 * t4;
t6 = [0, 0, 0, 0, t79, qJ(2) * t79, -0.2e1 * t36 * t63, 0.2e1 * (t36 ^ 2 - t38 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t36 + 0.2e1 * t38 * t62, 0.2e1 * qJD(2) * t38 - 0.2e1 * t36 * t62, 0.2e1 * t40, 0.2e1 * t10 * t4 + 0.2e1 * t70 * t22 + 0.2e1 * t9 * t3, -0.2e1 * t17 * t60 + 0.2e1 * t34 * t74, t71 * t17 * t78 + t15 * t58, 0.2e1 * t46 * t18 - 0.2e1 * t19 * t72, 0.2e1 * t47 * t18 + 0.2e1 * t19 * t73, 0.2e1 * t75, -0.2e1 * t2 * t19 + 0.2e1 * t51 * t14 + 0.2e1 * t9 * t73 + 0.2e1 * (t9 * t65 + t77) * t18, -0.2e1 * t1 * t19 + 0.2e1 * t52 * t14 + 0.2e1 * t9 * t72 + 0.2e1 * (-t9 * t66 + t76) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t40, 0, 0, 0, 0, 0, -t35 * t43 - t59 * t65, -t37 * t43 + t59 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, -t39 * t64, -t39 * t63, -t80, (-t68 * t3 + t67 * t4) * pkin(3), -t18 * t55 + t35 * t72, qJD(5) * t58 - t71 * t15, -t47, t46, 0, -t76 + t49 * t35 + (t35 * t9 + t48 * t37) * qJD(5), t77 + t49 * t37 + (-t48 * t35 + t37 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, t80, 0, 0, 0, 0, 0, t44, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t60, -0.2e1 * t55, 0, 0, 0, t35 * t61, t37 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t45, -t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, -t29 * t65, t29 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
