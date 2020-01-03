% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:54
% DurationCPUTime: 0.78s
% Computational Cost: add. (1093->110), mult. (2310->209), div. (0->0), fcn. (2082->8), ass. (0->80)
t43 = sin(pkin(9));
t47 = cos(qJ(3));
t45 = sin(qJ(3));
t79 = cos(pkin(9));
t66 = t79 * t45;
t33 = t43 * t47 + t66;
t28 = t33 * qJD(3);
t65 = t79 * t47;
t76 = t45 * qJD(3);
t29 = qJD(3) * t65 - t43 * t76;
t87 = t43 * t45;
t32 = -t65 + t87;
t96 = -t33 * t28 - t29 * t32;
t44 = sin(qJ(5));
t41 = t44 ^ 2;
t46 = cos(qJ(5));
t42 = t46 ^ 2;
t64 = qJD(5) * (t41 - t42);
t38 = sin(pkin(8)) * pkin(1) + pkin(6);
t80 = qJ(4) + t38;
t30 = t80 * t47;
t16 = t79 * t30 - t80 * t87;
t63 = qJD(3) * t80;
t24 = t47 * qJD(4) - t45 * t63;
t53 = t45 * qJD(4) + t47 * t63;
t48 = t79 * t24 - t43 * t53;
t40 = -cos(pkin(8)) * pkin(1) - pkin(2);
t34 = -t47 * pkin(3) + t40;
t50 = -t32 * pkin(4) + t33 * pkin(7) - t34;
t49 = t46 * t50;
t71 = pkin(3) * t76;
t52 = t28 * pkin(4) - t29 * pkin(7) + t71;
t78 = qJD(5) * t44;
t1 = qJD(5) * t49 + t16 * t78 - t44 * t52 - t46 * t48;
t4 = t46 * t16 - t44 * t50;
t2 = -qJD(5) * t4 - t44 * t48 + t46 * t52;
t3 = -t44 * t16 - t49;
t59 = t3 * t44 - t4 * t46;
t95 = t59 * qJD(5) + t1 * t44 - t2 * t46;
t94 = 0.2e1 * qJD(3);
t15 = t43 * t30 + t80 * t66;
t7 = t43 * t24 + t79 * t53;
t93 = t15 * t7;
t92 = t7 * t44;
t91 = t32 * t28;
t90 = t33 * t29;
t89 = t33 * t46;
t88 = t41 * t29;
t25 = t42 * t29;
t86 = t44 * t28;
t85 = t46 * t28;
t84 = t46 * t29;
t83 = t32 * t84 + t33 * t85;
t77 = qJD(5) * t46;
t75 = t47 * qJD(3);
t74 = 0.2e1 * t91;
t73 = t40 * t94;
t39 = -t79 * pkin(3) - pkin(4);
t72 = 0.2e1 * qJD(5) * t39;
t70 = t33 * t78;
t69 = t44 * t77;
t68 = t45 * t75;
t67 = -0.4e1 * t44 * t89;
t31 = t33 ^ 2;
t62 = t31 * t69;
t60 = -t3 * t46 - t4 * t44;
t58 = t15 * t28 + t7 * t32;
t37 = t43 * pkin(3) + pkin(7);
t56 = -t28 * t37 + t29 * t39;
t55 = t32 * t37 - t33 * t39;
t13 = t32 * t77 + t86;
t54 = t44 * t29 + t33 * t77;
t12 = t70 - t84;
t51 = t60 * qJD(5) - t1 * t46 - t2 * t44;
t18 = t33 * t25;
t17 = t33 * t88;
t11 = t32 * t78 - t85;
t10 = -t25 - t88;
t5 = t33 * t64 - t44 * t84;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t68, (-t45 ^ 2 + t47 ^ 2) * t94, 0, -0.2e1 * t68, 0, 0, t45 * t73, t47 * t73, 0, 0, 0.2e1 * t90, 0.2e1 * t96, 0, t74, 0, 0, 0.2e1 * t34 * t28 + 0.2e1 * t32 * t71, 0.2e1 * t34 * t29 + 0.2e1 * t33 * t71, 0.2e1 * t15 * t29 - 0.2e1 * t16 * t28 - 0.2e1 * t32 * t48 + 0.2e1 * t7 * t33, 0.2e1 * t16 * t48 + 0.2e1 * t34 * t71 + 0.2e1 * t93, 0.2e1 * t18 - 0.2e1 * t62, t29 * t67 + 0.2e1 * t31 * t64, -0.2e1 * t32 * t70 + 0.2e1 * t83, 0.2e1 * t17 + 0.2e1 * t62, -0.2e1 * t54 * t32 - 0.2e1 * t33 * t86, t74, 0.2e1 * t54 * t15 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t28 + 0.2e1 * t33 * t92, 0.2e1 * t1 * t32 - 0.2e1 * t12 * t15 - 0.2e1 * t4 * t28 + 0.2e1 * t7 * t89, 0.2e1 * t60 * t29 + 0.2e1 * t95 * t33, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t29 + t33 * t48 + t58, 0, 0, 0, 0, 0, 0, 0, t46 * t96 + t83, 0, -t59 * t29 + t51 * t33 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t90 + 0.2e1 * t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 + 0.2e1 * t18 + 0.2e1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t76, 0, -t38 * t75, t38 * t76, 0, 0, 0, 0, t29, 0, -t28, 0, -t7, -t48, (-t28 * t43 - t79 * t29) * pkin(3), (t48 * t43 - t7 * t79) * pkin(3), -t5, qJD(5) * t67 + t25 - t88, t13, t5, -t11, 0, -t7 * t46 + t56 * t44 + (t15 * t44 - t55 * t46) * qJD(5), t92 + t56 * t46 + (t15 * t46 + t55 * t44) * qJD(5), t51, t37 * t51 + t7 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, (-t79 * t28 + t29 * t43) * pkin(3), 0, 0, 0, 0, 0, 0, t11, t13, -t10, t28 * t39 + (t41 + t42) * t37 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t64, 0, -0.2e1 * t69, 0, 0, t44 * t72, t46 * t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t29, 0, t71, 0, 0, 0, 0, 0, 0, -t11, -t13, t10, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, -t54, t28, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t78, 0, -t37 * t77, t37 * t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
