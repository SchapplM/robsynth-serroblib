% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:44
% EndTime: 2021-01-15 23:52:47
% DurationCPUTime: 0.68s
% Computational Cost: add. (1689->111), mult. (3856->200), div. (0->0), fcn. (3524->6), ass. (0->75)
t89 = qJD(2) + qJD(3);
t88 = pkin(6) + pkin(7);
t53 = sin(qJ(2));
t87 = t53 * pkin(2);
t86 = cos(qJ(3));
t85 = cos(qJ(4));
t52 = sin(qJ(3));
t84 = qJD(3) * t52;
t51 = sin(qJ(4));
t83 = qJD(4) * t51;
t82 = t53 * qJD(2);
t54 = cos(qJ(2));
t81 = t54 * qJD(2);
t80 = t51 * t52 * pkin(2);
t79 = -0.2e1 * pkin(1) * qJD(2);
t78 = pkin(2) * t84;
t77 = pkin(3) * t83;
t76 = t85 * pkin(3);
t49 = -t54 * pkin(2) - pkin(1);
t75 = t85 * t52;
t74 = qJD(2) * t88;
t73 = t86 * qJD(3);
t72 = t85 * qJD(4);
t71 = pkin(2) * t73;
t70 = pkin(3) * t72;
t69 = t52 * t53 - t86 * t54;
t39 = t53 * t74;
t40 = t54 * t74;
t68 = t86 * t39 + t52 * t40;
t67 = t52 * t39 - t86 * t40;
t42 = t88 * t53;
t43 = t88 * t54;
t66 = t86 * t42 + t52 * t43;
t65 = t52 * t42 - t86 * t43;
t38 = t52 * t54 + t86 * t53;
t20 = -t38 * pkin(8) - t66;
t21 = -t69 * pkin(8) - t65;
t62 = t69 * qJD(3);
t58 = -t69 * qJD(2) - t62;
t55 = -t58 * pkin(8) + t42 * t84 - t43 * t73 + t67;
t57 = t89 * t38;
t56 = -t57 * pkin(8) - t42 * t73 - t43 * t84 - t68;
t3 = -t20 * t72 + t21 * t83 - t51 * t55 - t85 * t56;
t64 = t49 * t38;
t63 = t51 * t69;
t48 = t86 * pkin(2) + pkin(3);
t25 = -t48 * t72 - t85 * t71 + (qJD(3) + qJD(4)) * t80;
t61 = t85 * t69;
t33 = t69 * pkin(3) + t49;
t59 = (-t52 * t72 + (-t86 * t51 - t75) * qJD(3)) * pkin(2);
t24 = pkin(2) * t82 + t57 * pkin(3);
t4 = -t20 * t83 - t21 * t72 - t51 * t56 + t85 * t55;
t47 = t76 + pkin(4);
t46 = -0.2e1 * t70;
t45 = -0.2e1 * t77;
t36 = pkin(2) * t75 + t51 * t48;
t34 = t85 * t48 + pkin(4) - t80;
t32 = t85 * t38 - t63;
t31 = t51 * t38 + t61;
t26 = -t48 * t83 + t59;
t23 = 0.2e1 * t26;
t22 = 0.2e1 * t25;
t19 = -t70 + t25;
t18 = (-pkin(3) - t48) * t83 + t59;
t16 = t31 * pkin(4) + t33;
t15 = t65 * qJD(3) + t67;
t14 = t66 * qJD(3) + t68;
t10 = -qJD(4) * t63 + t38 * t72 + t51 * t58 + t85 * t57;
t9 = qJD(4) * t61 + t38 * t83 + t51 * t57 - t85 * t58;
t7 = -t31 * qJ(5) + t51 * t20 + t85 * t21;
t6 = -t32 * qJ(5) + t85 * t20 - t51 * t21;
t5 = t10 * pkin(4) + t24;
t2 = t9 * qJ(5) - t32 * qJD(5) + t4;
t1 = t10 * qJ(5) + t31 * qJD(5) + t3;
t8 = [0, 0, 0, 0.2e1 * t53 * t81, 0.2e1 * (-t53 ^ 2 + t54 ^ 2) * qJD(2), 0, 0, 0, t53 * t79, t54 * t79, 0.2e1 * t38 * t58, 0.2e1 * t89 * (-t38 ^ 2 + t69 ^ 2), 0, 0, 0, 0.2e1 * qJD(3) * t64 + 0.2e1 * (t69 * t87 + t64) * qJD(2), -0.2e1 * t49 * t62 + 0.2e1 * (t38 * t87 - t49 * t69) * qJD(2), -0.2e1 * t32 * t9, -0.2e1 * t32 * t10 + 0.2e1 * t9 * t31, 0, 0, 0, 0.2e1 * t33 * t10 + 0.2e1 * t24 * t31, 0.2e1 * t24 * t32 - 0.2e1 * t33 * t9, 0.2e1 * t16 * t10 + 0.2e1 * t5 * t31, -0.2e1 * t16 * t9 + 0.2e1 * t5 * t32, 0.2e1 * t1 * t31 - 0.2e1 * t7 * t10 - 0.2e1 * t2 * t32 + 0.2e1 * t6 * t9, -0.2e1 * t7 * t1 + 0.2e1 * t16 * t5 + 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t81, -t82, 0, -pkin(6) * t81, pkin(6) * t82, 0, 0, t58, -t57, 0, t15, t14, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, -t36 * t10 + t25 * t31 - t26 * t32 + t34 * t9, -t1 * t36 + t2 * t34 - t7 * t25 + t6 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t78, -0.2e1 * t71, 0, 0, 0, 0, 0, t23, t22, t23, t22, 0, -0.2e1 * t36 * t25 + 0.2e1 * t34 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t57, 0, t15, t14, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, t47 * t9 + (-t10 * t51 + (-t85 * t31 + t32 * t51) * qJD(4)) * pkin(3), t2 * t47 + (-t1 * t51 + (-t51 * t6 + t85 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t71, 0, 0, 0, 0, 0, t18, t19, t18, t19, 0, t26 * t47 + (-t25 * t51 + (-t34 * t51 + t85 * t36) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t46, t45, t46, 0, 0.2e1 * (t76 - t47) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, t4, t3, t2, t1, t9 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t25, t26, t25, 0, t26 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t70, -t77, -t70, 0, -pkin(4) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
