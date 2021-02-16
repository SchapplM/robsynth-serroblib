% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:57
% EndTime: 2021-01-15 20:30:03
% DurationCPUTime: 0.82s
% Computational Cost: add. (1199->143), mult. (2742->276), div. (0->0), fcn. (2422->6), ass. (0->82)
t52 = sin(qJ(2));
t54 = cos(qJ(2));
t87 = -qJ(3) - pkin(6);
t69 = qJD(2) * t87;
t28 = t54 * qJD(3) + t52 * t69;
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t59 = -t52 * qJD(3) + t54 * t69;
t13 = t50 * t28 + t49 * t59;
t34 = t49 * t52 - t50 * t54;
t35 = t49 * t54 + t50 * t52;
t44 = -t54 * pkin(2) - pkin(1);
t60 = -t34 * pkin(3) + t35 * pkin(7) - t44;
t98 = qJD(4) * t60 - t13;
t38 = t87 * t52;
t39 = t87 * t54;
t23 = t49 * t38 - t50 * t39;
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t86 = t53 * t23 - t51 * t60;
t29 = t35 * qJD(2);
t80 = t54 * qJD(2);
t81 = t52 * qJD(2);
t30 = -t49 * t81 + t50 * t80;
t76 = pkin(2) * t81;
t57 = t29 * pkin(3) - t30 * pkin(7) + t76;
t10 = t53 * t57;
t58 = t30 * qJ(5) + qJD(4) * t23 + t35 * qJD(5);
t84 = qJ(5) * t35;
t73 = qJD(4) * t84;
t55 = -t58 * t53 + t10 + (t73 + t98) * t51;
t95 = t29 * pkin(4);
t1 = t55 + t95;
t79 = -t51 * t57 + t98 * t53;
t2 = t58 * t51 + t53 * t73 + t79;
t70 = -t51 * t23 - t53 * t60;
t5 = t34 * pkin(4) - t53 * t84 + t70;
t6 = -t51 * t84 + t86;
t97 = -t1 * t53 + t2 * t51 + (t5 * t51 - t53 * t6) * qJD(4);
t96 = 0.2e1 * qJD(4);
t94 = t53 * pkin(4);
t93 = t35 * t51;
t92 = t35 * t53;
t91 = t51 * t29;
t90 = t51 * t30;
t89 = t53 * t29;
t88 = t53 * t30;
t47 = t51 ^ 2;
t48 = t53 ^ 2;
t85 = t47 - t48;
t42 = t49 * pkin(2) + pkin(7);
t83 = qJ(5) + t42;
t82 = qJD(4) * t51;
t46 = qJD(4) * t53;
t78 = -0.2e1 * pkin(1) * qJD(2);
t43 = -t50 * pkin(2) - pkin(3);
t77 = t43 * t96;
t75 = pkin(4) * t82;
t74 = t51 * t46;
t37 = t43 - t94;
t72 = -t37 + t94;
t71 = -0.4e1 * t51 * t92;
t12 = t49 * t28 - t50 * t59;
t22 = -t50 * t38 - t49 * t39;
t68 = t85 * qJD(4);
t65 = pkin(4) * t47 + t37 * t53;
t64 = -t29 * t42 + t30 * t43;
t63 = t34 * t42 - t35 * t43;
t62 = t34 * t46 + t91;
t20 = t35 * t46 + t90;
t61 = t35 * t82 - t88;
t33 = t35 ^ 2;
t32 = t83 * t53;
t31 = t83 * t51;
t25 = -t51 * qJD(5) - t83 * t46;
t24 = -t53 * qJD(5) + t83 * t82;
t17 = -t34 * t82 + t89;
t11 = pkin(4) * t93 + t22;
t7 = t20 * pkin(4) + t12;
t4 = -t86 * qJD(4) - t51 * t13 + t10;
t3 = t23 * t82 + t79;
t8 = [0, 0, 0, 0.2e1 * t52 * t80, 0.2e1 * (-t52 ^ 2 + t54 ^ 2) * qJD(2), 0, 0, 0, t52 * t78, t54 * t78, 0.2e1 * t44 * t29 + 0.2e1 * t34 * t76, 0.2e1 * t44 * t30 + 0.2e1 * t35 * t76, 0.2e1 * t12 * t35 - 0.2e1 * t13 * t34 + 0.2e1 * t22 * t30 - 0.2e1 * t23 * t29, 0.2e1 * t22 * t12 + 0.2e1 * t23 * t13 + 0.2e1 * t44 * t76, 0.2e1 * t48 * t35 * t30 - 0.2e1 * t33 * t74, t85 * t33 * t96 + t30 * t71, -0.2e1 * t61 * t34 + 0.2e1 * t35 * t89, -0.2e1 * t20 * t34 - 0.2e1 * t35 * t91, 0.2e1 * t34 * t29, 0.2e1 * t12 * t93 + 0.2e1 * t20 * t22 + 0.2e1 * t70 * t29 + 0.2e1 * t4 * t34, 0.2e1 * t12 * t92 - 0.2e1 * t61 * t22 - 0.2e1 * t86 * t29 + 0.2e1 * t3 * t34, 0.2e1 * t1 * t34 + 0.2e1 * t20 * t11 + 0.2e1 * t5 * t29 + 0.2e1 * t7 * t93, -0.2e1 * t61 * t11 + 0.2e1 * t2 * t34 - 0.2e1 * t6 * t29 + 0.2e1 * t7 * t92, 0.2e1 * (-t5 * t53 - t51 * t6) * t30 + 0.2e1 * t97 * t35, 0.2e1 * t5 * t1 + 0.2e1 * t11 * t7 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, t80, -t81, 0, -pkin(6) * t80, pkin(6) * t81, -t12, -t13, (-t29 * t49 - t30 * t50) * pkin(2), (-t12 * t50 + t13 * t49) * pkin(2), -t35 * t68 + t51 * t88, qJD(4) * t71 - t85 * t30, t62, t17, 0, -t12 * t53 + t64 * t51 + (t22 * t51 - t63 * t53) * qJD(4), t12 * t51 + t64 * t53 + (t22 * t53 + t63 * t51) * qJD(4), t37 * t90 + t25 * t34 - t31 * t29 - t7 * t53 + (t11 * t51 + t65 * t35) * qJD(4), t37 * t88 + t24 * t34 - t32 * t29 + t7 * t51 + (t11 * t53 + t72 * t93) * qJD(4), (-t25 * t35 + t30 * t31 - t2 + (-t32 * t35 - t5) * qJD(4)) * t53 + (t24 * t35 - t30 * t32 - t1 + (-t31 * t35 - t6) * qJD(4)) * t51, -t1 * t31 + t11 * t75 - t2 * t32 - t6 * t24 + t5 * t25 + t7 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, -0.2e1 * t68, 0, 0, 0, t51 * t77, t53 * t77, -0.2e1 * t72 * t82, t65 * t96, -0.2e1 * t24 * t53 - 0.2e1 * t25 * t51 + 0.2e1 * (t31 * t53 - t32 * t51) * qJD(4), -0.2e1 * t32 * t24 - 0.2e1 * t31 * t25 + 0.2e1 * t37 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, t76, 0, 0, 0, 0, 0, t17, -t62, t17, -t62, (-t47 - t48) * t30, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t51 + t25 * t53 + (t31 * t51 + t32 * t53) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t20, t29, t4, t3, t55 + 0.2e1 * t95, t2, t61 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t82, 0, -t42 * t46, t42 * t82, t25, t24, -pkin(4) * t46, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t46, -t82, -t46, 0, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t61, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t46, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
