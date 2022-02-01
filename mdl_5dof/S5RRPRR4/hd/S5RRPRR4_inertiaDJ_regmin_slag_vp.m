% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:32
% EndTime: 2022-01-20 10:48:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (345->70), mult. (851->111), div. (0->0), fcn. (692->8), ass. (0->70)
t85 = qJD(4) + qJD(5);
t55 = cos(qJ(4));
t84 = t55 * pkin(4);
t56 = cos(qJ(2));
t46 = t56 * pkin(1) + pkin(2);
t49 = sin(pkin(9));
t50 = cos(pkin(9));
t53 = sin(qJ(2));
t78 = t50 * t53;
t71 = pkin(1) * t78 + t49 * t46;
t25 = pkin(7) + t71;
t83 = -pkin(8) - t25;
t44 = t49 * pkin(2) + pkin(7);
t82 = -pkin(8) - t44;
t51 = sin(qJ(5));
t52 = sin(qJ(4));
t54 = cos(qJ(5));
t36 = t51 * t55 + t54 * t52;
t16 = t85 * t36;
t70 = pkin(1) * qJD(2);
t28 = (t49 * t56 + t78) * t70;
t68 = t52 * qJD(4);
t65 = pkin(4) * t68;
t19 = t28 + t65;
t79 = t49 * t53;
t59 = -pkin(1) * t79 + t50 * t46;
t24 = -pkin(3) - t59;
t22 = t24 - t84;
t77 = t51 * t52;
t35 = -t54 * t55 + t77;
t81 = t22 * t16 + t19 * t35;
t47 = t55 * qJD(4);
t69 = qJD(5) * t54;
t15 = -t54 * t47 - t55 * t69 + t85 * t77;
t80 = -t22 * t15 + t19 * t36;
t29 = (t50 * t56 - t79) * t70;
t76 = t52 * t29;
t75 = t55 * t29;
t74 = t24 * t47 + t28 * t52;
t45 = -t50 * pkin(2) - pkin(3);
t37 = t45 - t84;
t73 = t37 * t16 + t35 * t65;
t72 = -t37 * t15 + t36 * t65;
t67 = t53 * t70;
t66 = t56 * t70;
t64 = qJD(5) * t51 * pkin(4);
t63 = pkin(4) * t69;
t62 = t45 * t68;
t61 = t45 * t47;
t60 = t24 * t68 - t28 * t55;
t58 = qJD(4) * t83;
t57 = qJD(4) * t82;
t48 = t55 * pkin(8);
t39 = 0.2e1 * t52 * t47;
t34 = 0.2e1 * (-t52 ^ 2 + t55 ^ 2) * qJD(4);
t33 = t55 * t44 + t48;
t32 = t82 * t52;
t31 = t55 * t57;
t30 = t52 * t57;
t18 = t55 * t25 + t48;
t17 = t83 * t52;
t10 = -0.2e1 * t36 * t15;
t9 = t55 * t58 - t76;
t8 = t52 * t58 + t75;
t5 = -t51 * t30 + t54 * t31 + (-t32 * t51 - t33 * t54) * qJD(5);
t4 = -t54 * t30 - t51 * t31 + (-t32 * t54 + t33 * t51) * qJD(5);
t3 = 0.2e1 * t15 * t35 - 0.2e1 * t36 * t16;
t2 = -t51 * t8 + t54 * t9 + (-t17 * t51 - t18 * t54) * qJD(5);
t1 = -t51 * t9 - t54 * t8 + (-t17 * t54 + t18 * t51) * qJD(5);
t6 = [0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t66, -0.2e1 * t59 * t28 + 0.2e1 * t71 * t29, t39, t34, 0, 0, 0, 0.2e1 * t60, 0.2e1 * t74, t10, t3, 0, 0, 0, 0.2e1 * t81, 0.2e1 * t80; 0, 0, 0, 0, -t67, -t66, (-t28 * t50 + t29 * t49) * pkin(2), t39, t34, 0, 0, 0, t60 + t62, t61 + t74, t10, t3, 0, 0, 0, t73 + t81, t72 + t80; 0, 0, 0, 0, 0, 0, 0, t39, t34, 0, 0, 0, 0.2e1 * t62, 0.2e1 * t61, t10, t3, 0, 0, 0, 0.2e1 * t73, 0.2e1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t68, 0, -t25 * t47 - t76, t25 * t68 - t75, 0, 0, -t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t68, 0, -t44 * t47, t44 * t68, 0, 0, -t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t47, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
