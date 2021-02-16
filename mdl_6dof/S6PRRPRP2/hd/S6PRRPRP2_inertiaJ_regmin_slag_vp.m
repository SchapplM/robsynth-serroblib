% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:48
% EndTime: 2021-01-16 02:54:51
% DurationCPUTime: 0.67s
% Computational Cost: add. (684->104), mult. (1435->189), div. (0->0), fcn. (1739->10), ass. (0->67)
t41 = sin(qJ(3));
t38 = sin(pkin(6));
t57 = t38 * sin(qJ(2));
t61 = cos(pkin(6));
t70 = cos(qJ(3));
t23 = t61 * t41 + t70 * t57;
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t45 = -t41 * t57 + t61 * t70;
t9 = t37 * t23 - t39 * t45;
t78 = t9 ^ 2;
t34 = -t70 * pkin(3) - pkin(2);
t77 = 0.2e1 * t34;
t42 = cos(qJ(5));
t76 = -0.2e1 * t42;
t26 = t37 * t41 - t39 * t70;
t75 = t26 * pkin(5);
t74 = t37 * pkin(3);
t73 = t39 * pkin(3);
t40 = sin(qJ(5));
t72 = t9 * t40;
t71 = t9 * t42;
t32 = pkin(9) + t74;
t69 = t26 * t32;
t43 = cos(qJ(2));
t68 = t38 * t43;
t20 = t40 * t26;
t27 = t37 * t70 + t39 * t41;
t67 = t40 * t27;
t66 = t40 * t32;
t65 = t40 * t42;
t21 = t42 * t26;
t22 = t42 * t27;
t64 = t42 * t32;
t14 = t26 * pkin(4) - t27 * pkin(9) + t34;
t59 = t70 * pkin(8);
t29 = t70 * qJ(4) + t59;
t56 = (-pkin(8) - qJ(4)) * t41;
t18 = t39 * t29 + t37 * t56;
t4 = t40 * t14 + t42 * t18;
t35 = t40 ^ 2;
t36 = t42 ^ 2;
t63 = t35 + t36;
t62 = t26 * qJ(6);
t60 = 0.2e1 * t70;
t33 = -pkin(4) - t73;
t11 = t39 * t23 + t37 * t45;
t7 = t40 * t11 + t42 * t68;
t58 = -t7 * t26 + t9 * t67;
t55 = -t42 * t14 + t40 * t18;
t16 = t37 * t29 - t39 * t56;
t1 = t62 + t4;
t2 = t55 - t75;
t54 = t1 * t42 + t2 * t40;
t53 = t1 * t40 - t2 * t42;
t8 = t42 * t11 - t40 * t68;
t52 = t8 * t40 - t7 * t42;
t51 = t7 * t40 + t8 * t42;
t50 = t42 * pkin(5) + t40 * qJ(6);
t49 = pkin(5) * t40 - t42 * qJ(6);
t24 = t33 - t50;
t48 = t24 * t27 - t69;
t47 = t27 * t33 - t69;
t46 = t9 * t22 - t8 * t26;
t25 = t27 ^ 2;
t5 = t27 * t49 + t16;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 ^ 2 * t43 ^ 2 + t11 ^ 2 + t78, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t78; 0, 0, t68, -t57, 0, 0, 0, 0, 0, t70 * t68, -t41 * t68, -t26 * t68, -t27 * t68, -t11 * t26 + t9 * t27, t11 * t18 + t9 * t16 - t34 * t68, 0, 0, 0, 0, 0, t58, t46, t58, -t52 * t27, -t46, t8 * t1 + t7 * t2 + t9 * t5; 0, 1, 0, 0, t41 ^ 2, t41 * t60, 0, 0, 0, pkin(2) * t60, -0.2e1 * pkin(2) * t41, t26 * t77, t27 * t77, 0.2e1 * t16 * t27 - 0.2e1 * t18 * t26, t16 ^ 2 + t18 ^ 2 + t34 ^ 2, t36 * t25, -0.2e1 * t25 * t65, 0.2e1 * t26 * t22, -0.2e1 * t26 * t67, t26 ^ 2, 0.2e1 * t16 * t67 - 0.2e1 * t26 * t55, 0.2e1 * t16 * t22 - 0.2e1 * t4 * t26, -0.2e1 * t2 * t26 + 0.2e1 * t5 * t67, -0.2e1 * t53 * t27, 0.2e1 * t1 * t26 - 0.2e1 * t22 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t23, -t9, -t11, 0, (t11 * t37 - t39 * t9) * pkin(3), 0, 0, 0, 0, 0, -t71, t72, -t71, t51, -t72, t9 * t24 + t32 * t51; 0, 0, 0, 0, 0, 0, t41, t70, 0, -t41 * pkin(8), -t59, -t16, -t18, (-t26 * t37 - t27 * t39) * pkin(3), (-t16 * t39 + t18 * t37) * pkin(3), t40 * t22, (-t35 + t36) * t27, t20, t21, 0, -t16 * t42 + t40 * t47, t16 * t40 + t42 * t47, t40 * t48 - t5 * t42, t54, -t5 * t40 - t42 * t48, t5 * t24 + t32 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t73, -0.2e1 * t74, 0, (t37 ^ 2 + t39 ^ 2) * pkin(3) ^ 2, t35, 0.2e1 * t65, 0, 0, 0, t33 * t76, 0.2e1 * t33 * t40, t24 * t76, 0.2e1 * t63 * t32, -0.2e1 * t24 * t40, t32 ^ 2 * t63 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, t34, 0, 0, 0, 0, 0, t21, -t20, t21, -t63 * t27, t20, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, -t7 * pkin(5) + t8 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t67, t26, -t55, -t4, -t55 + 0.2e1 * t75, -t50 * t27, 0.2e1 * t62 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42, 0, -t66, -t64, -t66, -t49, t64, -t49 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, t42, 0, t40, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
