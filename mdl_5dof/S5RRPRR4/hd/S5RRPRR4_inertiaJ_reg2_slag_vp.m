% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:33
% EndTime: 2022-01-20 10:48:36
% DurationCPUTime: 0.67s
% Computational Cost: add. (451->71), mult. (818->121), div. (0->0), fcn. (822->8), ass. (0->57)
t44 = sin(pkin(9));
t69 = t44 * pkin(2);
t36 = pkin(7) + t69;
t47 = sin(qJ(4));
t42 = t47 ^ 2;
t50 = cos(qJ(4));
t43 = t50 ^ 2;
t57 = t42 + t43;
t59 = t57 * t36;
t55 = -t50 * pkin(4) - pkin(3);
t51 = cos(qJ(2));
t41 = t51 * pkin(1);
t39 = t41 + pkin(2);
t45 = cos(pkin(9));
t48 = sin(qJ(2));
t66 = t48 * pkin(1);
t58 = -t45 * t39 + t44 * t66;
t14 = t55 + t58;
t75 = 0.2e1 * t14;
t68 = t45 * pkin(2);
t29 = t55 - t68;
t74 = 0.2e1 * t29;
t73 = 0.2e1 * t47;
t72 = -0.2e1 * t50;
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t26 = t46 * t47 - t49 * t50;
t28 = t46 * t50 + t49 * t47;
t56 = t45 * t66;
t21 = t44 * t39 + t56;
t19 = pkin(7) + t21;
t12 = (-pkin(8) - t19) * t47;
t40 = t50 * pkin(8);
t64 = t50 * t19;
t13 = t40 + t64;
t3 = t49 * t12 - t46 * t13;
t4 = t46 * t12 + t49 * t13;
t71 = -t4 * t26 - t3 * t28;
t22 = (-pkin(8) - t36) * t47;
t63 = t50 * t36;
t23 = t40 + t63;
t8 = t49 * t22 - t46 * t23;
t9 = t46 * t22 + t49 * t23;
t70 = -t9 * t26 - t8 * t28;
t67 = t46 * pkin(4);
t65 = t49 * pkin(4);
t62 = t14 + t29;
t61 = t57 * t19;
t18 = -pkin(3) + t58;
t37 = -pkin(3) - t68;
t60 = t18 + t37;
t35 = t50 * t73;
t25 = t28 ^ 2;
t24 = t26 ^ 2;
t11 = -0.2e1 * t28 * t26;
t10 = (-t26 * t46 - t28 * t49) * pkin(4);
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t66, 0, (t48 ^ 2 + t51 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t58, -0.2e1 * t21, 0, t21 ^ 2 + t58 ^ 2, t42, t35, 0, t43, 0, 0, t18 * t72, t18 * t73, 0.2e1 * t61, t57 * t19 ^ 2 + t18 ^ 2, t25, t11, 0, t24, 0, 0, t26 * t75, t28 * t75, 0.2e1 * t71, t14 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t66, 0, 0, 0, 0, 0, 0, 0, 1, -t58 + t68, -t56 + (-pkin(2) - t39) * t44, 0, (t21 * t44 - t45 * t58) * pkin(2), t42, t35, 0, t43, 0, 0, -t60 * t50, t60 * t47, t59 + t61, t18 * t37 + t19 * t59, t25, t11, 0, t24, 0, 0, t62 * t26, t62 * t28, t70 + t71, t14 * t29 + t3 * t8 + t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t69, 0, (t44 ^ 2 + t45 ^ 2) * pkin(2) ^ 2, t42, t35, 0, t43, 0, 0, t37 * t72, t37 * t73, 0.2e1 * t59, t57 * t36 ^ 2 + t37 ^ 2, t25, t11, 0, t24, 0, 0, t26 * t74, t28 * t74, 0.2e1 * t70, t29 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t26 + t4 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t26 + t9 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t50, 0, -t47 * t19, -t64, 0, 0, 0, 0, t28, 0, -t26, 0, t3, -t4, t10, (t3 * t49 + t4 * t46) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t50, 0, -t47 * t36, -t63, 0, 0, 0, 0, t28, 0, -t26, 0, t8, -t9, t10, (t46 * t9 + t49 * t8) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t28, 0, (-t26 * t49 + t28 * t46) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t67, 0, (t46 ^ 2 + t49 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t67, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
