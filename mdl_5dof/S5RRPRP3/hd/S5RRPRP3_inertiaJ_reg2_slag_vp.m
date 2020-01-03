% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t45 = sin(pkin(8));
t43 = t45 ^ 2;
t46 = cos(pkin(8));
t44 = t46 ^ 2;
t57 = t43 + t44;
t58 = t57 * qJ(3);
t42 = t46 * pkin(7);
t31 = t46 * qJ(3) + t42;
t47 = sin(qJ(4));
t52 = (-pkin(7) - qJ(3)) * t45;
t65 = cos(qJ(4));
t18 = t47 * t31 - t65 * t52;
t20 = t65 * t31 + t47 * t52;
t61 = t47 * t45;
t27 = -t65 * t46 + t61;
t54 = t65 * t45;
t29 = t47 * t46 + t54;
t81 = t18 * t29 - t20 * t27;
t85 = 0.2e1 * t81;
t48 = sin(qJ(2));
t70 = t48 * pkin(1);
t37 = qJ(3) + t70;
t22 = t46 * t37 + t42;
t67 = -pkin(7) - t37;
t10 = t47 * t22 - t67 * t54;
t12 = t65 * t22 + t67 * t61;
t82 = t10 * t29 - t12 * t27;
t84 = 0.2e1 * t82;
t83 = t81 + t82;
t80 = t27 ^ 2;
t79 = 0.2e1 * t27;
t78 = -0.2e1 * t29;
t38 = -t46 * pkin(3) - pkin(2);
t49 = cos(qJ(2));
t69 = t49 * pkin(1);
t30 = t38 - t69;
t77 = 0.2e1 * t30;
t76 = 0.2e1 * t38;
t75 = 0.2e1 * t46;
t39 = -pkin(2) - t69;
t68 = pkin(2) - t39;
t13 = t27 * pkin(4) - t29 * qJ(5) + t38;
t8 = t13 - t69;
t66 = t13 + t8;
t62 = t29 * t27;
t60 = t30 + t38;
t59 = t57 * t37;
t56 = t10 ^ 2 + t12 ^ 2;
t55 = t18 ^ 2 + t20 ^ 2;
t53 = t10 * t18 + t12 * t20;
t34 = t45 * t75;
t25 = t29 ^ 2;
t17 = -0.2e1 * t62;
t16 = 0.2e1 * t62;
t14 = -pkin(4) * t29 - t27 * qJ(5);
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t70, 0, (t48 ^ 2 + t49 ^ 2) * pkin(1) ^ 2, t43, t34, 0, t44, 0, 0, -0.2e1 * t39 * t46, 0.2e1 * t39 * t45, 0.2e1 * t59, t57 * t37 ^ 2 + t39 ^ 2, t25, t17, 0, t80, 0, 0, t27 * t77, t29 * t77, t84, t30 ^ 2 + t56, t25, 0, t16, 0, 0, t80, t8 * t79, t84, t8 * t78, t8 ^ 2 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t70, 0, 0, t43, t34, 0, t44, 0, 0, t68 * t46, -t68 * t45, t58 + t59, -t39 * pkin(2) + t37 * t58, t25, t17, 0, t80, 0, 0, t60 * t27, t60 * t29, t83, t30 * t38 + t53, t25, 0, t16, 0, 0, t80, t66 * t27, t83, -t66 * t29, t8 * t13 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, t34, 0, t44, 0, 0, pkin(2) * t75, -0.2e1 * pkin(2) * t45, 0.2e1 * t58, t57 * qJ(3) ^ 2 + pkin(2) ^ 2, t25, t17, 0, t80, 0, 0, t27 * t76, t29 * t76, t85, t38 ^ 2 + t55, t25, 0, t16, 0, 0, t80, t13 * t79, t85, t13 * t78, t13 ^ 2 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, t39, 0, 0, 0, 0, 0, 0, t27, t29, 0, t30, 0, 0, 0, 0, 0, 0, t27, 0, -t29, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t27, t29, 0, t38, 0, 0, 0, 0, 0, 0, t27, 0, -t29, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, 0, -t10, -t12, 0, 0, 0, t29, 0, 0, t27, 0, -t10, t14, t12, -t10 * pkin(4) + t12 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, 0, -t18, -t20, 0, 0, 0, t29, 0, 0, t27, 0, -t18, t14, t20, -t18 * pkin(4) + t20 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
