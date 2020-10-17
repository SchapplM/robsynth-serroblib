% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:11
% DurationCPUTime: 0.73s
% Computational Cost: add. (518->79), mult. (908->139), div. (0->0), fcn. (989->6), ass. (0->60)
t42 = sin(qJ(4));
t70 = t42 * pkin(3);
t31 = pkin(8) + t70;
t41 = sin(qJ(5));
t36 = t41 ^ 2;
t44 = cos(qJ(5));
t38 = t44 ^ 2;
t59 = t36 + t38;
t77 = t59 * t31;
t43 = sin(qJ(3));
t45 = cos(qJ(4));
t46 = cos(qJ(3));
t18 = t42 * t46 + t45 * t43;
t20 = -t42 * t43 + t45 * t46;
t76 = (t18 * t42 + t20 * t45) * pkin(3);
t16 = t18 ^ 2;
t17 = t20 ^ 2;
t75 = t16 + t17;
t47 = -pkin(1) - pkin(6);
t65 = -pkin(7) + t47;
t23 = t65 * t43;
t56 = t65 * t46;
t8 = t42 * t23 - t45 * t56;
t74 = t8 ^ 2;
t29 = t43 * pkin(3) + qJ(2);
t73 = 0.2e1 * t29;
t72 = 0.2e1 * qJ(2);
t71 = t20 * pkin(4);
t69 = t45 * pkin(3);
t68 = t8 * t20;
t67 = t8 * t44;
t32 = -pkin(4) - t69;
t66 = pkin(4) - t32;
t64 = t20 * t32;
t14 = t41 * t20;
t63 = t41 * t44;
t62 = t44 * t20;
t60 = pkin(8) * t59;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t26 = t37 + t39;
t58 = -0.2e1 * t20 * t18;
t6 = t59 * t18;
t54 = -pkin(8) * t18 - t71;
t10 = t45 * t23 + t42 * t56;
t4 = t18 * pkin(4) - t20 * pkin(8) + t29;
t2 = -t41 * t10 + t44 * t4;
t3 = t44 * t10 + t41 * t4;
t1 = -t2 * t41 + t3 * t44;
t53 = t10 * t18 - t68;
t52 = -t18 * t31 + t64;
t48 = qJ(2) ^ 2;
t27 = 0.2e1 * t63;
t22 = t26 * t47;
t13 = t44 * t18;
t12 = t41 * t18;
t11 = t41 * t62;
t7 = t8 * t41;
t5 = (-t36 + t38) * t20;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t72, (pkin(1) ^ 2) + t48, t39, -0.2e1 * t46 * t43, 0, t37, 0, 0, t43 * t72, t46 * t72, -0.2e1 * t22, t26 * t47 ^ 2 + t48, t17, t58, 0, t16, 0, 0, t18 * t73, t20 * t73, -0.2e1 * t53, t10 ^ 2 + t29 ^ 2 + t74, t38 * t17, -0.2e1 * t17 * t63, 0.2e1 * t18 * t62, t36 * t17, t41 * t58, t16, 0.2e1 * t8 * t14 + 0.2e1 * t2 * t18, -0.2e1 * t3 * t18 + 0.2e1 * t8 * t62, 0.2e1 * (-t2 * t44 - t3 * t41) * t20, t2 ^ 2 + t3 ^ 2 + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t26, t22, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t53, 0, 0, 0, 0, 0, 0, -t75 * t41, -t75 * t44, 0, t1 * t18 - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t16 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t43, 0, t46 * t47, -t43 * t47, 0, 0, 0, 0, t20, 0, -t18, 0, -t8, -t10, -t76, (t10 * t42 - t45 * t8) * pkin(3), t11, t5, t12, -t11, t13, 0, t52 * t41 - t67, t52 * t44 + t7, t1, t1 * t31 + t8 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, t76, 0, 0, 0, 0, 0, 0, t62, -t14, t6, t18 * t77 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t70, 0, (t42 ^ 2 + t45 ^ 2) * pkin(3) ^ 2, t36, t27, 0, t38, 0, 0, -0.2e1 * t32 * t44, 0.2e1 * t32 * t41, 0.2e1 * t77, t59 * t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, -t8, -t10, 0, 0, t11, t5, t12, -t11, t13, 0, t54 * t41 - t67, t54 * t44 + t7, t1, -t8 * pkin(4) + t1 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t14, t6, pkin(8) * t6 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t70, 0, 0, t36, t27, 0, t38, 0, 0, t66 * t44, -t66 * t41, t60 + t77, -t32 * pkin(4) + pkin(8) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t27, 0, t38, 0, 0, 0.2e1 * pkin(4) * t44, -0.2e1 * pkin(4) * t41, 0.2e1 * t60, t59 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, -t14, t18, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * t31, -t44 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * pkin(8), -t44 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
