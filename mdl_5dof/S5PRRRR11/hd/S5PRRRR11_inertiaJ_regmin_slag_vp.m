% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:57
% DurationCPUTime: 0.42s
% Computational Cost: add. (307->61), mult. (876->100), div. (0->0), fcn. (984->8), ass. (0->49)
t32 = cos(pkin(5));
t53 = -0.2e1 * t32;
t52 = 0.2e1 * t32;
t51 = pkin(7) + pkin(8);
t35 = sin(qJ(3));
t50 = pkin(2) * t35;
t49 = t32 * pkin(3);
t48 = t32 * pkin(4);
t33 = sin(qJ(5));
t47 = t33 * pkin(4);
t34 = sin(qJ(4));
t46 = t34 * pkin(3);
t36 = cos(qJ(5));
t27 = t36 * pkin(4);
t31 = sin(pkin(5));
t38 = cos(qJ(3));
t25 = t31 * t38;
t37 = cos(qJ(4));
t43 = t31 * t35;
t15 = -t37 * t25 + t34 * t43;
t24 = t32 * t38 * pkin(2);
t12 = -t51 * t43 + t24 + t49;
t40 = t32 * t50;
t13 = t51 * t25 + t40;
t42 = t37 * t13;
t7 = t34 * t12 + t42;
t5 = -t15 * pkin(9) + t7;
t45 = t36 * t5;
t28 = t37 * pkin(3);
t29 = t31 ^ 2;
t44 = t29 * t38;
t41 = t31 * t52;
t39 = t36 * t46;
t16 = (t34 * t38 + t35 * t37) * t31;
t6 = t37 * t12 - t34 * t13;
t4 = -t16 * pkin(9) + t48 + t6;
t1 = -t33 * t5 + t36 * t4;
t26 = t28 + pkin(4);
t17 = t36 * t26 - t33 * t46;
t21 = (-pkin(3) * t38 - pkin(2)) * t31;
t2 = t33 * t4 + t45;
t30 = t32 ^ 2;
t20 = pkin(7) * t25 + t40;
t19 = -pkin(7) * t43 + t24;
t18 = t33 * t26 + t39;
t10 = t15 * pkin(4) + t21;
t9 = -t33 * t15 + t36 * t16;
t8 = t36 * t15 + t33 * t16;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, t29 * t35 ^ 2, 0.2e1 * t35 * t44, t35 * t41, t38 * t41, t30, 0.2e1 * pkin(2) * t44 + 0.2e1 * t19 * t32, -0.2e1 * t20 * t32 - 0.2e1 * t29 * t50, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t52, t15 * t53, t30, 0.2e1 * t21 * t15 + 0.2e1 * t6 * t32, 0.2e1 * t21 * t16 - 0.2e1 * t7 * t32, t9 ^ 2, -0.2e1 * t9 * t8, t9 * t52, t8 * t53, t30, 0.2e1 * t1 * t32 + 0.2e1 * t10 * t8, 0.2e1 * t10 * t9 - 0.2e1 * t2 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t43, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, t43, t25, t32, t19, -t20, 0, 0, t16, -t15, t32, t32 * t28 + t6, -t42 + (-t12 - t49) * t34, 0, 0, t9, -t8, t32, t17 * t32 + t1, -t18 * t32 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, -0.2e1 * t46, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t32, t6, -t7, 0, 0, t9, -t8, t32, t32 * t27 + t1, -t45 + (-t4 - t48) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t28, -t46, 0, 0, 0, 0, 1, t17 + t27, -t39 + (-pkin(4) - t26) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
