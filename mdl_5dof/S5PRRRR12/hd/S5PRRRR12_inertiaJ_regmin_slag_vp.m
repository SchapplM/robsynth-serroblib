% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:21
% DurationCPUTime: 0.27s
% Computational Cost: add. (269->56), mult. (790->92), div. (0->0), fcn. (885->12), ass. (0->58)
t46 = sin(qJ(4));
t64 = t46 * pkin(3);
t47 = sin(qJ(3));
t63 = t47 * pkin(2);
t51 = cos(qJ(3));
t38 = t51 * pkin(2);
t35 = t38 + pkin(3);
t50 = cos(qJ(4));
t20 = t50 * t35 - t46 * t63;
t19 = pkin(4) + t20;
t43 = cos(pkin(6));
t41 = sin(pkin(6));
t39 = t41 ^ 2;
t49 = cos(qJ(5));
t60 = t39 * t49;
t54 = t50 * t63;
t21 = -t46 * t35 - t54;
t36 = t41 * pkin(10);
t15 = -t21 + t36;
t45 = sin(qJ(5));
t58 = t43 * t49;
t7 = -t45 * t15 + t19 * t58;
t62 = t19 * t60 + t7 * t43;
t61 = t39 * t45;
t32 = t41 * t45;
t33 = t41 * t49;
t59 = t43 * t45;
t28 = t36 + t64;
t37 = t50 * pkin(3);
t34 = t37 + pkin(4);
t12 = -t45 * t28 + t34 * t58;
t57 = t12 * t43 + t34 * t60;
t22 = pkin(4) * t58 - pkin(10) * t32;
t56 = pkin(4) * t60 + t22 * t43;
t55 = 0.2e1 * t41 * t43;
t44 = cos(pkin(5));
t42 = sin(pkin(5));
t48 = sin(qJ(2));
t52 = cos(qJ(2));
t17 = (-t47 * t48 + t51 * t52) * t42;
t18 = (t47 * t52 + t48 * t51) * t42;
t9 = t50 * t17 - t46 * t18;
t53 = t41 * t44 + t43 * t9;
t40 = t43 ^ 2;
t31 = t39 * t45 ^ 2;
t27 = 0.2e1 * t45 * t60;
t26 = t49 * t55;
t25 = t45 * t55;
t23 = pkin(4) * t59 + pkin(10) * t33;
t13 = t49 * t28 + t34 * t59;
t10 = t46 * t17 + t50 * t18;
t8 = t49 * t15 + t19 * t59;
t5 = -t41 * t9 + t43 * t44;
t4 = t49 * t10 + t53 * t45;
t3 = -t45 * t10 + t53 * t49;
t2 = t5 * t32 - t4 * t43;
t1 = t3 * t43 - t5 * t33;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t42 * t52, -t42 * t48, 0, t17, -t18, 0, t9, -t10, 0, 0, 0, 0, 0, t1, t2; 0, 1, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t63, 1, 0.2e1 * t20, 0.2e1 * t21, t31, t27, t25, t26, t40, 0.2e1 * t62, -0.2e1 * t19 * t61 - 0.2e1 * t8 * t43; 0, 0, 0, 0, 0, t17, -t18, 0, t9, -t10, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 1, t38, -t63, 1, t20 + t37, -t54 + (-pkin(3) - t35) * t46, t31, t27, t25, t26, t40, t57 + t62, (-t13 - t8) * t43 + (-t19 - t34) * t61; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t64, t31, t27, t25, t26, t40, 0.2e1 * t57, -0.2e1 * t13 * t43 - 0.2e1 * t34 * t61; 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, t20, t21, t31, t27, t25, t26, t40, t56 + t62, (-t23 - t8) * t43 + (-pkin(4) - t19) * t61; 0, 0, 0, 0, 0, 0, 0, 1, t37, -t64, t31, t27, t25, t26, t40, t56 + t57, (-t13 - t23) * t43 + (-pkin(4) - t34) * t61; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t31, t27, t25, t26, t40, 0.2e1 * t56, -0.2e1 * pkin(4) * t61 - 0.2e1 * t23 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, t43, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, t43, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, t43, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
