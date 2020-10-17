% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (259->48), mult. (453->83), div. (0->0), fcn. (378->6), ass. (0->44)
t30 = sin(qJ(4));
t28 = t30 ^ 2;
t33 = cos(qJ(4));
t39 = t33 ^ 2 + t28;
t26 = cos(qJ(2)) * pkin(1);
t22 = t26 + pkin(2);
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t52 = sin(qJ(2)) * pkin(1);
t38 = t34 * t52;
t9 = -t22 * t31 - t38;
t7 = pkin(8) - t9;
t57 = t39 * t7;
t53 = t31 * pkin(2);
t20 = pkin(8) + t53;
t60 = t39 * t20;
t59 = -0.2e1 * t30;
t58 = -0.2e1 * t33;
t56 = pkin(3) * t30;
t55 = t30 * pkin(8);
t54 = t30 * t7;
t51 = t33 * pkin(8);
t50 = t33 * t7;
t41 = -t34 * t22 + t31 * t52;
t6 = -pkin(3) + t41;
t49 = t6 * t33;
t11 = -pkin(4) * t33 - qJ(5) * t30 - pkin(3);
t1 = t11 + t41;
t25 = t34 * pkin(2);
t10 = t11 - t25;
t48 = -t1 - t10;
t47 = -t1 - t11;
t21 = -t25 - pkin(3);
t46 = t21 * t33;
t45 = t30 * t20;
t44 = t33 * t20;
t43 = -t10 - t11;
t40 = pkin(8) * t39;
t12 = -pkin(4) * t30 + qJ(5) * t33;
t27 = pkin(3) * t33;
t18 = 0.2e1 * t30 * t33;
t16 = t21 * t30;
t4 = t6 * t30;
t2 = [1, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t52, 1, -0.2e1 * t41, 0.2e1 * t9, t28, t18, 0, 0, 0, -0.2e1 * t49, 0.2e1 * t4, t1 * t58, 0.2e1 * t57, t1 * t59, t39 * t7 ^ 2 + t1 ^ 2; 0, 0, 0, 1, t26, -t52, 1, t25 - t41, -t38 + (-pkin(2) - t22) * t31, t28, t18, 0, 0, 0, (-t21 - t6) * t33, t16 + t4, t48 * t33, t60 + t57, t48 * t30, t1 * t10 + t60 * t7; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t53, t28, t18, 0, 0, 0, -0.2e1 * t46, 0.2e1 * t16, t10 * t58, 0.2e1 * t60, t10 * t59, t39 * t20 ^ 2 + t10 ^ 2; 0, 0, 0, 0, 0, 0, 1, -t41, t9, t28, t18, 0, 0, 0, t27 - t49, t4 - t56, t47 * t33, t40 + t57, t47 * t30, pkin(8) * t57 + t1 * t11; 0, 0, 0, 0, 0, 0, 1, t25, -t53, t28, t18, 0, 0, 0, t27 - t46, t16 - t56, t43 * t33, t40 + t60, t43 * t30, pkin(8) * t60 + t10 * t11; 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, t18, 0, 0, 0, 0.2e1 * t27, -0.2e1 * t56, t11 * t58, 0.2e1 * t40, t11 * t59, t39 * pkin(8) ^ 2 + t11 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t54, -t50, -t54, t12, t50, t12 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t45, -t44, -t45, t12, t44, t12 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t55, -t51, -t55, t12, t51, t12 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
