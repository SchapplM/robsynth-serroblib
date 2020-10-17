% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:47:59
% EndTime: 2019-12-31 19:48:00
% DurationCPUTime: 0.33s
% Computational Cost: add. (245->63), mult. (478->120), div. (0->0), fcn. (491->6), ass. (0->51)
t37 = sin(pkin(8));
t38 = cos(pkin(8));
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t18 = -t40 * t37 + t42 * t38;
t28 = t37 * pkin(4) + qJ(3);
t61 = 0.2e1 * t28;
t41 = sin(qJ(2));
t60 = -0.2e1 * t41;
t59 = 0.2e1 * t41;
t43 = cos(qJ(2));
t58 = 0.2e1 * t43;
t57 = 0.2e1 * qJ(3);
t39 = -pkin(2) - qJ(4);
t56 = -pkin(7) + t39;
t17 = t42 * t37 + t40 * t38;
t55 = t17 * t41;
t54 = t37 * t43;
t53 = t38 * t43;
t48 = -t41 * qJ(3) - pkin(1);
t15 = t39 * t43 + t48;
t29 = t41 * pkin(6);
t24 = t41 * pkin(3) + t29;
t7 = t38 * t15 + t37 * t24;
t30 = t43 * pkin(6);
t25 = t43 * pkin(3) + t30;
t27 = t37 ^ 2 + t38 ^ 2;
t35 = t41 ^ 2;
t50 = t43 ^ 2 + t35;
t49 = t43 * qJ(3);
t20 = t38 * t24;
t6 = -t37 * t15 + t20;
t3 = t7 * t37 + t6 * t38;
t47 = -t41 * pkin(2) + t49;
t46 = t39 * t41 + t49;
t44 = qJ(3) ^ 2;
t23 = -t43 * pkin(2) + t48;
t22 = t56 * t38;
t21 = t56 * t37;
t16 = t27 * t39;
t14 = pkin(4) * t53 + t25;
t13 = t18 * t41;
t11 = t17 * t43;
t10 = t18 * t43;
t9 = t42 * t21 + t40 * t22;
t8 = -t40 * t21 + t42 * t22;
t5 = -pkin(7) * t53 + t7;
t4 = t41 * pkin(4) + t20 + (pkin(7) * t43 - t15) * t37;
t2 = t40 * t4 + t42 * t5;
t1 = t42 * t4 - t40 * t5;
t12 = [1, 0, 0, t35, t41 * t58, 0, 0, 0, pkin(1) * t58, pkin(1) * t60, 0.2e1 * t50 * pkin(6), t23 * t58, t23 * t60, t50 * pkin(6) ^ 2 + t23 ^ 2, 0.2e1 * t25 * t53 + 0.2e1 * t6 * t41, -0.2e1 * t25 * t54 - 0.2e1 * t7 * t41, (t37 * t6 - t38 * t7) * t58, t25 ^ 2 + t6 ^ 2 + t7 ^ 2, t11 ^ 2, 0.2e1 * t11 * t10, -t11 * t59, -t10 * t59, t35, 0.2e1 * t1 * t41 + 0.2e1 * t14 * t10, -0.2e1 * t14 * t11 - 0.2e1 * t2 * t41; 0, 0, 0, 0, 0, t41, t43, 0, -t29, -t30, t47, t29, t30, t47 * pkin(6), t25 * t37 + t46 * t38, t25 * t38 - t46 * t37, -t3, t25 * qJ(3) + t3 * t39, -t11 * t18, -t18 * t10 + t11 * t17, t13, -t55, 0, t28 * t10 + t14 * t17 + t8 * t41, -t28 * t11 + t14 * t18 - t9 * t41; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t57, pkin(2) ^ 2 + t44, t37 * t57, t38 * t57, -0.2e1 * t16, t27 * t39 ^ 2 + t44, t18 ^ 2, -0.2e1 * t18 * t17, 0, 0, 0, t17 * t61, t18 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, t29, t38 * t41, -t37 * t41, 0, t3, 0, 0, 0, 0, 0, t13, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -t27, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t54, 0, t25, 0, 0, 0, 0, 0, t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t38, 0, qJ(3), 0, 0, 0, 0, 0, t17, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t10, t41, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
