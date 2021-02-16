% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR2
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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:35
% EndTime: 2021-01-15 21:23:36
% DurationCPUTime: 0.32s
% Computational Cost: add. (465->45), mult. (906->95), div. (0->0), fcn. (1103->8), ass. (0->45)
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t26 = t36 * t40 - t37 * t43;
t27 = t36 * t43 + t37 * t40;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t16 = t42 * t26 + t39 * t27;
t34 = -t43 * pkin(2) - pkin(1);
t21 = t26 * pkin(3) + t34;
t53 = 0.2e1 * t16 * pkin(4) + 0.2e1 * t21;
t52 = 0.2e1 * t21;
t51 = 0.2e1 * t34;
t50 = 0.2e1 * t43;
t49 = t36 * pkin(2);
t48 = t37 * pkin(2);
t38 = sin(qJ(5));
t47 = t38 * pkin(4);
t33 = pkin(3) + t48;
t24 = t39 * t33 + t42 * t49;
t41 = cos(qJ(5));
t46 = t41 * t24;
t45 = -qJ(3) - pkin(6);
t29 = t45 * t40;
t30 = t45 * t43;
t18 = t37 * t29 + t36 * t30;
t13 = -t27 * pkin(7) + t18;
t19 = t36 * t29 - t37 * t30;
t14 = -t26 * pkin(7) + t19;
t5 = t42 * t13 - t39 * t14;
t23 = t42 * t33 - t39 * t49;
t22 = pkin(4) + t23;
t11 = t41 * t22 - t38 * t24;
t6 = -t39 * t13 - t42 * t14;
t35 = t41 * pkin(4);
t17 = -t39 * t26 + t42 * t27;
t12 = -t38 * t22 - t46;
t8 = -t38 * t16 + t41 * t17;
t7 = t41 * t16 + t38 * t17;
t4 = -t16 * pkin(8) - t6;
t3 = -t17 * pkin(8) + t5;
t2 = -t38 * t3 - t41 * t4;
t1 = t41 * t3 - t38 * t4;
t9 = [1, 0, 0, t40 ^ 2, t40 * t50, 0, 0, 0, pkin(1) * t50, -0.2e1 * pkin(1) * t40, t26 * t51, t27 * t51, -0.2e1 * t18 * t27 - 0.2e1 * t19 * t26, t18 ^ 2 + t19 ^ 2 + t34 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t52, t17 * t52, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t53, t8 * t53; 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(6), -t43 * pkin(6), t18, -t19, (-t26 * t36 - t27 * t37) * pkin(2), (t18 * t37 + t19 * t36) * pkin(2), 0, 0, t17, -t16, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t48, -0.2e1 * t49, 0, (t36 ^ 2 + t37 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t23, -0.2e1 * t24, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, t34, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t23, -t24, 0, 0, 0, 0, 1, t11 + t35, -t46 + (-pkin(4) - t22) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
