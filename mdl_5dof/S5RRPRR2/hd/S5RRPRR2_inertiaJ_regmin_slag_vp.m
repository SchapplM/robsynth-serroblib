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
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (453->44), mult. (876->90), div. (0->0), fcn. (1071->8), ass. (0->43)
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t26 = -t37 * t41 + t38 * t44;
t27 = t37 * t44 + t38 * t41;
t40 = sin(qJ(4));
t43 = cos(qJ(4));
t16 = -t43 * t26 + t40 * t27;
t35 = -t44 * pkin(2) - pkin(1);
t21 = -t26 * pkin(3) + t35;
t52 = 0.2e1 * t16 * pkin(4) + 0.2e1 * t21;
t51 = 0.2e1 * t21;
t50 = 0.2e1 * t44;
t49 = pkin(2) * t37;
t39 = sin(qJ(5));
t48 = t39 * pkin(4);
t34 = t38 * pkin(2) + pkin(3);
t24 = t40 * t34 + t43 * t49;
t42 = cos(qJ(5));
t47 = t42 * t24;
t46 = -qJ(3) - pkin(6);
t31 = t46 * t41;
t32 = t46 * t44;
t19 = t37 * t31 - t38 * t32;
t18 = t38 * t31 + t37 * t32;
t13 = -t27 * pkin(7) + t18;
t14 = t26 * pkin(7) + t19;
t5 = t43 * t13 - t40 * t14;
t23 = t43 * t34 - t40 * t49;
t22 = pkin(4) + t23;
t11 = t42 * t22 - t39 * t24;
t6 = -t40 * t13 - t43 * t14;
t36 = t42 * pkin(4);
t17 = t40 * t26 + t43 * t27;
t12 = -t39 * t22 - t47;
t8 = -t39 * t16 + t42 * t17;
t7 = t42 * t16 + t39 * t17;
t4 = -t16 * pkin(8) - t6;
t3 = -t17 * pkin(8) + t5;
t2 = -t39 * t3 - t42 * t4;
t1 = t42 * t3 - t39 * t4;
t9 = [1, 0, 0, t41 ^ 2, t41 * t50, 0, 0, 0, pkin(1) * t50, -0.2e1 * pkin(1) * t41, -0.2e1 * t18 * t27 + 0.2e1 * t19 * t26, t18 ^ 2 + t19 ^ 2 + t35 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t51, t17 * t51, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t52, t8 * t52; 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(6), -t44 * pkin(6), (t26 * t37 - t27 * t38) * pkin(2), (t18 * t38 + t19 * t37) * pkin(2), 0, 0, t17, -t16, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t37 ^ 2 + t38 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t23, -0.2e1 * t24, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t23, -t24, 0, 0, 0, 0, 1, t11 + t36, -t47 + (-pkin(4) - t22) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
