% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:42
% DurationCPUTime: 0.29s
% Computational Cost: add. (184->46), mult. (313->80), div. (0->0), fcn. (252->4), ass. (0->41)
t48 = sin(qJ(2)) * pkin(1);
t21 = pkin(7) + t48;
t33 = sin(qJ(3));
t29 = t33 ^ 2;
t35 = cos(qJ(3));
t42 = t35 ^ 2 + t29;
t44 = t42 * t21;
t23 = t35 * qJ(4);
t57 = -t33 * pkin(3) + t23;
t56 = -0.2e1 * t33;
t55 = -0.2e1 * t35;
t54 = 0.2e1 * t35;
t31 = pkin(3) + qJ(5);
t41 = -t33 * qJ(4) - pkin(2);
t4 = -t31 * t35 + t41;
t47 = cos(qJ(2)) * pkin(1);
t1 = t4 - t47;
t53 = -t1 - t4;
t16 = t33 * t21;
t49 = t33 * pkin(4);
t6 = t16 + t49;
t18 = t35 * t21;
t28 = t35 * pkin(4);
t7 = t18 + t28;
t52 = t6 * t33 + t7 * t35;
t25 = t33 * pkin(7);
t13 = t25 + t49;
t27 = t35 * pkin(7);
t14 = t27 + t28;
t51 = t13 * t33 + t14 * t35;
t22 = -pkin(2) - t47;
t46 = pkin(2) - t22;
t11 = -t35 * pkin(3) + t41;
t5 = t11 - t47;
t45 = t11 + t5;
t43 = t42 * pkin(7);
t38 = qJ(4) ^ 2;
t37 = 0.2e1 * qJ(4);
t19 = t33 * t54;
t10 = -t31 * t33 + t23;
t2 = [1, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t48, t29, t19, 0, 0, 0, t22 * t55, 0.2e1 * t22 * t33, 0.2e1 * t44, t5 * t54, t5 * t56, t42 * t21 ^ 2 + t5 ^ 2, 0.2e1 * t52, t1 * t56, t1 * t55, t1 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 1, t47, -t48, t29, t19, 0, 0, 0, t46 * t35, -t46 * t33, t43 + t44, t45 * t35, -t45 * t33, pkin(7) * t44 + t5 * t11, t51 + t52, t53 * t33, t53 * t35, t1 * t4 + t6 * t13 + t7 * t14; 0, 0, 0, 1, 0, 0, t29, t19, 0, 0, 0, pkin(2) * t54, pkin(2) * t56, 0.2e1 * t43, t11 * t54, t11 * t56, t42 * pkin(7) ^ 2 + t11 ^ 2, 0.2e1 * t51, t4 * t56, t4 * t55, t13 ^ 2 + t14 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, 0, -t16, -t18, t57, t16, t18, t57 * t21, t10, t7, -t6, t7 * qJ(4) - t6 * t31; 0, 0, 0, 0, 0, 0, 0, 0, t33, t35, 0, -t25, -t27, t57, t25, t27, t57 * pkin(7), t10, t14, -t13, t14 * qJ(4) - t13 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t37, pkin(3) ^ 2 + t38, 0, t37, 0.2e1 * t31, t31 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, t16, t33, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, t25, t33, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -1, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
