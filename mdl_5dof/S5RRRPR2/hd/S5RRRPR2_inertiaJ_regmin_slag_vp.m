% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:52
% EndTime: 2022-01-20 11:30:53
% DurationCPUTime: 0.21s
% Computational Cost: add. (171->40), mult. (314->62), div. (0->0), fcn. (294->8), ass. (0->33)
t33 = cos(qJ(5));
t40 = -0.2e1 * t33;
t26 = cos(qJ(2)) * pkin(1);
t24 = t26 + pkin(2);
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t38 = sin(qJ(2)) * pkin(1);
t14 = t34 * t24 - t31 * t38;
t11 = pkin(3) + t14;
t37 = t34 * t38;
t15 = t31 * t24 + t37;
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t5 = t28 * t11 + t29 * t15;
t39 = t31 * pkin(2);
t25 = t34 * pkin(2);
t23 = t25 + pkin(3);
t13 = t28 * t23 + t29 * t39;
t4 = t29 * t11 - t28 * t15;
t12 = t29 * t23 - t28 * t39;
t30 = sin(qJ(5));
t27 = t30 ^ 2;
t22 = -t29 * pkin(3) - pkin(4);
t21 = t28 * pkin(3) + pkin(8);
t20 = 0.2e1 * t30 * t33;
t17 = t22 * t30;
t10 = pkin(8) + t13;
t9 = -pkin(4) - t12;
t6 = t9 * t30;
t3 = pkin(8) + t5;
t2 = -pkin(4) - t4;
t1 = t2 * t30;
t7 = [1, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t38, 1, 0.2e1 * t14, -0.2e1 * t15, t4 ^ 2 + t5 ^ 2, t27, t20, 0, 0, 0, t2 * t40, 0.2e1 * t1; 0, 0, 0, 1, t26, -t38, 1, t14 + t25, -t37 + (-pkin(2) - t24) * t31, t4 * t12 + t5 * t13, t27, t20, 0, 0, 0, (-t2 - t9) * t33, t6 + t1; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t39, t12 ^ 2 + t13 ^ 2, t27, t20, 0, 0, 0, t9 * t40, 0.2e1 * t6; 0, 0, 0, 0, 0, 0, 1, t14, -t15, (t28 * t5 + t29 * t4) * pkin(3), t27, t20, 0, 0, 0, (-t2 - t22) * t33, t17 + t1; 0, 0, 0, 0, 0, 0, 1, t25, -t39, (t12 * t29 + t13 * t28) * pkin(3), t27, t20, 0, 0, 0, (-t22 - t9) * t33, t17 + t6; 0, 0, 0, 0, 0, 0, 1, 0, 0, (t28 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, t27, t20, 0, 0, 0, t22 * t40, 0.2e1 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t30 * t3, -t33 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t30 * t10, -t33 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t30 * t21, -t33 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
