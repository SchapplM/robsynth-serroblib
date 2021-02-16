% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:04
% EndTime: 2021-01-15 12:17:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (215->48), mult. (387->92), div. (0->0), fcn. (451->6), ass. (0->38)
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t11 = -t21 * t24 + t22 * t26;
t12 = t21 * t26 + t22 * t24;
t44 = (t11 * t22 + t12 * t21) * pkin(3);
t10 = t12 ^ 2;
t9 = t11 ^ 2;
t43 = t10 + t9;
t16 = t24 * pkin(3) + qJ(2);
t42 = 0.2e1 * t16;
t41 = 0.2e1 * qJ(2);
t40 = -pkin(1) - pkin(6);
t39 = t21 * pkin(3);
t38 = t22 * pkin(3);
t23 = sin(qJ(5));
t37 = t23 * t11;
t36 = t23 * t12;
t25 = cos(qJ(5));
t35 = t23 * t25;
t34 = t25 * t11;
t7 = t25 * t12;
t33 = -qJ(4) + t40;
t31 = t33 * t26;
t13 = t33 * t24;
t4 = t21 * t13 - t22 * t31;
t6 = t22 * t13 + t21 * t31;
t30 = t4 * t11 - t6 * t12;
t15 = pkin(7) + t39;
t17 = -pkin(4) - t38;
t29 = t11 * t17 - t12 * t15;
t20 = t25 ^ 2;
t19 = t23 ^ 2;
t3 = t12 * pkin(4) - t11 * pkin(7) + t16;
t2 = t23 * t3 + t25 * t6;
t1 = -t23 * t6 + t25 * t3;
t5 = [1, 0, 0, -2 * pkin(1), t41, (pkin(1) ^ 2) + qJ(2) ^ 2, t26 ^ 2, -0.2e1 * t26 * t24, 0, 0, 0, t24 * t41, t26 * t41, t12 * t42, t11 * t42, 0.2e1 * t30, t16 ^ 2 + t4 ^ 2 + t6 ^ 2, t20 * t9, -0.2e1 * t9 * t35, 0.2e1 * t11 * t7, -0.2e1 * t11 * t36, t10, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t37, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t34; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t30, 0, 0, 0, 0, 0, -t43 * t23, -t43 * t25; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, t26 * t40, -t24 * t40, -t4, -t6, -t44, (t21 * t6 - t22 * t4) * pkin(3), t23 * t34, (-t19 + t20) * t11, t36, t7, 0, t29 * t23 - t4 * t25, t4 * t23 + t29 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, t11, -t12, 0, t44, 0, 0, 0, 0, 0, t34, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t38, -0.2e1 * t39, 0, (t21 ^ 2 + t22 ^ 2) * pkin(3) ^ 2, t19, 0.2e1 * t35, 0, 0, 0, -0.2e1 * t17 * t25, 0.2e1 * t17 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t16, 0, 0, 0, 0, 0, t7, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t37, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25, 0, -t23 * t15, -t25 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
