% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:06:26
% EndTime: 2021-01-15 12:06:27
% DurationCPUTime: 0.26s
% Computational Cost: add. (220->47), mult. (427->99), div. (0->0), fcn. (487->8), ass. (0->38)
t22 = sin(pkin(9));
t24 = cos(pkin(9));
t27 = sin(qJ(3));
t38 = cos(qJ(3));
t12 = t22 * t27 - t24 * t38;
t43 = t12 ^ 2;
t25 = cos(pkin(8));
t19 = -t25 * pkin(1) - pkin(2);
t15 = -t38 * pkin(3) + t19;
t42 = 0.2e1 * t15;
t41 = 0.2e1 * t27;
t40 = t22 * pkin(3);
t39 = t24 * pkin(3);
t26 = sin(qJ(5));
t8 = t26 * t12;
t14 = t22 * t38 + t24 * t27;
t37 = t26 * t14;
t28 = cos(qJ(5));
t36 = t26 * t28;
t35 = t28 * t14;
t23 = sin(pkin(8));
t34 = t23 * pkin(1) + pkin(6);
t17 = pkin(7) + t40;
t18 = -pkin(4) - t39;
t33 = -t12 * t17 + t14 * t18;
t32 = (-qJ(4) - t34) * t27;
t31 = t38 * t34;
t21 = t28 ^ 2;
t20 = t26 ^ 2;
t11 = t14 ^ 2;
t10 = t38 * qJ(4) + t31;
t9 = t28 * t12;
t6 = t24 * t10 + t22 * t32;
t4 = t22 * t10 - t24 * t32;
t3 = t12 * pkin(4) - t14 * pkin(7) + t15;
t2 = t26 * t3 + t28 * t6;
t1 = -t26 * t6 + t28 * t3;
t5 = [1, 0, 0, (t23 ^ 2 + t25 ^ 2) * pkin(1) ^ 2, t27 ^ 2, t38 * t41, 0, 0, 0, -0.2e1 * t19 * t38, t19 * t41, t12 * t42, t14 * t42, -0.2e1 * t6 * t12 + 0.2e1 * t4 * t14, t15 ^ 2 + t4 ^ 2 + t6 ^ 2, t21 * t11, -0.2e1 * t11 * t36, 0.2e1 * t12 * t35, -0.2e1 * t12 * t37, t43, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t37, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t12 + t6 * t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t27, t38, 0, -t27 * t34, -t31, -t4, -t6, (-t12 * t22 - t14 * t24) * pkin(3), (t22 * t6 - t24 * t4) * pkin(3), t26 * t35, (-t20 + t21) * t14, t8, t9, 0, t33 * t26 - t4 * t28, t4 * t26 + t33 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t27, -t12, -t14, 0, (-t12 * t24 + t14 * t22) * pkin(3), 0, 0, 0, 0, 0, -t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, -0.2e1 * t40, 0, (t22 ^ 2 + t24 ^ 2) * pkin(3) ^ 2, t20, 0.2e1 * t36, 0, 0, 0, -0.2e1 * t18 * t28, 0.2e1 * t18 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14, 0, t15, 0, 0, 0, 0, 0, t9, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t37, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, -t26 * t17, -t28 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
