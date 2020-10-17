% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR16_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:34
% DurationCPUTime: 0.55s
% Computational Cost: add. (213->58), mult. (350->91), div. (0->0), fcn. (309->4), ass. (0->45)
t30 = sin(qJ(3));
t24 = t30 ^ 2;
t32 = cos(qJ(3));
t26 = t32 ^ 2;
t14 = t24 + t26;
t9 = t30 * pkin(3) - qJ(4) * t32 + qJ(2);
t49 = -0.2e1 * t9;
t48 = 0.2e1 * qJ(2);
t47 = 0.2e1 * qJ(4);
t29 = sin(qJ(5));
t16 = t29 * t30;
t46 = t29 * t32;
t31 = cos(qJ(5));
t45 = t31 * t29;
t17 = t31 * t30;
t44 = t32 * t30;
t33 = -pkin(3) - pkin(7);
t43 = t32 * t33;
t34 = -pkin(1) - pkin(6);
t42 = t32 * t34;
t41 = t14 * t34 ^ 2;
t23 = t29 ^ 2;
t25 = t31 ^ 2;
t13 = t23 + t25;
t40 = t30 * qJ(4);
t39 = 0.2e1 * t44;
t38 = t29 * t17;
t11 = (pkin(4) - t34) * t32;
t4 = pkin(7) * t30 + t9;
t2 = t11 * t31 - t29 * t4;
t3 = t11 * t29 + t31 * t4;
t1 = t2 * t31 + t29 * t3;
t12 = pkin(3) * t32 + t40;
t37 = t40 - t43;
t36 = qJ(2) ^ 2;
t35 = qJ(4) ^ 2;
t20 = t30 * t34;
t18 = t31 * t32;
t15 = -0.2e1 * t44;
t10 = -pkin(4) * t30 + t20;
t8 = t14 * t34;
t7 = t13 * t33;
t6 = t13 * t32;
t5 = 0.2e1 * t8;
t19 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t48, (pkin(1) ^ 2) + t36, t26, t15, 0, t24, 0, 0, t30 * t48, t32 * t48, -t5, t36 + t41, 0, 0, 0, t26, t15, t24, -t5, t30 * t49, t32 * t49, t9 ^ 2 + t41, t23 * t24, 0.2e1 * t24 * t45, t29 * t39, t25 * t24, t31 * t39, t26, -0.2e1 * t10 * t17 + 0.2e1 * t2 * t32, 0.2e1 * t10 * t16 - 0.2e1 * t3 * t32, 0.2e1 * (-t2 * t29 + t3 * t31) * t30, t10 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t14, t8, 0, 0, 0, 0, 0, 0, -t14, 0, 0, t8, 0, 0, 0, 0, 0, 0, -t14 * t31, t14 * t29, 0, -t1 * t32 + t10 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t26 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t30, 0, t42, -t20, 0, 0, 0, -t32, t30, 0, 0, 0, -t12, -t42, t20, t12 * t34, t38, (-t23 + t25) * t30, t18, -t38, -t46, 0, t10 * t29 - t37 * t31, t10 * t31 + t37 * t29, -t1, t10 * qJ(4) + t1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, t12, 0, 0, 0, 0, 0, 0, t16, t17, t6, -t13 * t43 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t47, pkin(3) ^ 2 + t35, t25, -0.2e1 * t45, 0, t23, 0, 0, t29 * t47, t31 * t47, -0.2e1 * t7, t13 * t33 ^ 2 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, -t42, 0, 0, 0, 0, 0, 0, t18, -t46, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t13, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t17, t32, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, 0, t31 * t33, -t29 * t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t19;
