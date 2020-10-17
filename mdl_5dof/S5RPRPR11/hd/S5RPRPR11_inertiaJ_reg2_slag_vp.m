% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:55
% EndTime: 2019-12-31 18:27:56
% DurationCPUTime: 0.46s
% Computational Cost: add. (361->56), mult. (680->98), div. (0->0), fcn. (782->6), ass. (0->38)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t33 = sin(qJ(3));
t45 = cos(qJ(3));
t19 = t33 * t30 - t45 * t31;
t21 = t45 * t30 + t33 * t31;
t26 = -t31 * pkin(2) - pkin(1);
t37 = t21 * qJ(4) - t26;
t46 = -pkin(3) - pkin(4);
t4 = t46 * t19 + t37;
t50 = 0.2e1 * t4;
t49 = t19 ^ 2;
t48 = 0.2e1 * t26;
t47 = 0.2e1 * t31;
t44 = t21 * t19;
t43 = pkin(6) + qJ(2);
t28 = t30 ^ 2;
t29 = t31 ^ 2;
t42 = t28 + t29;
t39 = t43 * t31;
t40 = t43 * t30;
t13 = t33 * t39 + t45 * t40;
t15 = -t33 * t40 + t45 * t39;
t41 = t13 ^ 2 + t15 ^ 2;
t38 = 0.2e1 * t13 * t21 - 0.2e1 * t15 * t19;
t36 = -t21 * pkin(7) + t13;
t34 = cos(qJ(5));
t32 = sin(qJ(5));
t24 = t34 * qJ(4) + t32 * t46;
t22 = t32 * qJ(4) - t34 * t46;
t17 = t21 ^ 2;
t11 = t32 * t19 + t34 * t21;
t9 = -t34 * t19 + t32 * t21;
t8 = t19 * pkin(3) - t37;
t6 = t19 * pkin(7) + t15;
t3 = t32 * t36 + t34 * t6;
t1 = t32 * t6 - t34 * t36;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, t30 * t47, 0, t29, 0, 0, pkin(1) * t47, -0.2e1 * pkin(1) * t30, 0.2e1 * t42 * qJ(2), t42 * qJ(2) ^ 2 + pkin(1) ^ 2, t17, -0.2e1 * t44, 0, t49, 0, 0, t19 * t48, t21 * t48, t38, t26 ^ 2 + t41, t17, 0, 0.2e1 * t44, 0, 0, t49, 0.2e1 * t8 * t19, t38, -0.2e1 * t8 * t21, t8 ^ 2 + t41, t11 ^ 2, -0.2e1 * t11 * t9, 0, t9 ^ 2, 0, 0, t9 * t50, t11 * t50, 0.2e1 * t1 * t11 - 0.2e1 * t3 * t9, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t30, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t19, t21, 0, t26, 0, 0, 0, 0, 0, 0, t19, 0, -t21, t8, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, 0, -t13, -t15, 0, 0, 0, t21, 0, 0, t19, 0, -t13, -pkin(3) * t21 - t19 * qJ(4), t15, -t13 * pkin(3) + t15 * qJ(4), 0, 0, -t11, 0, t9, 0, t1, t3, t22 * t11 - t24 * t9, t1 * t22 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t22, 0.2e1 * t24, 0, t22 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t11 - t32 * t9, -t1 * t34 + t3 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t34, t32, 0, -t22 * t34 + t24 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, -t1, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t22, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
