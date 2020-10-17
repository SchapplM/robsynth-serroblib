% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR4
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:50
% DurationCPUTime: 0.45s
% Computational Cost: add. (272->58), mult. (335->80), div. (0->0), fcn. (280->6), ass. (0->39)
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t37 = -pkin(2) - pkin(3);
t13 = t32 * qJ(3) + t31 * t37;
t10 = -pkin(7) + t13;
t33 = sin(qJ(5));
t29 = t33 ^ 2;
t35 = cos(qJ(5));
t30 = t35 ^ 2;
t44 = t29 + t30;
t42 = t44 * t10;
t50 = t44 * t31;
t49 = -0.2e1 * t33;
t48 = 0.2e1 * t35;
t36 = cos(qJ(2));
t46 = t36 * pkin(1);
t23 = pkin(2) + t46;
t17 = -pkin(3) - t23;
t15 = t32 * t17;
t34 = sin(qJ(2));
t26 = t34 * pkin(1);
t20 = t26 + qJ(3);
t4 = t31 * t20 - t15;
t2 = pkin(4) + t4;
t22 = t32 * t37;
t11 = t31 * qJ(3) - t22;
t9 = pkin(4) + t11;
t47 = t2 + t9;
t45 = t32 * t35;
t6 = t31 * t17 + t32 * t20;
t3 = -pkin(7) + t6;
t43 = t44 * t3;
t39 = 0.2e1 * pkin(2);
t38 = 0.2e1 * qJ(3);
t28 = t32 ^ 2;
t27 = t31 ^ 2;
t19 = t32 * t33;
t18 = t33 * t48;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t26, 0, (t34 ^ 2 + t36 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t23, 0, 0.2e1 * t20, t20 ^ 2 + t23 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t4, 0.2e1 * t6, 0, t4 ^ 2 + t6 ^ 2, t29, t18, 0, t30, 0, 0, t2 * t48, t2 * t49, -0.2e1 * t43, t44 * t3 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t26, 0, 0, 0, 0, 0, 1, 0, 0, t39 + t46, 0, t38 + t26, t23 * pkin(2) + t20 * qJ(3), 0, 0, 0, 0, 0, 1, -t15 - t22 + (qJ(3) + t20) * t31, t13 + t6, 0, t4 * t11 + t6 * t13, t29, t18, 0, t30, 0, 0, t47 * t35, -t47 * t33, -t43 - t42, t2 * t9 + t3 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t39, 0, t38, pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t13, 0, t11 ^ 2 + t13 ^ 2, t29, t18, 0, t30, 0, 0, t9 * t48, t9 * t49, -0.2e1 * t42, t44 * t10 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t23, 0, 0, 0, 0, 0, 0, -t32, t31, 0, t6 * t31 - t4 * t32, 0, 0, 0, 0, 0, 0, -t45, t19, -t50, -t2 * t32 + t3 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t32, t31, 0, -t11 * t32 + t13 * t31, 0, 0, 0, 0, 0, 0, -t45, t19, -t50, t10 * t50 - t9 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 + t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t27 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t35, 0, -t33 * t3, -t35 * t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t35, 0, -t33 * t10, -t35 * t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t31, -t35 * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
