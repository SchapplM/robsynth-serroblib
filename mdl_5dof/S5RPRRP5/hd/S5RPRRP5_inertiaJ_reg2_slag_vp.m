% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:59
% DurationCPUTime: 0.54s
% Computational Cost: add. (188->35), mult. (366->64), div. (0->0), fcn. (299->6), ass. (0->40)
t31 = sin(qJ(4));
t27 = t31 ^ 2;
t33 = cos(qJ(4));
t28 = t33 ^ 2;
t18 = t27 + t28;
t30 = cos(pkin(8));
t47 = t30 * pkin(1);
t22 = pkin(2) + t47;
t32 = sin(qJ(3));
t34 = cos(qJ(3));
t29 = sin(pkin(8));
t48 = t29 * pkin(1);
t14 = t32 * t22 + t34 * t48;
t12 = pkin(7) + t14;
t49 = t18 * t12;
t53 = -0.2e1 * t31;
t52 = -0.2e1 * t33;
t51 = t49 * pkin(7);
t50 = t18 * t12 ^ 2;
t46 = t31 * pkin(7);
t45 = t33 * pkin(7);
t13 = t34 * t22 - t32 * t48;
t11 = -pkin(3) - t13;
t44 = pkin(3) - t11;
t37 = t33 * pkin(4) + t31 * qJ(5);
t16 = -pkin(3) - t37;
t3 = -t13 + t16;
t43 = -t16 - t3;
t42 = t31 * t12;
t41 = t31 * t33;
t40 = t33 * t12;
t39 = t18 * pkin(7) ^ 2;
t38 = t18 * pkin(7);
t17 = -t31 * pkin(4) + t33 * qJ(5);
t21 = -0.2e1 * t41;
t20 = 0.2e1 * t41;
t15 = 0.2e1 * t38;
t2 = 0.2e1 * t49;
t1 = t38 + t49;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t48, 0, (t29 ^ 2 + t30 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t14, 0, t13 ^ 2 + t14 ^ 2, t27, t20, 0, t28, 0, 0, t11 * t52, 0.2e1 * t11 * t31, t2, t11 ^ 2 + t50, t27, 0, t21, 0, 0, t28, t3 * t52, t2, t3 * t53, t3 ^ 2 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, -t14, 0, 0, t27, t20, 0, t28, 0, 0, t44 * t33, -t44 * t31, t1, -t11 * pkin(3) + t51, t27, 0, t21, 0, 0, t28, t43 * t33, t1, t43 * t31, t3 * t16 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, t20, 0, t28, 0, 0, 0.2e1 * pkin(3) * t33, pkin(3) * t53, t15, pkin(3) ^ 2 + t39, t27, 0, t21, 0, 0, t28, t16 * t52, t15, t16 * t53, t16 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t42, -t40, 0, 0, 0, t31, 0, 0, -t33, 0, -t42, t17, t40, t17 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t31, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t33, 0, -t46, -t45, 0, 0, 0, t31, 0, 0, -t33, 0, -t46, t17, t45, t17 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
