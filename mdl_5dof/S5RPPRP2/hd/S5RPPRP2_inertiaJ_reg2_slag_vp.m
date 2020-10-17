% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:25
% DurationCPUTime: 0.37s
% Computational Cost: add. (208->32), mult. (385->62), div. (0->0), fcn. (426->6), ass. (0->34)
t26 = cos(pkin(8));
t28 = sin(qJ(4));
t24 = sin(pkin(8));
t38 = cos(qJ(4));
t33 = t38 * t24;
t16 = t28 * t26 + t33;
t13 = t16 ^ 2;
t37 = t28 * t24;
t14 = -t38 * t26 + t37;
t44 = t14 ^ 2;
t45 = t13 + t44;
t27 = cos(pkin(7));
t40 = t27 * pkin(1);
t21 = -pkin(2) - t40;
t17 = -t26 * pkin(3) + t21;
t43 = 0.2e1 * t17;
t42 = 0.2e1 * t24;
t25 = sin(pkin(7));
t41 = t25 * pkin(1);
t19 = qJ(3) + t41;
t39 = pkin(6) + t19;
t8 = t16 * t14;
t22 = t24 ^ 2;
t23 = t26 ^ 2;
t36 = t22 + t23;
t10 = t39 * t26;
t5 = t28 * t10 + t39 * t33;
t7 = t38 * t10 - t39 * t37;
t35 = t5 ^ 2 + t7 ^ 2;
t34 = t5 * t14 + t7 * t16;
t31 = -t14 * pkin(4) + t16 * qJ(5);
t30 = -0.2e1 * t7 * t14 + 0.2e1 * t5 * t16;
t3 = t17 - t31;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t41, 0, (t25 ^ 2 + t27 ^ 2) * pkin(1) ^ 2, t22, t26 * t42, 0, t23, 0, 0, -0.2e1 * t21 * t26, t21 * t42, 0.2e1 * t36 * t19, t36 * t19 ^ 2 + t21 ^ 2, t13, -0.2e1 * t8, 0, t44, 0, 0, t14 * t43, t16 * t43, t30, t17 ^ 2 + t35, t13, 0, 0.2e1 * t8, 0, 0, t44, 0.2e1 * t3 * t14, t30, -0.2e1 * t3 * t16, t3 ^ 2 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, 0, t21, 0, 0, 0, 0, 0, 0, t14, t16, 0, t17, 0, 0, 0, 0, 0, 0, t14, 0, -t16, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, -t5, -t7, 0, 0, 0, t16, 0, 0, t14, 0, -t5, -pkin(4) * t16 - t14 * qJ(5), t7, -t5 * pkin(4) + t7 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t16, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, t16, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
