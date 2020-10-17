% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.43s
% Computational Cost: add. (165->41), mult. (327->68), div. (0->0), fcn. (347->4), ass. (0->32)
t24 = sin(pkin(7));
t22 = t24 ^ 2;
t25 = cos(pkin(7));
t23 = t25 ^ 2;
t42 = t22 + t23;
t26 = sin(qJ(4));
t27 = cos(qJ(4));
t12 = t24 * t26 + t25 * t27;
t41 = t12 ^ 2;
t40 = 0.2e1 * t12;
t39 = -0.2e1 * t24;
t14 = t24 * t27 - t25 * t26;
t38 = t14 * t12;
t37 = t24 * t25;
t36 = -pkin(6) + qJ(2);
t35 = t42 * qJ(2) ^ 2;
t17 = t36 * t25;
t31 = t36 * t24;
t5 = t26 * t17 - t27 * t31;
t7 = t27 * t17 + t26 * t31;
t34 = t5 ^ 2 + t7 ^ 2;
t33 = t7 * t26 - t5 * t27;
t32 = -t26 * t12 - t27 * t14;
t30 = t24 * qJ(3) + pkin(1);
t29 = -0.2e1 * t7 * t12 + 0.2e1 * t5 * t14;
t9 = (pkin(2) + pkin(3)) * t25 + t30;
t18 = t26 ^ 2 + t27 ^ 2;
t16 = -t25 * pkin(2) - t30;
t15 = 0.2e1 * t42 * qJ(2);
t10 = t14 ^ 2;
t1 = t12 * pkin(4) - t14 * qJ(5) + t9;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, 0.2e1 * t37, 0, t23, 0, 0, 0.2e1 * pkin(1) * t25, pkin(1) * t39, t15, pkin(1) ^ 2 + t35, t22, 0, -0.2e1 * t37, 0, 0, t23, -0.2e1 * t16 * t25, t15, t16 * t39, t16 ^ 2 + t35, t10, -0.2e1 * t38, 0, t41, 0, 0, t9 * t40, 0.2e1 * t9 * t14, t29, t9 ^ 2 + t34, t10, 0, 0.2e1 * t38, 0, 0, t41, t1 * t40, t29, -0.2e1 * t1 * t14, t1 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t25, 0, -t24, t16, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, -t9, 0, 0, 0, 0, 0, 0, -t12, 0, t14, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t24 * qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, t32, t33, 0, 0, 0, 0, 0, 0, 0, t32, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, 0, -t5, -t7, 0, 0, 0, t14, 0, 0, t12, 0, -t5, -t14 * pkin(4) - t12 * qJ(5), t7, -t5 * pkin(4) + t7 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t26, t27 * pkin(4) + t26 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
