% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t21 = sin(qJ(3));
t22 = sin(qJ(2));
t23 = cos(qJ(3));
t24 = cos(qJ(2));
t8 = t21 * t22 - t23 * t24;
t38 = t8 ^ 2;
t17 = -t24 * pkin(2) - pkin(1);
t37 = 0.2e1 * t17;
t36 = 0.2e1 * t24;
t35 = -pkin(6) - pkin(5);
t10 = t21 * t24 + t23 * t22;
t34 = t10 * t8;
t33 = t23 * pkin(2);
t19 = t22 ^ 2;
t20 = t24 ^ 2;
t32 = t19 + t20;
t12 = t35 * t24;
t30 = t35 * t22;
t4 = -t21 * t12 - t23 * t30;
t6 = -t23 * t12 + t21 * t30;
t31 = t4 ^ 2 + t6 ^ 2;
t29 = 0.2e1 * t4 * t10 - 0.2e1 * t6 * t8;
t26 = 2 * pkin(3);
t25 = 2 * qJ(4);
t18 = t21 * pkin(2);
t15 = pkin(3) + t33;
t13 = t18 + qJ(4);
t7 = t10 ^ 2;
t2 = t8 * pkin(3) - t10 * qJ(4) + t17;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t19, t22 * t36, 0, t20, 0, 0, pkin(1) * t36, -0.2e1 * pkin(1) * t22, 0.2e1 * t32 * pkin(5), t32 * pkin(5) ^ 2 + pkin(1) ^ 2, t7, -0.2e1 * t34, 0, t38, 0, 0, t8 * t37, t10 * t37, t29, t17 ^ 2 + t31, t7, 0, 0.2e1 * t34, 0, 0, t38, 0.2e1 * t2 * t8, t29, -0.2e1 * t2 * t10, t2 ^ 2 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t24, 0, -t22 * pkin(5), -t24 * pkin(5), 0, 0, 0, 0, t10, 0, -t8, 0, -t4, -t6, (-t10 * t23 - t21 * t8) * pkin(2), (t21 * t6 - t23 * t4) * pkin(2), 0, t10, 0, 0, t8, 0, -t4, -t15 * t10 - t13 * t8, t6, t6 * t13 - t4 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t18, 0, (t21 ^ 2 + t23 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t15, 0, 0.2e1 * t13, t13 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, 0, -t4, -t6, 0, 0, 0, t10, 0, 0, t8, 0, -t4, -pkin(3) * t10 - t8 * qJ(4), t6, -t4 * pkin(3) + t6 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t18, 0, 0, 0, 0, 0, 1, 0, 0, t26 + t33, 0, t25 + t18, t15 * pkin(3) + t13 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t26, 0, t25, pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
