% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t14 = sin(qJ(4));
t10 = t14 ^ 2;
t16 = cos(qJ(4));
t12 = t16 ^ 2;
t21 = t10 + t12;
t28 = t21 * pkin(6);
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t18 = -pkin(1) - pkin(2);
t7 = t17 * qJ(2) + t15 * t18;
t4 = -pkin(6) + t7;
t20 = t21 * t4;
t27 = -0.2e1 * t14;
t26 = 0.2e1 * t16;
t5 = t15 * qJ(2) - t17 * t18;
t3 = pkin(3) + t5;
t25 = pkin(3) + t3;
t24 = t14 * t16;
t23 = t17 * t14;
t22 = t17 * t16;
t1 = t21 * t15;
t13 = t17 ^ 2;
t11 = t15 ^ 2;
t8 = 0.2e1 * t24;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t5, 0.2e1 * t7, 0, t5 ^ 2 + t7 ^ 2, t10, t8, 0, t12, 0, 0, t3 * t26, t3 * t27, -0.2e1 * t20, t21 * t4 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t17, t15, 0, t7 * t15 - t5 * t17, 0, 0, 0, 0, 0, 0, -t22, t23, -t1, t4 * t1 - t3 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t11 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t5, -t7, 0, 0, -t10, -0.2e1 * t24, 0, -t12, 0, 0, -t25 * t16, t25 * t14, t20 - t28, -t3 * pkin(3) + pkin(6) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, t1, t17 * pkin(3) + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t10, t8, 0, t12, 0, 0, pkin(3) * t26, pkin(3) * t27, 0.2e1 * t28, t21 * pkin(6) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t16, 0, -t14 * t4, -t16 * t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t15, -t16 * t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t16, 0, -t14 * pkin(6), -t16 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
