% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t19 = sin(qJ(3));
t32 = t19 * pkin(2);
t12 = pkin(6) + t32;
t18 = sin(qJ(4));
t16 = t18 ^ 2;
t21 = cos(qJ(4));
t17 = t21 ^ 2;
t27 = t16 + t17;
t36 = t27 * t12;
t20 = sin(qJ(2));
t22 = cos(qJ(3));
t23 = cos(qJ(2));
t4 = t19 * t20 - t22 * t23;
t35 = t4 ^ 2;
t34 = 0.2e1 * t21;
t31 = t22 * pkin(2);
t30 = t4 * t21;
t13 = -pkin(3) - t31;
t29 = pkin(3) - t13;
t28 = pkin(6) * t27;
t6 = t19 * t23 + t22 * t20;
t1 = t27 * t6;
t9 = t18 * t34;
t3 = t6 ^ 2;
t2 = t4 * t18;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 ^ 2 + t23 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 + t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t3 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, 0, (t19 * t6 - t22 * t4) * pkin(2), 0, 0, 0, 0, 0, 0, -t30, t2, t1, t4 * t13 + t36 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t32, 0, (t19 ^ 2 + t22 ^ 2) * pkin(2) ^ 2, t16, t9, 0, t17, 0, 0, -0.2e1 * t13 * t21, 0.2e1 * t13 * t18, 0.2e1 * t36, t27 * t12 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t6, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t2, t1, -t4 * pkin(3) + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t32, 0, 0, t16, t9, 0, t17, 0, 0, t29 * t21, -t29 * t18, t28 + t36, -t13 * pkin(3) + pkin(6) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t16, t9, 0, t17, 0, 0, pkin(3) * t34, -0.2e1 * pkin(3) * t18, 0.2e1 * t28, t27 * pkin(6) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t6, -t21 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t21, 0, -t18 * t12, -t21 * t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t21, 0, -t18 * pkin(6), -t21 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
