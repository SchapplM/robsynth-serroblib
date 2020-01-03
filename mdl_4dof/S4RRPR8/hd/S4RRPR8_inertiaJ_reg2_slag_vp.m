% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t23 = sin(qJ(2));
t20 = t23 ^ 2;
t25 = cos(qJ(2));
t21 = t25 ^ 2;
t36 = t20 + t21;
t28 = t23 * qJ(3) + pkin(1);
t33 = pkin(2) + pkin(3);
t4 = t33 * t25 + t28;
t35 = 0.2e1 * t4;
t34 = -0.2e1 * t23;
t32 = t23 * pkin(5);
t31 = t23 * t25;
t30 = t36 * pkin(5) ^ 2;
t29 = (pkin(5) - pkin(6)) * t23;
t27 = -t23 * pkin(2) + t25 * qJ(3);
t24 = cos(qJ(4));
t22 = sin(qJ(4));
t19 = t25 * pkin(5);
t14 = -t25 * pkin(6) + t19;
t13 = -t25 * pkin(2) - t28;
t12 = t24 * qJ(3) - t22 * t33;
t10 = t22 * qJ(3) + t24 * t33;
t9 = 0.2e1 * t36 * pkin(5);
t7 = -t25 * t22 + t23 * t24;
t5 = t23 * t22 + t25 * t24;
t3 = t24 * t14 + t22 * t29;
t1 = t22 * t14 - t24 * t29;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, 0.2e1 * t31, 0, t21, 0, 0, 0.2e1 * pkin(1) * t25, pkin(1) * t34, t9, pkin(1) ^ 2 + t30, t20, 0, -0.2e1 * t31, 0, 0, t21, -0.2e1 * t13 * t25, t9, t13 * t34, t13 ^ 2 + t30, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t35, t7 * t35, 0.2e1 * t1 * t7 - 0.2e1 * t3 * t5, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t25, 0, -t32, -t19, 0, 0, 0, t23, 0, 0, -t25, 0, -t32, t27, t19, t27 * pkin(5), 0, 0, -t7, 0, t5, 0, t1, t3, t10 * t7 - t12 * t5, t1 * t10 + t3 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t10, 0.2e1 * t12, 0, t10 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t5 - t24 * t7, -t1 * t24 + t3 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t24, t22, 0, -t10 * t24 + t12 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, -t1, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t10, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
