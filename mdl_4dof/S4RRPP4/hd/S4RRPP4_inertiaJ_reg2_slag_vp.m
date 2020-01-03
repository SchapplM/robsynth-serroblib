% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t12 = sin(qJ(2));
t10 = t12 ^ 2;
t13 = cos(qJ(2));
t11 = t13 ^ 2;
t27 = t10 + t11;
t26 = -0.2e1 * t12;
t25 = 0.2e1 * t13;
t14 = pkin(2) + pkin(3);
t24 = t27 * pkin(5) ^ 2;
t23 = t12 * pkin(5);
t22 = t12 * t13;
t21 = t13 * qJ(3);
t20 = t12 * qJ(3) + pkin(1);
t19 = -t12 * pkin(2) + t21;
t17 = qJ(3) ^ 2;
t16 = 0.2e1 * qJ(3);
t9 = t13 * pkin(5);
t6 = -0.2e1 * t22;
t5 = -t13 * qJ(4) + t9;
t4 = (pkin(5) - qJ(4)) * t12;
t3 = -t13 * pkin(2) - t20;
t2 = 0.2e1 * t27 * pkin(5);
t1 = t14 * t13 + t20;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t10, 0.2e1 * t22, 0, t11, 0, 0, pkin(1) * t25, pkin(1) * t26, t2, pkin(1) ^ 2 + t24, t10, 0, t6, 0, 0, t11, -0.2e1 * t3 * t13, t2, t3 * t26, t3 ^ 2 + t24, t10, t6, 0, t11, 0, 0, t1 * t25, 0.2e1 * t1 * t12, -0.2e1 * t4 * t12 - 0.2e1 * t5 * t13, t1 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t13, 0, -t23, -t9, 0, 0, 0, t12, 0, 0, -t13, 0, -t23, t19, t9, t19 * pkin(5), 0, 0, -t12, 0, t13, 0, -t4, t5, t14 * t12 - t21, t5 * qJ(3) - t4 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t16, pkin(2) ^ 2 + t17, 0, 0, 0, 0, 0, 1, 0.2e1 * t14, t16, 0, t14 ^ 2 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -1, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t12, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
