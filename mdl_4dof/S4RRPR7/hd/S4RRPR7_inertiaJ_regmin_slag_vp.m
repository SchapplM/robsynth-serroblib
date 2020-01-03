% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t24 = cos(qJ(2));
t33 = 0.2e1 * t24;
t21 = sin(qJ(4));
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t22 = sin(qJ(2));
t9 = t19 * t22 - t20 * t24;
t32 = t21 * t9;
t10 = t19 * t24 + t20 * t22;
t31 = t21 * t10;
t23 = cos(qJ(4));
t30 = t21 * t23;
t29 = t23 * t10;
t28 = -qJ(3) - pkin(5);
t16 = -t24 * pkin(2) - pkin(1);
t27 = t28 * t22;
t14 = t19 * pkin(2) + pkin(6);
t15 = -t20 * pkin(2) - pkin(3);
t26 = t10 * t15 - t14 * t9;
t18 = t23 ^ 2;
t17 = t21 ^ 2;
t12 = t28 * t24;
t8 = t10 ^ 2;
t7 = t23 * t9;
t6 = -t20 * t12 + t19 * t27;
t4 = -t19 * t12 - t20 * t27;
t3 = t9 * pkin(3) - t10 * pkin(6) + t16;
t2 = t21 * t3 + t23 * t6;
t1 = -t21 * t6 + t23 * t3;
t5 = [1, 0, 0, t22 ^ 2, t22 * t33, 0, 0, 0, pkin(1) * t33, -0.2e1 * pkin(1) * t22, 0.2e1 * t4 * t10 - 0.2e1 * t6 * t9, t16 ^ 2 + t4 ^ 2 + t6 ^ 2, t18 * t8, -0.2e1 * t8 * t30, 0.2e1 * t9 * t29, -0.2e1 * t9 * t31, t9 ^ 2, 0.2e1 * t1 * t9 + 0.2e1 * t4 * t31, -0.2e1 * t2 * t9 + 0.2e1 * t4 * t29; 0, 0, 0, 0, 0, t22, t24, 0, -t22 * pkin(5), -t24 * pkin(5), (-t10 * t20 - t19 * t9) * pkin(2), (t19 * t6 - t20 * t4) * pkin(2), t21 * t29, (-t17 + t18) * t10, t32, t7, 0, t26 * t21 - t4 * t23, t4 * t21 + t26 * t23; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t19 ^ 2 + t20 ^ 2) * pkin(2) ^ 2, t17, 0.2e1 * t30, 0, 0, 0, -0.2e1 * t15 * t23, 0.2e1 * t15 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, t7, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t31, t9, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, -t21 * t14, -t23 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
