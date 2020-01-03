% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR10
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
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t15 = sin(qJ(2));
t33 = -0.2e1 * t15;
t17 = cos(qJ(2));
t32 = 0.2e1 * t17;
t31 = 2 * qJ(3);
t18 = -pkin(2) - pkin(6);
t14 = sin(qJ(4));
t30 = t14 * t15;
t29 = t14 * t17;
t28 = t15 * t17;
t16 = cos(qJ(4));
t27 = t16 * t14;
t26 = t16 * t17;
t11 = t15 ^ 2;
t13 = t17 ^ 2;
t25 = t11 + t13;
t24 = t17 * qJ(3);
t23 = -0.2e1 * t28;
t22 = -t15 * qJ(3) - pkin(1);
t21 = -t15 * pkin(2) + t24;
t20 = t15 * t18 + t24;
t12 = t16 ^ 2;
t10 = t14 ^ 2;
t9 = t17 * pkin(5);
t8 = t15 * pkin(5);
t7 = t16 * t15;
t6 = t17 * pkin(3) + t9;
t5 = t15 * pkin(3) + t8;
t4 = -t17 * pkin(2) + t22;
t3 = t18 * t17 + t22;
t2 = t14 * t5 + t16 * t3;
t1 = -t14 * t3 + t16 * t5;
t19 = [1, 0, 0, t11, 0.2e1 * t28, 0, 0, 0, pkin(1) * t32, pkin(1) * t33, 0.2e1 * t25 * pkin(5), t4 * t32, t4 * t33, t25 * pkin(5) ^ 2 + t4 ^ 2, t10 * t13, 0.2e1 * t13 * t27, t14 * t23, t16 * t23, t11, 0.2e1 * t1 * t15 + 0.2e1 * t6 * t26, -0.2e1 * t2 * t15 - 0.2e1 * t6 * t29; 0, 0, 0, 0, 0, t15, t17, 0, -t8, -t9, t21, t8, t9, t21 * pkin(5), -t14 * t26, (t10 - t12) * t17, t7, -t30, 0, t6 * t14 + t20 * t16, -t20 * t14 + t6 * t16; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t31, pkin(2) ^ 2 + (qJ(3) ^ 2), t12, -0.2e1 * t27, 0, 0, 0, t14 * t31, t16 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, t8, 0, 0, 0, 0, 0, t7, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t26, t15, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t16 * t18, -t14 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t19;
