% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP3
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
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t20 = cos(pkin(8));
t15 = -t20 * pkin(1) - pkin(2);
t26 = cos(qJ(3));
t13 = -t26 * pkin(3) + t15;
t30 = 0.2e1 * t13;
t22 = sin(qJ(3));
t29 = 0.2e1 * t22;
t21 = sin(qJ(4));
t25 = cos(qJ(4));
t10 = t21 * t22 - t25 * t26;
t28 = t10 * pkin(4);
t27 = t21 * pkin(3);
t19 = sin(pkin(8));
t14 = t19 * pkin(1) + pkin(6);
t7 = (-pkin(7) - t14) * t22;
t24 = t26 * t14;
t8 = t26 * pkin(7) + t24;
t3 = -t21 * t8 + t25 * t7;
t4 = -t21 * t7 - t25 * t8;
t18 = t25 * pkin(3);
t17 = t18 + pkin(4);
t12 = t21 * t26 + t25 * t22;
t9 = t12 ^ 2;
t5 = t13 + t28;
t2 = -t10 * qJ(5) - t4;
t1 = -t12 * qJ(5) + t3;
t6 = [1, 0, 0, (t19 ^ 2 + t20 ^ 2) * pkin(1) ^ 2, t22 ^ 2, t26 * t29, 0, 0, 0, -0.2e1 * t15 * t26, t15 * t29, t9, -0.2e1 * t12 * t10, 0, 0, 0, t10 * t30, t12 * t30, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t10 + t2 * t12; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9; 0, 0, 0, 0, 0, 0, t22, t26, 0, -t22 * t14, -t24, 0, 0, t12, -t10, 0, t3, t4, -t10 * t27 - t17 * t12, t1 * t17 + t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t22, 0, 0, 0, 0, 0, -t10, -t12, 0, -t10 * t17 + t12 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, -0.2e1 * t27, 0, t21 ^ 2 * pkin(3) ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, t3, t4, -pkin(4) * t12, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t18, -t27, 0, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
