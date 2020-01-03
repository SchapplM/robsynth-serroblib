% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP2
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
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = cos(qJ(4));
t31 = 0.2e1 * t22;
t19 = cos(pkin(8));
t14 = t19 * pkin(1) + pkin(2);
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t18 = sin(pkin(8));
t29 = pkin(1) * t18;
t7 = t23 * t14 - t21 * t29;
t5 = -pkin(3) - t7;
t30 = pkin(3) - t5;
t20 = sin(qJ(4));
t28 = pkin(4) * t20;
t27 = t22 * pkin(4);
t26 = t22 * pkin(7);
t8 = -t21 * t14 - t23 * t29;
t6 = pkin(7) - t8;
t25 = t22 * t6;
t17 = t20 ^ 2;
t16 = t22 * qJ(5);
t15 = -pkin(3) - t27;
t13 = t20 * t31;
t11 = t16 + t26;
t10 = (-qJ(5) - pkin(7)) * t20;
t9 = t11 * t22;
t4 = t5 - t27;
t3 = t16 + t25;
t2 = (-qJ(5) - t6) * t20;
t1 = t3 * t22;
t12 = [1, 0, 0, (t18 ^ 2 + t19 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t7, 0.2e1 * t8, t17, t13, 0, 0, 0, -0.2e1 * t5 * t22, 0.2e1 * t5 * t20, -0.2e1 * t2 * t20 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t22 + t3 * t20; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 + t17; 0, 0, 0, 0, 1, t7, t8, t17, t13, 0, 0, 0, t30 * t22, -t30 * t20, t1 + t9 + (-t10 - t2) * t20, t2 * t10 + t3 * t11 + t4 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t10 + t20 * t11; 0, 0, 0, 0, 1, 0, 0, t17, t13, 0, 0, 0, pkin(3) * t31, -0.2e1 * pkin(3) * t20, -0.2e1 * t10 * t20 + 0.2e1 * t9, t10 ^ 2 + t11 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t20 * t6, -t25, -t28, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t20 * pkin(7), -t26, -t28, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t12;
