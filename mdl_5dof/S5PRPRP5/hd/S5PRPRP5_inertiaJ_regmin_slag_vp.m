% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t18 = sin(qJ(4));
t27 = cos(qJ(4));
t8 = t27 * t16 + t18 * t17;
t31 = -0.2e1 * t8;
t11 = -t17 * pkin(3) - pkin(2);
t30 = 0.2e1 * t11;
t20 = cos(qJ(2));
t22 = -t18 * t16 + t27 * t17;
t29 = t20 * t22;
t28 = t20 * t8;
t26 = pkin(6) + qJ(3);
t25 = t16 ^ 2 + t17 ^ 2;
t24 = t26 * t16;
t23 = t25 * qJ(3);
t19 = sin(qJ(2));
t15 = t20 ^ 2;
t9 = t26 * t17;
t5 = t22 * t19;
t4 = t8 * t19;
t3 = -t18 * t24 + t27 * t9;
t2 = t18 * t9 + t27 * t24;
t1 = -pkin(4) * t22 - t8 * qJ(5) + t11;
t6 = [1, 0, 0, 0, 0, 0, 0, t25 * t19 ^ 2 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 ^ 2 + t5 ^ 2 + t15; 0, 0, t20, -t19, t20 * t17, -t20 * t16, t25 * t19, t20 * pkin(2) + t19 * t23, 0, 0, 0, 0, 0, t29, -t28, t29, t22 * t5 + t4 * t8, t28, -t20 * t1 + t4 * t2 + t5 * t3; 0, 1, 0, 0, 0.2e1 * pkin(2) * t17, -0.2e1 * pkin(2) * t16, 0.2e1 * t23, t25 * qJ(3) ^ 2 + pkin(2) ^ 2, t8 ^ 2, -t22 * t31, 0, 0, 0, -t22 * t30, t8 * t30, -0.2e1 * t1 * t22, 0.2e1 * t2 * t8 + 0.2e1 * t22 * t3, t1 * t31, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, -t17, t16, 0, -pkin(2), 0, 0, 0, 0, 0, -t22, t8, -t22, 0, -t8, t1; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, -t4, 0, t5, -t4 * pkin(4) + t5 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t22, 0, -t2, -t3, -t2, -pkin(4) * t8 + qJ(5) * t22, t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
