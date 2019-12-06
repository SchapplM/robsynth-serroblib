% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t18 = cos(pkin(9));
t15 = -t18 * pkin(2) - pkin(3);
t23 = cos(qJ(4));
t30 = -0.2e1 * t23 * pkin(4) + 0.2e1 * t15;
t20 = sin(qJ(4));
t29 = 0.2e1 * t20;
t19 = sin(qJ(5));
t28 = t19 * pkin(4);
t22 = cos(qJ(5));
t27 = t22 * pkin(4);
t17 = sin(pkin(9));
t14 = t17 * pkin(2) + pkin(6);
t26 = pkin(7) + t14;
t11 = t19 * t23 + t22 * t20;
t10 = t19 * t20 - t22 * t23;
t24 = cos(qJ(2));
t21 = sin(qJ(2));
t9 = t17 * t24 + t18 * t21;
t7 = t17 * t21 - t18 * t24;
t6 = t26 * t23;
t5 = t26 * t20;
t4 = t19 * t5 - t22 * t6;
t3 = -t19 * t6 - t22 * t5;
t2 = t10 * t9;
t1 = t11 * t9;
t8 = [1, 0, 0, 0, t7 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t24, -t21, (t17 * t9 - t18 * t7) * pkin(2), 0, 0, 0, 0, 0, -t7 * t23, t7 * t20, 0, 0, 0, 0, 0, t7 * t10, t7 * t11; 0, 1, 0, 0, (t17 ^ 2 + t18 ^ 2) * pkin(2) ^ 2, t20 ^ 2, t23 * t29, 0, 0, 0, -0.2e1 * t15 * t23, t15 * t29, t11 ^ 2, -0.2e1 * t11 * t10, 0, 0, 0, t10 * t30, t11 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t9, -t23 * t9, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, t20, t23, 0, -t20 * t14, -t23 * t14, 0, 0, t11, -t10, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t20, 0, 0, 0, 0, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
