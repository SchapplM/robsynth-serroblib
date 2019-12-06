% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t18 = sin(pkin(8));
t20 = cos(pkin(8));
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t4 = t18 * t22 - t20 * t24;
t30 = t4 ^ 2;
t14 = -t20 * pkin(2) - pkin(3);
t19 = cos(pkin(9));
t29 = -0.2e1 * t19 * pkin(4) + 0.2e1 * t14;
t11 = t18 * pkin(2) + qJ(4);
t28 = pkin(6) + t11;
t17 = sin(pkin(9));
t27 = t17 ^ 2 + t19 ^ 2;
t26 = t27 * t11;
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t8 = t23 * t17 + t21 * t19;
t6 = t21 * t17 - t23 * t19;
t7 = t18 * t24 + t20 * t22;
t3 = t7 ^ 2;
t2 = t28 * t19;
t1 = t28 * t17;
t5 = [1, 0, 0, 0, t3 + t30, 0, 0, 0, t27 * t3 + t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, t24, -t22, (t18 * t7 - t20 * t4) * pkin(2), -t4 * t19, t4 * t17, t27 * t7, t4 * t14 + t7 * t26, 0, 0, 0, 0, 0, t4 * t6, t4 * t8; 0, 1, 0, 0, (t18 ^ 2 + t20 ^ 2) * pkin(2) ^ 2, -0.2e1 * t14 * t19, 0.2e1 * t14 * t17, 0.2e1 * t26, t27 * t11 ^ 2 + t14 ^ 2, t8 ^ 2, -0.2e1 * t8 * t6, 0, 0, 0, t6 * t29, t8 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t19, t17, 0, t14, 0, 0, 0, 0, 0, t6, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8 * t7, t6 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t6, 0, -t23 * t1 - t21 * t2, t21 * t1 - t23 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
