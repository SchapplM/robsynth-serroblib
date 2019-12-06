% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x13]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t15 = cos(qJ(5));
t19 = 0.2e1 * t15;
t11 = cos(pkin(9));
t9 = sin(pkin(9));
t18 = t11 ^ 2 + t9 ^ 2;
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t17 = t16 * t11 - t14 * t9;
t4 = t14 * t11 + t16 * t9;
t13 = sin(qJ(5));
t12 = cos(pkin(8));
t10 = sin(pkin(8));
t8 = t12 ^ 2;
t6 = t10 ^ 2;
t2 = t17 * t10;
t1 = t4 * t10;
t3 = [1, t6 + t8, t18 * t6 + t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t1, -t2, 0, 0, 0, 0, 0, -t1 * t15, t1 * t13; 0, 0, 0, 0, t17, -t4, 0, 0, 0, 0, 0, t17 * t15, -t17 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, t13 ^ 2, t13 * t19, 0, 0, 0, pkin(4) * t19, -0.2e1 * pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t15 - t13 * t2, t12 * t13 - t15 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t4, -t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13; 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * pkin(6), -t15 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
