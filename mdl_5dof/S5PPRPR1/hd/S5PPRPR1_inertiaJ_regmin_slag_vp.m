% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t14 = cos(pkin(9));
t24 = -0.2e1 * t14 * pkin(4) - (2 * pkin(3));
t23 = pkin(6) + qJ(4);
t12 = sin(pkin(9));
t22 = t12 ^ 2 + t14 ^ 2;
t21 = t22 * qJ(4);
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t4 = t18 * t12 + t16 * t14;
t2 = t16 * t12 - t18 * t14;
t19 = cos(qJ(3));
t17 = sin(qJ(3));
t15 = cos(pkin(8));
t13 = sin(pkin(8));
t7 = t23 * t14;
t6 = t23 * t12;
t5 = t19 * t13 + t17 * t15;
t3 = t17 * t13 - t19 * t15;
t1 = [1, t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, t22 * t5 ^ 2 + t3 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t3, -t5, -t3 * t14, t3 * t12, t22 * t5, -t3 * pkin(3) + t5 * t21, 0, 0, 0, 0, 0, t3 * t2, t3 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t14, -0.2e1 * pkin(3) * t12, 0.2e1 * t21, t22 * qJ(4) ^ 2 + (pkin(3) ^ 2), t4 ^ 2, -0.2e1 * t4 * t2, 0, 0, 0, t2 * t24, t4 * t24; 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t14, t12, 0, -pkin(3), 0, 0, 0, 0, 0, t2, t4; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t5, t2 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t2, 0, -t16 * t7 - t18 * t6, t16 * t6 - t18 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
