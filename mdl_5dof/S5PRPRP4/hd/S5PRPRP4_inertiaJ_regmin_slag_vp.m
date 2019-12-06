% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP4
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
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t15 = sin(qJ(2));
t17 = cos(qJ(2));
t3 = t12 * t15 - t13 * t17;
t31 = t3 ^ 2;
t14 = sin(qJ(4));
t30 = 0.2e1 * t14;
t16 = cos(qJ(4));
t29 = -0.2e1 * t16;
t5 = t12 * t17 + t13 * t15;
t28 = t14 * t5;
t8 = t12 * pkin(2) + pkin(6);
t27 = t14 * t8;
t26 = t16 * t5;
t25 = t16 * t8;
t24 = t3 * t14;
t23 = t3 * t16;
t10 = t14 ^ 2;
t22 = t16 ^ 2 + t10;
t9 = -t13 * pkin(2) - pkin(3);
t21 = t22 * t8;
t20 = t16 * pkin(4) + t14 * qJ(5);
t19 = -t14 * pkin(4) + t16 * qJ(5);
t2 = t5 ^ 2;
t1 = -t20 + t9;
t4 = [1, 0, 0, 0, t2 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t2 + t31; 0, 0, t17, -t15, (t12 * t5 - t13 * t3) * pkin(2), 0, 0, 0, 0, 0, -t23, t24, -t23, t22 * t5, -t24, t3 * t1 + t5 * t21; 0, 1, 0, 0, (t12 ^ 2 + t13 ^ 2) * pkin(2) ^ 2, t10, t16 * t30, 0, 0, 0, t9 * t29, t9 * t30, t1 * t29, 0.2e1 * t21, -0.2e1 * t1 * t14, t22 * t8 ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t26, -t28, 0, t26, t19 * t5; 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t27, -t25, -t27, t19, t25, t19 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, t16, 0, t14, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
