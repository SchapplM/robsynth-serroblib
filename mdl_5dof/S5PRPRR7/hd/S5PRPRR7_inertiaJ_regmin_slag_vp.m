% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = sin(qJ(4));
t22 = 0.2e1 * t12 * pkin(4) + (2 * qJ(3));
t21 = 2 * qJ(3);
t11 = sin(qJ(5));
t20 = t11 * pkin(4);
t14 = cos(qJ(5));
t19 = t14 * pkin(4);
t15 = cos(qJ(4));
t5 = t11 * t15 + t14 * t12;
t18 = t11 * t12 - t14 * t15;
t17 = -pkin(2) - pkin(6);
t16 = cos(qJ(2));
t13 = sin(qJ(2));
t10 = t15 * t17;
t8 = -t15 * pkin(7) + t10;
t7 = (-pkin(7) + t17) * t12;
t4 = t5 * t16;
t3 = t18 * t16;
t2 = -t11 * t8 - t14 * t7;
t1 = -t11 * t7 + t14 * t8;
t6 = [1, 0, 0, 0, 0, 0, t13 ^ 2 + t16 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t16, -t13, -t16, t13, t16 * pkin(2) + t13 * qJ(3), 0, 0, 0, 0, 0, t13 * t12, t13 * t15, 0, 0, 0, 0, 0, t13 * t5, -t13 * t18; 0, 1, 0, 0, -0.2e1 * pkin(2), t21, pkin(2) ^ 2 + (qJ(3) ^ 2), t15 ^ 2, -0.2e1 * t15 * t12, 0, 0, 0, t12 * t21, t15 * t21, t18 ^ 2, 0.2e1 * t18 * t5, 0, 0, 0, t5 * t22, -t18 * t22; 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t16, t12 * t16, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t12, 0, t10, -t12 * t17, 0, 0, -t18, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t12, 0, 0, 0, 0, 0, -t18, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
