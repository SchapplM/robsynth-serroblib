% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = sin(qJ(4));
t25 = -0.2e1 * t9;
t11 = cos(qJ(4));
t24 = 0.2e1 * t11;
t10 = sin(qJ(3));
t12 = cos(qJ(3));
t7 = sin(pkin(8));
t8 = cos(pkin(8));
t2 = t10 * t7 - t12 * t8;
t23 = t2 * t9;
t22 = t9 * pkin(6);
t3 = t10 * t8 + t12 * t7;
t21 = t9 * t3;
t5 = t9 ^ 2;
t20 = t11 ^ 2 + t5;
t19 = t11 * pkin(6);
t18 = t11 * t3;
t17 = t2 * t11;
t16 = t20 * t3;
t15 = -t9 * pkin(4) + t11 * qJ(5);
t14 = t11 * pkin(4) + t9 * qJ(5);
t4 = -pkin(3) - t14;
t1 = [1, t7 ^ 2 + t8 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t3 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t17, t23, -t17, t16, -t23, pkin(6) * t16 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, t5, t9 * t24, 0, 0, 0, pkin(3) * t24, pkin(3) * t25, -0.2e1 * t4 * t11, 0.2e1 * t20 * pkin(6), t4 * t25, t20 * pkin(6) ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t18, -t21, 0, t18, t15 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t9, t11, 0, t9, t14; 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, -t22, -t19, -t22, t15, t19, t15 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
