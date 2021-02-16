% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRP1
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
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:48:17
% EndTime: 2021-01-15 14:48:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (60->27), mult. (124->45), div. (0->0), fcn. (149->6), ass. (0->22)
t14 = cos(qJ(4));
t23 = 0.2e1 * t14;
t12 = sin(qJ(4));
t8 = t12 ^ 2;
t22 = t14 ^ 2 + t8;
t10 = sin(pkin(8));
t11 = cos(pkin(8));
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t4 = t15 * t10 + t13 * t11;
t21 = t12 * t4;
t20 = t14 * pkin(4);
t19 = t14 * t4;
t3 = t13 * t10 - t15 * t11;
t18 = t3 * t14;
t17 = -qJ(5) - pkin(6);
t5 = t17 * t12;
t6 = t17 * t14;
t16 = -t5 * t12 - t6 * t14;
t7 = -pkin(3) - t20;
t1 = t3 * t12;
t2 = [1, t10 ^ 2 + t11 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t4 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t18, t1, -t18, t1, t22 * t4, t16 * t4 + t3 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t6 + t14 * t5; 0, 0, 1, 0, 0, t8, t12 * t23, 0, 0, 0, pkin(3) * t23, -0.2e1 * pkin(3) * t12, -0.2e1 * t7 * t14, 0.2e1 * t7 * t12, 0.2e1 * t16, t5 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t19, -t21, -t19, 0, -pkin(4) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, t14, -t12, 0, t20; 0, 0, 0, 0, 0, 0, 0, t12, t14, 0, -t12 * pkin(6), -t14 * pkin(6), t5, t6, -t12 * pkin(4), t5 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t12, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
