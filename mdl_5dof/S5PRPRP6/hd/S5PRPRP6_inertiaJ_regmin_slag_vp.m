% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:18
% EndTime: 2019-12-05 15:41:19
% DurationCPUTime: 0.16s
% Computational Cost: add. (57->30), mult. (100->42), div. (0->0), fcn. (93->4), ass. (0->22)
t17 = -pkin(2) - pkin(6);
t15 = cos(qJ(4));
t10 = t15 ^ 2;
t13 = sin(qJ(4));
t5 = t13 ^ 2 + t10;
t2 = t5 * t17;
t23 = -0.2e1 * t15;
t22 = 2 * qJ(3);
t16 = cos(qJ(2));
t21 = t13 * t16;
t20 = t13 * t17;
t14 = sin(qJ(2));
t19 = t14 * t15;
t18 = t15 * t16;
t4 = t15 * pkin(4) + t13 * qJ(5);
t11 = t16 ^ 2;
t9 = t14 ^ 2;
t7 = t15 * t17;
t6 = t14 * t13;
t3 = t13 * pkin(4) - t15 * qJ(5) + qJ(3);
t1 = t5 * t16;
t8 = [1, 0, 0, 0, 0, 0, t9 + t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t11 + t9; 0, 0, t16, -t14, -t16, t14, t16 * pkin(2) + t14 * qJ(3), 0, 0, 0, 0, 0, t6, t19, t6, t1, -t19, t14 * t3 - t16 * t2; 0, 1, 0, 0, -0.2e1 * pkin(2), t22, pkin(2) ^ 2 + (qJ(3) ^ 2), t10, t13 * t23, 0, 0, 0, t13 * t22, t15 * t22, 0.2e1 * t3 * t13, -0.2e1 * t2, t3 * t23, t5 * t17 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, t2; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t21, -t18, 0, -t21, -t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, t7, -t20, t7, -t4, t20, t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, t15, 0, t13, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
