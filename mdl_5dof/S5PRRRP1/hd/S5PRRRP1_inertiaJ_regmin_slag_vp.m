% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:12
% EndTime: 2019-12-05 16:40:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (72->32), mult. (134->51), div. (0->0), fcn. (118->4), ass. (0->24)
t16 = cos(qJ(4));
t25 = 0.2e1 * t16;
t14 = sin(qJ(4));
t24 = pkin(4) * t14;
t23 = sin(qJ(3)) * pkin(2);
t22 = t16 * pkin(4);
t21 = t16 * pkin(7);
t9 = pkin(7) + t23;
t20 = t16 * t9;
t19 = cos(qJ(3)) * pkin(2);
t10 = -pkin(3) - t19;
t18 = pkin(3) - t10;
t11 = -pkin(3) - t22;
t13 = t14 ^ 2;
t12 = t16 * qJ(5);
t8 = t14 * t25;
t7 = t12 + t21;
t6 = (-qJ(5) - pkin(7)) * t14;
t5 = t11 - t19;
t4 = t7 * t16;
t3 = t12 + t20;
t2 = (-qJ(5) - t9) * t14;
t1 = t3 * t16;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t3 + t16 * t2; 0, 1, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t23, t13, t8, 0, 0, 0, -0.2e1 * t10 * t16, 0.2e1 * t10 * t14, -0.2e1 * t2 * t14 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t7 + t16 * t6; 0, 0, 0, 0, 1, t19, -t23, t13, t8, 0, 0, 0, t18 * t16, -t18 * t14, t1 + t4 + (-t2 - t6) * t14, t5 * t11 + t2 * t6 + t3 * t7; 0, 0, 0, 0, 1, 0, 0, t13, t8, 0, 0, 0, pkin(3) * t25, -0.2e1 * pkin(3) * t14, -0.2e1 * t6 * t14 + 0.2e1 * t4, t11 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * t9, -t20, -t24, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * pkin(7), -t21, -t24, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t15;
