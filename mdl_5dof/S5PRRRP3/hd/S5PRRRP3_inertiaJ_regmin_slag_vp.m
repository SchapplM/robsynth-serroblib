% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP3
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
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (113->32), mult. (231->56), div. (0->0), fcn. (266->4), ass. (0->23)
t22 = cos(qJ(3));
t15 = -t22 * pkin(3) - pkin(2);
t25 = 0.2e1 * t15;
t17 = sin(qJ(4));
t18 = sin(qJ(3));
t21 = cos(qJ(4));
t7 = t17 * t18 - t21 * t22;
t24 = t7 * pkin(4);
t23 = t17 * pkin(3);
t20 = 0.2e1 * t22;
t19 = t22 * pkin(6);
t11 = (-pkin(6) - pkin(7)) * t18;
t12 = t22 * pkin(7) + t19;
t3 = t21 * t11 - t17 * t12;
t4 = -t17 * t11 - t21 * t12;
t16 = t21 * pkin(3);
t14 = t16 + pkin(4);
t9 = t17 * t22 + t21 * t18;
t6 = t9 ^ 2;
t5 = t15 + t24;
t2 = -t7 * qJ(5) - t4;
t1 = -t9 * qJ(5) + t3;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t1 + t9 * t2; 0, 1, 0, 0, t18 ^ 2, t18 * t20, 0, 0, 0, pkin(2) * t20, -0.2e1 * pkin(2) * t18, t6, -0.2e1 * t9 * t7, 0, 0, 0, t7 * t25, t9 * t25, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t18, 0, 0, 0, 0, 0, -t7, -t9, 0, -t7 * t14 + t9 * t23; 0, 0, 0, 0, 0, 0, t18, t22, 0, -t18 * pkin(6), -t19, 0, 0, t9, -t7, 0, t3, t4, -t14 * t9 - t7 * t23, t1 * t14 + t2 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t16, -0.2e1 * t23, 0, t17 ^ 2 * pkin(3) ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t7, 0, t3, t4, -pkin(4) * t9, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t16, -t23, 0, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
