% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:55
% DurationCPUTime: 0.12s
% Computational Cost: add. (41->17), mult. (76->32), div. (0->0), fcn. (94->6), ass. (0->14)
t10 = sin(pkin(8));
t6 = t10 * pkin(2) + qJ(4);
t18 = 0.2e1 * t6;
t11 = cos(pkin(8));
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t1 = t10 * t13 - t11 * t15;
t3 = t10 * t15 + t11 * t13;
t17 = t1 ^ 2 + t3 ^ 2;
t9 = -t11 * pkin(2) - pkin(3);
t14 = cos(qJ(5));
t12 = sin(qJ(5));
t5 = -pkin(6) + t9;
t2 = [1, 0, 0, 0, t17, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, t15, -t13, (-t1 * t11 + t10 * t3) * pkin(2), t1, t3, t1 * t9 + t3 * t6, 0, 0, 0, 0, 0, t3 * t12, t3 * t14; 0, 1, 0, 0, (t10 ^ 2 + t11 ^ 2) * pkin(2) ^ 2, 0.2e1 * t9, t18, t6 ^ 2 + t9 ^ 2, t14 ^ 2, -0.2e1 * t14 * t12, 0, 0, 0, t12 * t18, t14 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t1, -t12 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, 0, t14 * t5, -t12 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t2;
