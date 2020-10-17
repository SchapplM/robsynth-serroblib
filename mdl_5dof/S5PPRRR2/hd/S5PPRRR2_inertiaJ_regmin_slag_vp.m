% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (57->26), mult. (134->44), div. (0->0), fcn. (187->8), ass. (0->24)
t19 = cos(qJ(4));
t25 = -0.2e1 * t19 * pkin(4) - (2 * pkin(3));
t24 = 0.2e1 * t19;
t23 = pkin(6) + pkin(7);
t15 = sin(qJ(5));
t22 = t15 * pkin(4);
t18 = cos(qJ(5));
t21 = t18 * pkin(4);
t16 = sin(qJ(4));
t8 = t15 * t19 + t18 * t16;
t7 = t15 * t16 - t18 * t19;
t20 = cos(qJ(3));
t17 = sin(qJ(3));
t14 = cos(pkin(9));
t13 = sin(pkin(9));
t10 = t23 * t19;
t9 = t23 * t16;
t6 = t20 * t13 + t17 * t14;
t5 = t17 * t13 - t20 * t14;
t4 = -t18 * t10 + t15 * t9;
t3 = -t15 * t10 - t18 * t9;
t2 = t7 * t6;
t1 = t8 * t6;
t11 = [1, t13 ^ 2 + t14 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t19, t5 * t16, 0, 0, 0, 0, 0, t5 * t7, t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, t16 ^ 2, t16 * t24, 0, 0, 0, pkin(3) * t24, -0.2e1 * pkin(3) * t16, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t25, t8 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t6, -t19 * t6, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t16, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * pkin(6), -t19 * pkin(6), 0, 0, t8, -t7, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
