% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP3
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
% MM_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:41
% DurationCPUTime: 0.14s
% Computational Cost: add. (68->26), mult. (131->47), div. (0->0), fcn. (145->6), ass. (0->23)
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t4 = t13 * t16 - t14 * t18;
t26 = t4 ^ 2;
t15 = sin(qJ(4));
t25 = 0.2e1 * t15;
t6 = t13 * t18 + t14 * t16;
t24 = t15 * t6;
t17 = cos(qJ(4));
t23 = t17 * pkin(4);
t9 = t13 * pkin(2) + pkin(6);
t22 = qJ(5) + t9;
t11 = t15 ^ 2;
t21 = t17 ^ 2 + t11;
t10 = -t14 * pkin(2) - pkin(3);
t1 = t22 * t15;
t2 = t22 * t17;
t20 = t1 * t15 + t2 * t17;
t7 = t10 - t23;
t3 = t6 ^ 2;
t5 = [1, 0, 0, 0, t3 + t26, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t3 + t26; 0, 0, t18, -t16, (t13 * t6 - t14 * t4) * pkin(2), 0, 0, 0, 0, 0, -t4 * t17, t4 * t15, t21 * t6, t20 * t6 + t4 * t7; 0, 1, 0, 0, (t13 ^ 2 + t14 ^ 2) * pkin(2) ^ 2, t11, t17 * t25, 0, 0, 0, -0.2e1 * t10 * t17, t10 * t25, 0.2e1 * t20, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t17 + t2 * t15; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t17 * t6, 0, -pkin(4) * t24; 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, -t15 * t9, -t17 * t9, -t15 * pkin(4), -t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
