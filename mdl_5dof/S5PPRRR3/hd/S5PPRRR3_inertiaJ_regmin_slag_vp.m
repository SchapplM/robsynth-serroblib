% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (58->28), mult. (146->54), div. (0->0), fcn. (201->8), ass. (0->30)
t21 = cos(qJ(4));
t31 = -0.2e1 * t21 * pkin(4) - (2 * pkin(3));
t30 = 0.2e1 * t21;
t29 = pkin(6) + pkin(7);
t17 = sin(qJ(5));
t28 = t17 * pkin(4);
t20 = cos(qJ(5));
t27 = t20 * pkin(4);
t15 = sin(pkin(9));
t19 = sin(qJ(3));
t26 = t19 * t15;
t25 = t21 * t19;
t18 = sin(qJ(4));
t22 = cos(qJ(3));
t24 = t22 * t18;
t23 = t22 * t21;
t10 = t17 * t21 + t20 * t18;
t9 = t17 * t18 - t20 * t21;
t16 = cos(pkin(9));
t12 = t29 * t21;
t11 = t29 * t18;
t8 = t15 * t23 - t16 * t18;
t7 = -t15 * t24 - t16 * t21;
t6 = t9 * t19;
t5 = t10 * t19;
t4 = t17 * t11 - t20 * t12;
t3 = -t20 * t11 - t17 * t12;
t2 = -t17 * t7 - t20 * t8;
t1 = -t17 * t8 + t20 * t7;
t13 = [1, t15 ^ 2 + t16 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t26, -t22 * t15, 0, 0, 0, 0, 0, -t15 * t25, t18 * t26, 0, 0, 0, 0, 0, t9 * t26, t10 * t26; 0, 0, 0, t22, -t19, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, -t22 * t9, -t22 * t10; 0, 0, 1, 0, 0, t18 ^ 2, t18 * t30, 0, 0, 0, pkin(3) * t30, -0.2e1 * pkin(3) * t18, t10 ^ 2, -0.2e1 * t10 * t9, 0, 0, 0, t9 * t31, t10 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t19, -t25, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, t18, t21, 0, -t18 * pkin(6), -t21 * pkin(6), 0, 0, t10, -t9, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t13;
