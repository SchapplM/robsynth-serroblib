% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x13]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (25->22), mult. (62->37), div. (0->0), fcn. (89->8), ass. (0->21)
t13 = cos(qJ(5));
t22 = 0.2e1 * t13;
t7 = sin(pkin(9));
t8 = sin(pkin(8));
t21 = t7 * t8;
t9 = cos(pkin(9));
t20 = t8 * t9;
t19 = t7 ^ 2 + t9 ^ 2;
t11 = sin(qJ(5));
t12 = sin(qJ(4));
t18 = t11 * t12;
t17 = t13 * t12;
t14 = cos(qJ(4));
t16 = t14 * t11;
t15 = t14 * t13;
t10 = cos(pkin(8));
t6 = t10 ^ 2;
t4 = t8 ^ 2;
t2 = -t10 * t12 + t14 * t20;
t1 = t10 * t14 + t12 * t20;
t3 = [1, t4 + t6, t19 * t4 + t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t1, -t2, 0, 0, 0, 0, 0, -t1 * t13, t1 * t11; 0, 0, 0, 0, -t12 * t7, -t14 * t7, 0, 0, 0, 0, 0, -t7 * t17, t7 * t18; 0, 0, 0, 0, t14, -t12, 0, 0, 0, 0, 0, t15, -t16; 0, 0, 0, 1, 0, 0, t11 ^ 2, t11 * t22, 0, 0, 0, pkin(4) * t22, -0.2e1 * pkin(4) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t11 + t13 * t21, -t11 * t21 - t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t9 - t7 * t16, t11 * t9 - t7 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t17; 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, 0, -t11 * pkin(6), -t13 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
