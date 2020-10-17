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
% MM_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (47->22), mult. (101->40), div. (0->0), fcn. (117->6), ass. (0->19)
t13 = cos(qJ(4));
t20 = 0.2e1 * t13;
t11 = sin(qJ(4));
t7 = t11 ^ 2;
t19 = t13 ^ 2 + t7;
t10 = cos(pkin(8));
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t9 = sin(pkin(8));
t3 = t12 * t10 + t14 * t9;
t18 = t11 * t3;
t17 = t13 * pkin(4);
t16 = -qJ(5) - pkin(6);
t4 = t16 * t11;
t5 = t16 * t13;
t15 = -t4 * t11 - t5 * t13;
t6 = -pkin(3) - t17;
t2 = -t14 * t10 + t12 * t9;
t1 = [1, t10 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t3 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t2 * t13, t2 * t11, t19 * t3, t15 * t3 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t5 + t13 * t4; 0, 0, 1, 0, 0, t7, t11 * t20, 0, 0, 0, pkin(3) * t20, -0.2e1 * pkin(3) * t11, 0.2e1 * t15, t4 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t13 * t3, 0, -pkin(4) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11, 0, t17; 0, 0, 0, 0, 0, 0, 0, t11, t13, 0, -t11 * pkin(6), -t13 * pkin(6), -t11 * pkin(4), t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
