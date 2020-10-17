% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRR1
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
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.15s
% Computational Cost: add. (51->24), mult. (109->34), div. (0->0), fcn. (153->8), ass. (0->23)
t16 = cos(qJ(5));
t23 = 0.2e1 * t16;
t17 = cos(qJ(4));
t20 = t17 * pkin(3);
t9 = -pkin(4) - t20;
t22 = pkin(4) - t9;
t14 = sin(qJ(4));
t21 = t14 * pkin(3);
t11 = sin(pkin(9));
t12 = cos(pkin(9));
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t5 = -t15 * t11 + t18 * t12;
t6 = t18 * t11 + t15 * t12;
t2 = t14 * t6 - t17 * t5;
t19 = t2 * t16;
t13 = sin(qJ(5));
t10 = t13 ^ 2;
t8 = pkin(7) + t21;
t7 = t13 * t23;
t3 = t14 * t5 + t17 * t6;
t1 = t2 * t13;
t4 = [1, t11 ^ 2 + t12 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t5, -t6, 0, -t2, -t3, 0, 0, 0, 0, 0, -t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 1, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t21, t10, t7, 0, 0, 0, -0.2e1 * t9 * t16, 0.2e1 * t9 * t13; 0, 0, 0, 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, t20, -t21, t10, t7, 0, 0, 0, t22 * t16, -t22 * t13; 0, 0, 0, 0, 0, 1, 0, 0, t10, t7, 0, 0, 0, pkin(4) * t23, -0.2e1 * pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t3, -t16 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * t8, -t16 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * pkin(7), -t16 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
