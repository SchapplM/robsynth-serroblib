% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (61->26), mult. (130->36), div. (0->0), fcn. (124->6), ass. (0->25)
t16 = sin(qJ(5));
t28 = pkin(4) * t16;
t17 = sin(qJ(4));
t27 = t17 * pkin(3);
t26 = sin(qJ(3)) * pkin(2);
t19 = cos(qJ(5));
t13 = cos(qJ(3)) * pkin(2);
t11 = t13 + pkin(3);
t20 = cos(qJ(4));
t22 = -t20 * t11 + t17 * t26;
t2 = -pkin(4) + t22;
t25 = t2 * t19;
t12 = t20 * pkin(3);
t10 = -t12 - pkin(4);
t24 = t10 * t19;
t23 = t20 * t26;
t5 = -t17 * t11 - t23;
t15 = t16 ^ 2;
t14 = pkin(4) * t19;
t9 = pkin(8) + t27;
t8 = 0.2e1 * t16 * t19;
t7 = t10 * t16;
t3 = pkin(8) - t5;
t1 = t2 * t16;
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t26, 1, -0.2e1 * t22, 0.2e1 * t5, t15, t8, 0, 0, 0, -0.2e1 * t25, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t13, -t26, 1, t12 - t22, -t23 + (-pkin(3) - t11) * t17, t15, t8, 0, 0, 0, (-t10 - t2) * t19, t7 + t1; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t12, -0.2e1 * t27, t15, t8, 0, 0, 0, -0.2e1 * t24, 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, -t22, t5, t15, t8, 0, 0, 0, t14 - t25, t1 - t28; 0, 0, 0, 0, 0, 0, 0, 1, t12, -t27, t15, t8, 0, 0, 0, t14 - t24, t7 - t28; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, t8, 0, 0, 0, 0.2e1 * t14, -0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * t3, -t19 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * t9, -t19 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * pkin(8), -t19 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
