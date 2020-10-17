% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:23
% EndTime: 2019-05-05 13:32:24
% DurationCPUTime: 0.26s
% Computational Cost: add. (88->38), mult. (150->64), div. (0->0), fcn. (149->6), ass. (0->32)
t17 = sin(pkin(9));
t8 = t17 * pkin(1) + qJ(3);
t34 = t8 ^ 2;
t18 = cos(pkin(9));
t10 = -t18 * pkin(1) - pkin(2);
t6 = qJ(4) - t10;
t33 = 0.2e1 * t6;
t21 = cos(qJ(6));
t32 = 0.2e1 * t21;
t22 = cos(qJ(5));
t16 = t22 ^ 2;
t5 = -pkin(7) + t8;
t31 = t16 * t5;
t19 = sin(qJ(6));
t20 = sin(qJ(5));
t11 = t19 * t20;
t30 = t19 * t21;
t29 = t19 * t22;
t28 = t21 * t20;
t12 = t21 * t22;
t27 = t22 * t20;
t14 = t20 ^ 2;
t26 = -t14 - t16;
t25 = -0.2e1 * t27;
t24 = -pkin(5) * t22 - pkin(8) * t20;
t15 = t21 ^ 2;
t13 = t19 ^ 2;
t4 = 0.2e1 * t8;
t3 = t20 * pkin(5) - t22 * pkin(8) + t6;
t2 = t19 * t3 + t5 * t28;
t1 = -t5 * t11 + t21 * t3;
t7 = [1, 0, 0 (t17 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, 0.2e1 * t10, t4, t10 ^ 2 + t34, t4, t33, t6 ^ 2 + t34, t16, t25, 0, 0, 0, t20 * t33, t22 * t33, t15 * t16, -0.2e1 * t16 * t30, t27 * t32, t19 * t25, t14, 0.2e1 * t1 * t20 - 0.2e1 * t19 * t31, -0.2e1 * t2 * t20 - 0.2e1 * t21 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t10, 0, -1, -t6, 0, 0, 0, 0, 0, -t20, -t22, 0, 0, 0, 0, 0, -t28, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t19, t26 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, t22 * t5, -t20 * t5, t19 * t12 (-t13 + t15) * t22, t11, t28, 0, t5 * t12 + t24 * t19, t24 * t21 - t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22, 0, 0, 0, 0, 0, -t28, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, 0, 0, 0, 0, t12, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t13, 0.2e1 * t30, 0, 0, 0, pkin(5) * t32, -0.2e1 * pkin(5) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t29, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t21, 0, -t19 * pkin(8), -t21 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t7;
