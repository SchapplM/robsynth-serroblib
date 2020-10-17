% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:46:23
% EndTime: 2019-05-05 13:46:24
% DurationCPUTime: 0.33s
% Computational Cost: add. (151->54), mult. (214->89), div. (0->0), fcn. (225->6), ass. (0->38)
t19 = sin(pkin(9));
t20 = cos(pkin(9));
t25 = -pkin(1) - pkin(2);
t12 = t20 * qJ(2) + t19 * t25;
t7 = -qJ(4) + t12;
t40 = -0.2e1 * t7;
t23 = cos(qJ(6));
t39 = 0.2e1 * t23;
t24 = cos(qJ(5));
t18 = t24 ^ 2;
t21 = sin(qJ(6));
t38 = t18 * t21;
t37 = t18 * t23;
t22 = sin(qJ(5));
t36 = t21 * t22;
t35 = t21 * t23;
t34 = t21 * t24;
t33 = t22 * t20;
t32 = t23 * t22;
t31 = t23 * t24;
t30 = t24 * t20;
t29 = t24 * t22;
t16 = t22 ^ 2;
t28 = t16 + t18;
t27 = -0.2e1 * t29;
t10 = t19 * qJ(2) - t20 * t25;
t9 = pkin(3) + t10;
t26 = pkin(5) * t24 + pkin(8) * t22;
t17 = t23 ^ 2;
t15 = t21 ^ 2;
t13 = t19 ^ 2 + t20 ^ 2;
t6 = pkin(7) + t9;
t5 = -t21 * t19 + t20 * t32;
t4 = t23 * t19 + t21 * t33;
t3 = -t22 * pkin(5) + t24 * pkin(8) + t7;
t2 = t21 * t3 + t6 * t32;
t1 = t23 * t3 - t6 * t36;
t8 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t10, 0.2e1 * t12, t10 ^ 2 + t12 ^ 2, -0.2e1 * t9, t40, t7 ^ 2 + t9 ^ 2, t18, t27, 0, 0, 0, t22 * t40, t24 * t40, t17 * t18, -0.2e1 * t18 * t35, t29 * t39, t21 * t27, t16, -0.2e1 * t1 * t22 + 0.2e1 * t6 * t38, 0.2e1 * t2 * t22 + 0.2e1 * t6 * t37; 0, 0, 0, -1, 0, -pkin(1), -t20, t19, -t10 * t20 + t12 * t19, t20, -t19, t7 * t19 - t9 * t20, 0, 0, 0, 0, 0, -t19 * t22, -t24 * t19, 0, 0, 0, 0, 0, -t20 * t38 - t4 * t22, -t20 * t37 - t5 * t22; 0, 0, 0, 0, 0, 1, 0, 0, t13, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t21, t28 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t22, 0, t24 * t6, -t22 * t6, -t21 * t31 (t15 - t17) * t24, -t36, -t32, 0, t26 * t21 + t6 * t31, t26 * t23 - t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t33, 0, 0, 0, 0, 0, -t23 * t30, t21 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t24, 0, 0, 0, 0, 0, -t32, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0, 0, 0, 0, t31, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, 0.2e1 * t35, 0, 0, 0, pkin(5) * t39, -0.2e1 * pkin(5) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t34, -t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, -t21 * pkin(8), -t23 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
