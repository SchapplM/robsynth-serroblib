% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:36
% EndTime: 2019-05-05 15:50:37
% DurationCPUTime: 0.48s
% Computational Cost: add. (218->64), mult. (390->90), div. (0->0), fcn. (477->6), ass. (0->48)
t32 = sin(qJ(5));
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t45 = cos(qJ(5));
t16 = t32 * t35 + t45 * t33;
t13 = t16 ^ 2;
t15 = t32 * t33 - t45 * t35;
t14 = t15 ^ 2;
t52 = -t13 - t14;
t51 = -0.2e1 * t15;
t50 = 0.2e1 * t16;
t29 = (pkin(1) + qJ(3));
t49 = 2 * t29;
t48 = t32 * pkin(4);
t34 = cos(qJ(6));
t28 = -pkin(7) + qJ(2);
t18 = (-pkin(8) + t28) * t33;
t21 = t35 * t28;
t40 = -t35 * pkin(8) + t21;
t6 = t32 * t18 - t45 * t40;
t47 = t6 * t34;
t41 = t45 * pkin(4);
t24 = -t41 - pkin(5);
t46 = pkin(5) - t24;
t31 = sin(qJ(6));
t11 = t15 * t31;
t44 = t15 * t34;
t9 = t31 * t16;
t43 = t31 * t34;
t10 = t34 * t16;
t42 = t15 * t50;
t19 = t33 * pkin(4) + t29;
t39 = pkin(5) * t15 - pkin(9) * t16;
t23 = pkin(9) + t48;
t38 = -t15 * t24 - t16 * t23;
t37 = (qJ(2) ^ 2);
t36 = 2 * qJ(2);
t27 = t34 ^ 2;
t26 = t31 ^ 2;
t20 = 0.2e1 * t43;
t8 = t15 * t43;
t7 = t45 * t18 + t32 * t40;
t5 = (t26 - t27) * t15;
t4 = t6 * t31;
t3 = t16 * pkin(5) + t15 * pkin(9) + t19;
t2 = t31 * t3 + t34 * t7;
t1 = t34 * t3 - t31 * t7;
t12 = [1, 0, 0, -2 * pkin(1), t36, pkin(1) ^ 2 + t37, t36, t49, t29 ^ 2 + t37, t35 ^ 2, -0.2e1 * t35 * t33, 0, 0, 0, t33 * t49, t35 * t49, t14, t42, 0, 0, 0, t19 * t50, t19 * t51, t27 * t14, -0.2e1 * t14 * t43, t10 * t51, t31 * t42, t13, 0.2e1 * t1 * t16 - 0.2e1 * t6 * t11, -0.2e1 * t2 * t16 - 0.2e1 * t6 * t44; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t29, 0, 0, 0, 0, 0, -t33, -t35, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, -t10, t9; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t31, t52 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, t21, -t33 * t28, 0, 0, -t15, -t16, 0, -t6, -t7, -t8, t5, t9, t10, 0, t31 * t38 - t47, t34 * t38 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t44, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t48, t26, t20, 0, 0, 0, -0.2e1 * t24 * t34, 0.2e1 * t24 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, -t6, -t7, -t8, t5, t9, t10, 0, t31 * t39 - t47, t34 * t39 + t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t44, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t41, -t48, t26, t20, 0, 0, 0, t46 * t34, -t46 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t26, t20, 0, 0, 0, 0.2e1 * pkin(5) * t34, -0.2e1 * pkin(5) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t11, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34, 0, -t31 * t23, -t34 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34, 0, -t31 * pkin(9), -t34 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
