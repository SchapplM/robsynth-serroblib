% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:36:14
% EndTime: 2019-05-05 13:36:16
% DurationCPUTime: 0.33s
% Computational Cost: add. (209->44), mult. (374->79), div. (0->0), fcn. (441->8), ass. (0->41)
t30 = sin(pkin(10));
t32 = cos(pkin(10));
t35 = sin(qJ(5));
t44 = cos(qJ(5));
t14 = t44 * t30 + t35 * t32;
t11 = t14 ^ 2;
t13 = t35 * t30 - t44 * t32;
t12 = t13 ^ 2;
t50 = -t11 - t12;
t31 = sin(pkin(9));
t20 = t31 * pkin(1) + qJ(3);
t49 = t20 ^ 2;
t48 = -0.2e1 * t13;
t47 = 0.2e1 * t14;
t46 = 0.2e1 * t20;
t33 = cos(pkin(9));
t23 = -t33 * pkin(1) - pkin(2);
t19 = -qJ(4) + t23;
t45 = -pkin(7) + t19;
t34 = sin(qJ(6));
t43 = t13 * t34;
t36 = cos(qJ(6));
t42 = t13 * t36;
t7 = t34 * t14;
t41 = t34 * t36;
t8 = t36 * t14;
t40 = t30 ^ 2 + t32 ^ 2;
t39 = t13 * t47;
t16 = t30 * pkin(4) + t20;
t38 = pkin(5) * t13 - pkin(8) * t14;
t29 = t36 ^ 2;
t28 = t34 ^ 2;
t10 = t45 * t32;
t9 = t45 * t30;
t6 = t40 * t19;
t5 = t35 * t10 + t44 * t9;
t4 = -t44 * t10 + t35 * t9;
t3 = t14 * pkin(5) + t13 * pkin(8) + t16;
t2 = t34 * t3 + t36 * t5;
t1 = t36 * t3 - t34 * t5;
t15 = [1, 0, 0 (t31 ^ 2 + t33 ^ 2) * pkin(1) ^ 2, 0.2e1 * t23, t46, t23 ^ 2 + t49, t30 * t46, t32 * t46, -0.2e1 * t6, t19 ^ 2 * t40 + t49, t12, t39, 0, 0, 0, t16 * t47, t16 * t48, t29 * t12, -0.2e1 * t12 * t41, t8 * t48, t34 * t39, t11, 0.2e1 * t1 * t14 - 0.2e1 * t4 * t43, -0.2e1 * t2 * t14 - 0.2e1 * t4 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t23, 0, 0, -t40, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t34, t50 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t20, 0, 0, 0, 0, 0, t14, -t13, 0, 0, 0, 0, 0, t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t4, -t5, -t13 * t41 (t28 - t29) * t13, t7, t8, 0, t34 * t38 - t4 * t36, t4 * t34 + t36 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, -t8, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, 0, 0, 0, 0, -t42, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t28, 0.2e1 * t41, 0, 0, 0, 0.2e1 * pkin(5) * t36, -0.2e1 * pkin(5) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t43, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, -t34 * pkin(8), -t36 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
