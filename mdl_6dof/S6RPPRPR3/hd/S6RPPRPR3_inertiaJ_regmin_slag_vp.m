% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:50
% EndTime: 2019-05-05 14:09:51
% DurationCPUTime: 0.36s
% Computational Cost: add. (291->56), mult. (488->98), div. (0->0), fcn. (574->8), ass. (0->39)
t30 = sin(pkin(10));
t32 = cos(pkin(10));
t36 = cos(qJ(4));
t50 = sin(qJ(4));
t14 = -t30 * t36 - t32 * t50;
t16 = -t30 * t50 + t32 * t36;
t53 = (t14 * t30 - t16 * t32) * pkin(4);
t12 = t14 ^ 2;
t13 = t16 ^ 2;
t52 = t12 + t13;
t31 = sin(pkin(9));
t18 = t31 * pkin(1) + qJ(3);
t51 = 0.2e1 * t18;
t34 = sin(qJ(6));
t9 = t34 * t14;
t49 = t34 * t16;
t35 = cos(qJ(6));
t48 = t34 * t35;
t47 = t35 * t16;
t17 = t50 * pkin(4) + t18;
t33 = cos(pkin(9));
t25 = -t33 * pkin(1) - pkin(2);
t44 = -pkin(7) + t25;
t39 = t50 * t44;
t11 = -t50 * qJ(5) + t39;
t40 = (-qJ(5) + t44) * t36;
t4 = t30 * t11 - t32 * t40;
t6 = t32 * t11 + t30 * t40;
t43 = t6 * t14 + t4 * t16;
t23 = t30 * pkin(4) + pkin(8);
t24 = -t32 * pkin(4) - pkin(5);
t42 = t14 * t23 + t16 * t24;
t29 = t35 ^ 2;
t28 = t34 ^ 2;
t10 = t35 * t14;
t3 = -t14 * pkin(5) - t16 * pkin(8) + t17;
t2 = t34 * t3 + t35 * t6;
t1 = t35 * t3 - t34 * t6;
t5 = [1, 0, 0 (t31 ^ 2 + t33 ^ 2) * pkin(1) ^ 2, 0.2e1 * t25, t51, t18 ^ 2 + t25 ^ 2, t36 ^ 2, -0.2e1 * t36 * t50, 0, 0, 0, t50 * t51, t36 * t51, 0.2e1 * t43, t17 ^ 2 + t4 ^ 2 + t6 ^ 2, t29 * t13, -0.2e1 * t13 * t48, -0.2e1 * t14 * t47, 0.2e1 * t14 * t49, t12, -0.2e1 * t1 * t14 + 0.2e1 * t4 * t49, 0.2e1 * t2 * t14 + 0.2e1 * t4 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t14 + t6 * t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t25, 0, 0, 0, 0, 0, 0, 0, -t52, -t43, 0, 0, 0, 0, 0, -t52 * t34, -t52 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t50, 0, t36 * t44, -t39, t53 (t30 * t6 - t32 * t4) * pkin(4), t34 * t47 (-t28 + t29) * t16, -t9, -t10, 0, t42 * t34 - t4 * t35, t4 * t34 + t42 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t36, 0 (t14 * t32 + t16 * t30) * pkin(4), 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t50, 0, -t53, 0, 0, 0, 0, 0, t47, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t30 ^ 2 + t32 ^ 2) * pkin(4) ^ 2, t28, 0.2e1 * t48, 0, 0, 0, -0.2e1 * t24 * t35, 0.2e1 * t24 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, -t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t49, -t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t35, 0, -t34 * t23, -t35 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
