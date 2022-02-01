% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:58
% EndTime: 2022-01-23 09:16:59
% DurationCPUTime: 0.27s
% Computational Cost: add. (318->56), mult. (707->105), div. (0->0), fcn. (791->8), ass. (0->46)
t36 = cos(pkin(8));
t51 = -0.2e1 * t36;
t50 = 0.2e1 * t36;
t49 = t36 * pkin(4);
t37 = sin(qJ(5));
t48 = t37 * pkin(4);
t39 = cos(qJ(5));
t47 = t39 * pkin(4);
t33 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t24 = t40 * t33 + t38 * t35;
t34 = sin(pkin(8));
t19 = t24 * t34;
t26 = -pkin(2) * t36 - t34 * qJ(3) - pkin(1);
t22 = t35 * t26;
t44 = t35 * t34;
t11 = -pkin(6) * t44 + t22 + (-qJ(2) * t33 - pkin(3)) * t36;
t43 = qJ(2) * t36;
t17 = t33 * t26 + t35 * t43;
t45 = t33 * t34;
t15 = -pkin(6) * t45 + t17;
t7 = t38 * t11 + t40 * t15;
t5 = -t19 * pkin(7) + t7;
t46 = t39 * t5;
t30 = t34 * qJ(2);
t25 = pkin(3) * t45 + t30;
t31 = t34 ^ 2;
t42 = t31 * qJ(2);
t23 = -t38 * t33 + t40 * t35;
t20 = t23 * t34;
t6 = t40 * t11 - t38 * t15;
t4 = -t20 * pkin(7) - t49 + t6;
t1 = -t37 * t5 + t39 * t4;
t41 = qJ(2) ^ 2;
t32 = t36 ^ 2;
t29 = t31 * t41;
t16 = -t33 * t43 + t22;
t14 = t37 * t23 + t39 * t24;
t13 = t39 * t23 - t37 * t24;
t12 = t19 * pkin(4) + t25;
t9 = -t37 * t19 + t39 * t20;
t8 = t39 * t19 + t37 * t20;
t2 = t37 * t4 + t46;
t3 = [1, 0, 0, pkin(1) * t50, 0.2e1 * (t31 + t32) * qJ(2), pkin(1) ^ 2 + t32 * t41 + t29, -0.2e1 * t16 * t36 + 0.2e1 * t33 * t42, 0.2e1 * t17 * t36 + 0.2e1 * t35 * t42, t16 ^ 2 + t17 ^ 2 + t29, t20 ^ 2, -0.2e1 * t20 * t19, t20 * t51, t19 * t50, t32, 0.2e1 * t25 * t19 - 0.2e1 * t6 * t36, 0.2e1 * t25 * t20 + 0.2e1 * t7 * t36, t9 ^ 2, -0.2e1 * t9 * t8, t9 * t51, t8 * t50, t32, -0.2e1 * t1 * t36 + 0.2e1 * t12 * t8, 0.2e1 * t12 * t9 + 0.2e1 * t2 * t36; 0, 0, 0, -t36, 0, -pkin(1), -t35 * t36, t33 * t36, t16 * t35 + t17 * t33, 0, 0, 0, 0, 0, -t23 * t36, t24 * t36, 0, 0, 0, 0, 0, -t13 * t36, t14 * t36; 0, 0, 0, 0, 0, 1, 0, 0, t33 ^ 2 + t35 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t45, t44, t30, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t36, t6, -t7, 0, 0, t9, -t8, -t36, -t36 * t47 + t1, -t46 + (-t4 + t49) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t47, -0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t36, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t47, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
