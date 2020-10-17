% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:53
% EndTime: 2019-12-31 18:32:54
% DurationCPUTime: 0.29s
% Computational Cost: add. (202->46), mult. (406->89), div. (0->0), fcn. (470->6), ass. (0->40)
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(3));
t41 = cos(qJ(3));
t13 = t26 * t23 - t41 * t24;
t14 = t41 * t23 + t26 * t24;
t18 = -t24 * pkin(2) - pkin(1);
t30 = -t14 * qJ(4) + t18;
t6 = t13 * pkin(3) + t30;
t45 = -0.2e1 * t6;
t44 = 0.2e1 * t18;
t43 = 0.2e1 * qJ(4);
t42 = pkin(3) + pkin(7);
t40 = t14 * t13;
t25 = sin(qJ(5));
t39 = t25 * t13;
t38 = t25 * t14;
t27 = cos(qJ(5));
t37 = t27 * t13;
t10 = t27 * t14;
t36 = t27 * t25;
t35 = pkin(6) + qJ(2);
t34 = t23 ^ 2 + t24 ^ 2;
t33 = qJ(4) * t13;
t32 = 0.2e1 * t40;
t31 = t14 * t42 + t33;
t15 = t35 * t23;
t16 = t35 * t24;
t7 = t41 * t15 + t26 * t16;
t8 = -t26 * t15 + t41 * t16;
t22 = t27 ^ 2;
t21 = t25 ^ 2;
t12 = t14 ^ 2;
t11 = t13 ^ 2;
t5 = -t13 * pkin(4) + t8;
t4 = t14 * pkin(4) + t7;
t3 = t42 * t13 + t30;
t2 = t25 * t4 + t27 * t3;
t1 = -t25 * t3 + t27 * t4;
t9 = [1, 0, 0, 0.2e1 * pkin(1) * t24, -0.2e1 * pkin(1) * t23, 0.2e1 * t34 * qJ(2), t34 * qJ(2) ^ 2 + pkin(1) ^ 2, t12, -0.2e1 * t40, 0, 0, 0, t13 * t44, t14 * t44, -0.2e1 * t8 * t13 + 0.2e1 * t7 * t14, t13 * t45, t14 * t45, t6 ^ 2 + t7 ^ 2 + t8 ^ 2, t21 * t11, 0.2e1 * t11 * t36, t25 * t32, t27 * t32, t12, 0.2e1 * t1 * t14 - 0.2e1 * t5 * t37, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t39; 0, 0, 0, -t24, t23, 0, -pkin(1), 0, 0, 0, 0, 0, t13, t14, 0, -t13, -t14, t6, 0, 0, 0, 0, 0, -t38, -t10; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t7, -t8, -pkin(3) * t14 - t33, t7, t8, -t7 * pkin(3) + t8 * qJ(4), t13 * t36, (-t21 + t22) * t13, t10, -t38, 0, t5 * t25 - t31 * t27, t31 * t25 + t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t43, pkin(3) ^ 2 + qJ(4) ^ 2, t22, -0.2e1 * t36, 0, 0, 0, t25 * t43, t27 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, t7, 0, 0, 0, 0, 0, t10, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t37, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, -t27 * t42, t25 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
