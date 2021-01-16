% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:27
% EndTime: 2021-01-15 18:08:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (217->69), mult. (409->121), div. (0->0), fcn. (392->6), ass. (0->42)
t24 = sin(qJ(3));
t42 = 0.2e1 * t24;
t25 = cos(qJ(4));
t41 = pkin(3) * t25;
t23 = sin(qJ(4));
t40 = t23 * pkin(4);
t21 = sin(pkin(8));
t10 = t21 * pkin(1) + pkin(6);
t39 = t10 * t23;
t17 = t23 ^ 2;
t38 = t17 * t24;
t37 = t23 * t24;
t36 = t23 * t25;
t26 = cos(qJ(3));
t35 = t23 * t26;
t14 = t25 * t24;
t15 = t25 * t26;
t34 = t26 * t10;
t33 = -qJ(5) - pkin(7);
t32 = qJ(5) * t24;
t31 = t26 * t42;
t30 = t25 * t34;
t22 = cos(pkin(8));
t11 = -t22 * pkin(1) - pkin(2);
t7 = -t26 * pkin(3) - t24 * pkin(7) + t11;
t5 = t25 * t7;
t29 = -t25 * t32 + t5;
t8 = t33 * t23;
t9 = t33 * t25;
t28 = -t8 * t23 - t9 * t25;
t20 = t26 ^ 2;
t19 = t25 ^ 2;
t18 = t24 ^ 2;
t16 = -t25 * pkin(4) - pkin(3);
t13 = t19 * t24;
t12 = t19 * t18;
t6 = (t10 + t40) * t24;
t4 = t23 * t7 + t30;
t3 = -t23 * t34 + t5;
t2 = t30 + (t7 - t32) * t23;
t1 = (-pkin(4) - t39) * t26 + t29;
t27 = [1, 0, 0, (t21 ^ 2 + t22 ^ 2) * pkin(1) ^ 2, t18, t31, 0, 0, 0, -0.2e1 * t11 * t26, t11 * t42, t12, -0.2e1 * t18 * t36, -0.2e1 * t24 * t15, t23 * t31, t20, 0.2e1 * t18 * t39 - 0.2e1 * t3 * t26, 0.2e1 * t18 * t10 * t25 + 0.2e1 * t4 * t26, -0.2e1 * t1 * t26 + 0.2e1 * t37 * t6, 0.2e1 * t14 * t6 + 0.2e1 * t2 * t26, (-t1 * t25 - t2 * t23) * t42, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t26 + (-t1 * t23 + t2 * t25) * t24; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t18 + t12 + t20; 0, 0, 0, 0, 0, 0, t24, t26, 0, -t24 * t10, -t34, t23 * t14, t13 - t38, -t35, -t15, 0, -t10 * t14 + (-pkin(3) * t24 + pkin(7) * t26) * t23, pkin(7) * t15 + (t39 - t41) * t24, t16 * t37 - t6 * t25 - t8 * t26, t14 * t16 + t6 * t23 - t9 * t26, (-t24 * t8 + t2) * t25 + (t24 * t9 - t1) * t23, t1 * t8 + t6 * t16 - t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, 0, 0, 0, 0, t15, -t35, t15, -t35, t13 + t38, -t26 * t16 + t24 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t17, 0.2e1 * t36, 0, 0, 0, 0.2e1 * t41, -0.2e1 * pkin(3) * t23, -0.2e1 * t16 * t25, 0.2e1 * t16 * t23, 0.2e1 * t28, t16 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t37, -t26, t3, -t4, (-0.2e1 * pkin(4) - t39) * t26 + t29, -t2, -pkin(4) * t14, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t14, -t37, -t14, 0, -pkin(4) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25, 0, -t23 * pkin(7), -t25 * pkin(7), t8, t9, -t40, t8 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t14, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t23, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t27;
