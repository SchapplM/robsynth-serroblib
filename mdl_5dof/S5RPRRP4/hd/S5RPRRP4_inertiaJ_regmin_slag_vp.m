% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP4
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
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:25
% DurationCPUTime: 0.28s
% Computational Cost: add. (260->56), mult. (560->110), div. (0->0), fcn. (589->6), ass. (0->40)
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t25 = sin(pkin(8));
t30 = cos(qJ(3));
t37 = t30 * t25;
t28 = sin(qJ(3));
t39 = t28 * t25;
t12 = -t27 * t39 + t29 * t37;
t44 = -0.2e1 * t12;
t43 = -0.2e1 * t25;
t26 = cos(pkin(8));
t42 = 0.2e1 * t26;
t41 = t27 * pkin(3);
t22 = t29 * pkin(3);
t17 = -t26 * pkin(2) - t25 * pkin(6) - pkin(1);
t36 = t30 * t26;
t32 = qJ(2) * t36;
t8 = t32 + (-pkin(7) * t25 + t17) * t28;
t40 = t29 * t8;
t38 = t28 * t26;
t14 = pkin(3) * t39 + t25 * qJ(2);
t23 = t25 ^ 2;
t24 = t26 ^ 2;
t35 = t23 + t24;
t34 = qJ(2) * t28;
t33 = t23 * qJ(2);
t13 = t30 * t17;
t6 = -pkin(7) * t37 + t13 + (-pkin(3) - t34) * t26;
t3 = -t27 * t8 + t29 * t6;
t4 = t27 * t6 + t40;
t16 = t27 * t30 + t29 * t28;
t20 = t22 + pkin(4);
t15 = -t27 * t28 + t29 * t30;
t11 = t16 * t25;
t10 = t28 * t17 + t32;
t9 = -t26 * t34 + t13;
t7 = t11 * pkin(4) + t14;
t2 = -t11 * qJ(5) + t4;
t1 = -t26 * pkin(4) - t12 * qJ(5) + t3;
t5 = [1, 0, 0, pkin(1) * t42, pkin(1) * t43, 0.2e1 * t35 * qJ(2), t35 * qJ(2) ^ 2 + pkin(1) ^ 2, t30 ^ 2 * t23, -0.2e1 * t30 * t23 * t28, t36 * t43, 0.2e1 * t25 * t38, t24, -0.2e1 * t9 * t26 + 0.2e1 * t28 * t33, 0.2e1 * t10 * t26 + 0.2e1 * t30 * t33, t12 ^ 2, t11 * t44, t26 * t44, t11 * t42, t24, 0.2e1 * t14 * t11 - 0.2e1 * t3 * t26, 0.2e1 * t14 * t12 + 0.2e1 * t4 * t26, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t11, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, -t26, t25, 0, -pkin(1), 0, 0, 0, 0, 0, -t36, t38, 0, 0, 0, 0, 0, -t15 * t26, t16 * t26, -t16 * t11 - t15 * t12, t1 * t15 + t2 * t16; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t39, -t26, t9, -t10, 0, 0, t12, -t11, -t26, -t26 * t22 + t3, -t40 + (t26 * pkin(3) - t6) * t27, -t11 * t41 - t20 * t12, t1 * t20 + t2 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, 0, 0, 0, 0, t15, -t16, 0, t15 * t20 + t16 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t22, -0.2e1 * t41, 0, t27 ^ 2 * pkin(3) ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t26, t3, -t4, -pkin(4) * t12, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, 0, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t22, -t41, 0, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
