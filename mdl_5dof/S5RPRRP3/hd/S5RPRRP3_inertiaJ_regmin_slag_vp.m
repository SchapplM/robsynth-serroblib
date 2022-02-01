% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP3
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:20
% DurationCPUTime: 0.22s
% Computational Cost: add. (195->43), mult. (354->67), div. (0->0), fcn. (398->6), ass. (0->30)
t21 = cos(pkin(8));
t15 = -t21 * pkin(1) - pkin(2);
t28 = cos(qJ(3));
t13 = -t28 * pkin(3) + t15;
t22 = sin(qJ(4));
t23 = sin(qJ(3));
t27 = cos(qJ(4));
t10 = t22 * t23 - t27 * t28;
t30 = t10 * pkin(4);
t5 = t13 + t30;
t33 = 0.2e1 * t5;
t32 = 0.2e1 * t13;
t31 = 0.2e1 * t23;
t29 = t22 * pkin(3);
t20 = sin(pkin(8));
t14 = t20 * pkin(1) + pkin(6);
t7 = (-pkin(7) - t14) * t23;
t26 = t28 * t14;
t8 = t28 * pkin(7) + t26;
t3 = -t22 * t8 + t27 * t7;
t4 = -t22 * t7 - t27 * t8;
t24 = 0.2e1 * pkin(4);
t19 = t27 * pkin(3);
t18 = -0.2e1 * t29;
t17 = t19 + pkin(4);
t12 = t22 * t28 + t27 * t23;
t9 = t12 ^ 2;
t2 = -t10 * qJ(5) - t4;
t1 = -t12 * qJ(5) + t3;
t6 = [1, 0, 0, (t20 ^ 2 + t21 ^ 2) * pkin(1) ^ 2, t23 ^ 2, t28 * t31, 0, 0, 0, -0.2e1 * t15 * t28, t15 * t31, t9, -0.2e1 * t12 * t10, 0, 0, 0, t10 * t32, t12 * t32, t10 * t33, t12 * t33, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t10 + t2 * t12; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9; 0, 0, 0, 0, 0, 0, t23, t28, 0, -t23 * t14, -t26, 0, 0, t12, -t10, 0, t3, t4, t1, -t2, -t10 * t29 - t17 * t12, t1 * t17 + t2 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t23, 0, 0, 0, 0, 0, -t10, -t12, -t10, -t12, 0, -t10 * t17 + t12 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t19, t18, 0.2e1 * t17, t18, 0, t22 ^ 2 * pkin(3) ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, t3, t4, t1, -t2, -t12 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, -t10, -t12, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, -t29, t24 + t19, -t29, 0, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
