% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:01
% DurationCPUTime: 0.28s
% Computational Cost: add. (354->45), mult. (682->88), div. (0->0), fcn. (790->6), ass. (0->38)
t33 = cos(qJ(2));
t25 = -t33 * pkin(2) - pkin(1);
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t32 = cos(qJ(3));
t34 = t29 * t30 - t32 * t33;
t13 = pkin(3) * t34 + t25;
t41 = 0.2e1 * t13;
t40 = 0.2e1 * t25;
t39 = 0.2e1 * t33;
t38 = pkin(6) + pkin(7);
t28 = sin(qJ(4));
t37 = t28 * pkin(3);
t36 = t29 * pkin(2);
t31 = cos(qJ(4));
t35 = t31 * t36;
t20 = t38 * t30;
t21 = t38 * t33;
t11 = -t32 * t20 - t29 * t21;
t18 = t29 * t33 + t32 * t30;
t7 = -t18 * pkin(8) + t11;
t12 = t29 * t20 - t32 * t21;
t8 = -pkin(8) * t34 - t12;
t3 = -t28 * t8 + t31 * t7;
t27 = t32 * pkin(2);
t24 = t27 + pkin(3);
t15 = t31 * t24 - t28 * t36;
t4 = -t28 * t7 - t31 * t8;
t26 = t31 * pkin(3);
t23 = t26 + pkin(4);
t16 = t28 * t24 + t35;
t14 = pkin(4) + t15;
t10 = t31 * t18 - t28 * t34;
t9 = t28 * t18 + t31 * t34;
t5 = t9 * pkin(4) + t13;
t2 = -t9 * qJ(5) - t4;
t1 = -t10 * qJ(5) + t3;
t6 = [1, 0, 0, t30 ^ 2, t30 * t39, 0, 0, 0, pkin(1) * t39, -0.2e1 * pkin(1) * t30, t18 ^ 2, -0.2e1 * t18 * t34, 0, 0, 0, t34 * t40, t18 * t40, t10 ^ 2, -0.2e1 * t10 * t9, 0, 0, 0, t9 * t41, t10 * t41, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t30, t33, 0, -t30 * pkin(6), -t33 * pkin(6), 0, 0, t18, -t34, 0, t11, t12, 0, 0, t10, -t9, 0, t3, t4, -t14 * t10 - t16 * t9, t1 * t14 + t2 * t16; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, -0.2e1 * t36, 0, 0, 0, 0, 1, 0.2e1 * t15, -0.2e1 * t16, 0, t14 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t34, 0, t11, t12, 0, 0, t10, -t9, 0, t3, t4, -t23 * t10 - t37 * t9, t1 * t23 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t36, 0, 0, 0, 0, 1, t15 + t26, -t35 + (-pkin(3) - t24) * t28, 0, t14 * t23 + t16 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t37, 0, t28 ^ 2 * pkin(3) ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t3, t4, -pkin(4) * t10, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t15, -t16, 0, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t37, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
