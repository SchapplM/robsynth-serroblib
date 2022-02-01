% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:46
% DurationCPUTime: 0.19s
% Computational Cost: add. (157->46), mult. (290->78), div. (0->0), fcn. (272->6), ass. (0->34)
t18 = sin(pkin(8));
t35 = -0.2e1 * t18;
t22 = sin(qJ(4));
t19 = sin(pkin(7));
t9 = t19 * pkin(1) + qJ(3);
t34 = t22 * t9;
t14 = t18 ^ 2;
t23 = cos(qJ(4));
t33 = t14 * t23;
t32 = t22 * t18;
t20 = cos(pkin(8));
t10 = t22 * t20;
t12 = t23 * t18;
t31 = t23 * t20;
t15 = t20 ^ 2;
t30 = t14 + t15;
t16 = t22 ^ 2;
t17 = t23 ^ 2;
t29 = t16 + t17;
t28 = qJ(5) * t18;
t27 = t9 * t31;
t21 = cos(pkin(7));
t13 = -t21 * pkin(1) - pkin(2);
t7 = -t20 * pkin(3) - t18 * pkin(6) + t13;
t5 = t23 * t7;
t26 = -t23 * t28 + t5;
t1 = (-pkin(4) - t34) * t20 + t26;
t2 = t27 + (t7 - t28) * t22;
t25 = t1 * t23 + t2 * t22;
t11 = t17 * t14;
t6 = (pkin(4) * t22 + t9) * t18;
t4 = t22 * t7 + t27;
t3 = -t9 * t10 + t5;
t8 = [1, 0, 0, (t19 ^ 2 + t21 ^ 2) * pkin(1) ^ 2, -0.2e1 * t13 * t20, 0.2e1 * t30 * t9, t30 * t9 ^ 2 + t13 ^ 2, t11, -0.2e1 * t22 * t33, t31 * t35, 0.2e1 * t18 * t10, t15, 0.2e1 * t14 * t34 - 0.2e1 * t3 * t20, 0.2e1 * t4 * t20 + 0.2e1 * t9 * t33, -0.2e1 * t1 * t20 + 0.2e1 * t6 * t32, 0.2e1 * t6 * t12 + 0.2e1 * t2 * t20, t25 * t35, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t20 + (-t1 * t22 + t2 * t23) * t18; 0, 0, 0, 1, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t14 + t11 + t15; 0, 0, 0, 0, -t20, 0, t13, 0, 0, 0, 0, 0, -t31, t10, -t31, t10, -t29 * t18, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t32, -t20, t3, -t4, (-0.2e1 * pkin(4) - t34) * t20 + t26, -t2, -pkin(4) * t12, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t12, -t32, -t12, 0, -pkin(4) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t23, -t22, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t12, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
