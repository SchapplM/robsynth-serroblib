% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:17
% EndTime: 2019-12-31 18:04:19
% DurationCPUTime: 0.46s
% Computational Cost: add. (304->58), mult. (601->102), div. (0->0), fcn. (690->6), ass. (0->37)
t29 = sin(pkin(8));
t27 = t29 ^ 2;
t30 = cos(pkin(8));
t28 = t30 ^ 2;
t45 = t27 + t28;
t37 = t29 * qJ(3) + pkin(1);
t12 = (pkin(2) + pkin(3)) * t30 + t37;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t14 = t29 * t32 + t30 * t34;
t8 = t14 * pkin(4) + t12;
t44 = 0.2e1 * t8;
t43 = 0.2e1 * t12;
t42 = -0.2e1 * t29;
t31 = sin(qJ(5));
t41 = t31 * pkin(4);
t33 = cos(qJ(5));
t40 = t33 * pkin(4);
t39 = t29 * t30;
t38 = t45 * qJ(2) ^ 2;
t25 = t29 * qJ(2);
t21 = -t29 * pkin(6) + t25;
t22 = (-pkin(6) + qJ(2)) * t30;
t9 = t34 * t21 - t32 * t22;
t10 = t32 * t21 + t34 * t22;
t20 = -t30 * pkin(2) - t37;
t19 = t31 * t34 + t33 * t32;
t18 = -t31 * t32 + t33 * t34;
t17 = 0.2e1 * t45 * qJ(2);
t16 = t29 * t34 - t30 * t32;
t7 = -t31 * t14 + t33 * t16;
t5 = t33 * t14 + t31 * t16;
t4 = -t14 * pkin(7) + t10;
t3 = -t16 * pkin(7) + t9;
t2 = t31 * t3 + t33 * t4;
t1 = t33 * t3 - t31 * t4;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, 0.2e1 * t39, 0, t28, 0, 0, 0.2e1 * pkin(1) * t30, pkin(1) * t42, t17, pkin(1) ^ 2 + t38, t27, 0, -0.2e1 * t39, 0, 0, t28, -0.2e1 * t20 * t30, t17, t20 * t42, t20 ^ 2 + t38, t16 ^ 2, -0.2e1 * t16 * t14, 0, t14 ^ 2, 0, 0, t14 * t43, t16 * t43, -0.2e1 * t10 * t14 - 0.2e1 * t9 * t16, t10 ^ 2 + t12 ^ 2 + t9 ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t44, t7 * t44, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t30, 0, -t29, t20, 0, 0, 0, 0, 0, 0, -t14, -t16, 0, -t12, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t14 - t34 * t16, t10 * t32 + t9 * t34, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t7 - t19 * t5, t1 * t18 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t34 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t9, -t10, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, (-t31 * t5 - t33 * t7) * pkin(4), (t1 * t33 + t2 * t31) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, (t18 * t33 + t19 * t31) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t41, 0, (t31 ^ 2 + t33 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
