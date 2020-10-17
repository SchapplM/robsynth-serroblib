% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:38
% DurationCPUTime: 0.30s
% Computational Cost: add. (317->52), mult. (656->108), div. (0->0), fcn. (716->6), ass. (0->36)
t29 = sin(qJ(3));
t30 = sin(qJ(2));
t32 = cos(qJ(3));
t33 = cos(qJ(2));
t14 = t29 * t30 - t32 * t33;
t23 = -pkin(2) * t33 - pkin(1);
t10 = pkin(3) * t14 + t23;
t44 = 0.2e1 * t10;
t43 = 0.2e1 * t23;
t42 = 0.2e1 * t33;
t41 = -pkin(6) - pkin(5);
t28 = sin(qJ(4));
t40 = t28 * pkin(3);
t39 = t29 * pkin(2);
t26 = t30 ^ 2;
t27 = t33 ^ 2;
t38 = t26 + t27;
t31 = cos(qJ(4));
t37 = t31 * t39;
t18 = t41 * t30;
t19 = t41 * t33;
t8 = t32 * t18 + t19 * t29;
t25 = t32 * pkin(2);
t22 = t25 + pkin(3);
t11 = t31 * t22 - t28 * t39;
t9 = t18 * t29 - t19 * t32;
t24 = t31 * pkin(3);
t16 = t29 * t33 + t30 * t32;
t12 = t22 * t28 + t37;
t7 = -t14 * t28 + t16 * t31;
t5 = t31 * t14 + t16 * t28;
t4 = -pkin(7) * t14 + t9;
t3 = -pkin(7) * t16 + t8;
t2 = t28 * t3 + t31 * t4;
t1 = -t28 * t4 + t3 * t31;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t26, t30 * t42, 0, t27, 0, 0, pkin(1) * t42, -0.2e1 * pkin(1) * t30, 0.2e1 * t38 * pkin(5), t38 * pkin(5) ^ 2 + pkin(1) ^ 2, t16 ^ 2, -0.2e1 * t16 * t14, 0, t14 ^ 2, 0, 0, t14 * t43, t16 * t43, -0.2e1 * t14 * t9 - 0.2e1 * t16 * t8, t23 ^ 2 + t8 ^ 2 + t9 ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t44, t7 * t44, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * pkin(5), -t33 * pkin(5), 0, 0, 0, 0, t16, 0, -t14, 0, t8, -t9, (-t14 * t29 - t16 * t32) * pkin(2), (t29 * t9 + t32 * t8) * pkin(2), 0, 0, t7, 0, -t5, 0, t1, -t2, -t11 * t7 - t12 * t5, t1 * t11 + t12 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t39, 0, (t29 ^ 2 + t32 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t12, 0, t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t8, -t9, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, (-t28 * t5 - t31 * t7) * pkin(3), (t1 * t31 + t2 * t28) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t39, 0, 0, 0, 0, 0, 0, 0, 1, t11 + t24, -t37 + (-pkin(3) - t22) * t28, 0, (t11 * t31 + t12 * t28) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, -0.2e1 * t40, 0, (t28 ^ 2 + t31 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, -t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
