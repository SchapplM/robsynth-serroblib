% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:24
% DurationCPUTime: 0.34s
% Computational Cost: add. (114->41), mult. (289->70), div. (0->0), fcn. (356->8), ass. (0->32)
t21 = sin(pkin(8));
t20 = sin(pkin(9));
t22 = cos(pkin(9));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t9 = t20 * t27 + t22 * t25;
t3 = t9 * t21;
t38 = t3 ^ 2;
t7 = t20 * t25 - t22 * t27;
t37 = t7 ^ 2;
t24 = sin(qJ(5));
t36 = 0.2e1 * t24;
t35 = t3 * t7;
t34 = t20 * pkin(3);
t33 = t22 * pkin(3);
t16 = t24 ^ 2;
t26 = cos(qJ(5));
t18 = t26 ^ 2;
t32 = t16 + t18;
t31 = t25 ^ 2 + t27 ^ 2;
t12 = pkin(6) + t34;
t30 = t32 * t12;
t23 = cos(pkin(8));
t5 = t7 * t21;
t1 = -t23 * t26 + t24 * t5;
t2 = -t23 * t24 - t26 * t5;
t29 = -t1 * t24 + t2 * t26;
t15 = t23 ^ 2;
t14 = t21 ^ 2;
t13 = -pkin(4) - t33;
t6 = t9 ^ 2;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t14 + t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t15 + t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t9 + t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t9 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t6 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t21, -t27 * t21, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t5, 0, (-t20 * t5 - t22 * t3) * pkin(3), 0, 0, 0, 0, 0, 0, -t3 * t26, t3 * t24, t29, t29 * t12 + t3 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, 0, (t20 * t9 - t22 * t7) * pkin(3), 0, 0, 0, 0, 0, 0, -t7 * t26, t7 * t24, t32 * t9, t7 * t13 + t9 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t34, 0, (t20 ^ 2 + t22 ^ 2) * pkin(3) ^ 2, t16, t26 * t36, 0, t18, 0, 0, -0.2e1 * t13 * t26, t13 * t36, 0.2e1 * t30, t32 * t12 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t26 + t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t9, -t26 * t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t26, 0, -t24 * t12, -t26 * t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t4;
