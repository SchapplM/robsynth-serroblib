% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:51
% DurationCPUTime: 0.20s
% Computational Cost: add. (200->43), mult. (387->94), div. (0->0), fcn. (441->8), ass. (0->35)
t23 = sin(pkin(9));
t25 = cos(pkin(9));
t28 = sin(qJ(3));
t38 = cos(qJ(3));
t12 = t23 * t28 - t25 * t38;
t40 = t12 ^ 2;
t39 = 0.2e1 * t28;
t27 = sin(qJ(5));
t8 = t27 * t12;
t14 = t23 * t38 + t25 * t28;
t37 = t27 * t14;
t29 = cos(qJ(5));
t36 = t27 * t29;
t35 = t29 * t14;
t26 = cos(pkin(8));
t20 = -t26 * pkin(1) - pkin(2);
t24 = sin(pkin(8));
t18 = t24 * pkin(1) + pkin(6);
t34 = t38 * t18;
t33 = (-qJ(4) - t18) * t28;
t17 = t23 * pkin(3) + pkin(7);
t19 = -t25 * pkin(3) - pkin(4);
t32 = -t12 * t17 + t14 * t19;
t15 = -t38 * pkin(3) + t20;
t22 = t29 ^ 2;
t21 = t27 ^ 2;
t11 = t14 ^ 2;
t10 = t38 * qJ(4) + t34;
t9 = t29 * t12;
t6 = t25 * t10 + t23 * t33;
t4 = t23 * t10 - t25 * t33;
t3 = t12 * pkin(4) - t14 * pkin(7) + t15;
t2 = t27 * t3 + t29 * t6;
t1 = -t27 * t6 + t29 * t3;
t5 = [1, 0, 0, (t24 ^ 2 + t26 ^ 2) * pkin(1) ^ 2, t28 ^ 2, t38 * t39, 0, 0, 0, -0.2e1 * t20 * t38, t20 * t39, -0.2e1 * t6 * t12 + 0.2e1 * t4 * t14, t15 ^ 2 + t4 ^ 2 + t6 ^ 2, t22 * t11, -0.2e1 * t11 * t36, 0.2e1 * t12 * t35, -0.2e1 * t12 * t37, t40, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t37, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t12 + t6 * t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t28, t38, 0, -t28 * t18, -t34, (-t12 * t23 - t14 * t25) * pkin(3), (t23 * t6 - t25 * t4) * pkin(3), t27 * t35, (-t21 + t22) * t14, t8, t9, 0, t32 * t27 - t4 * t29, t4 * t27 + t32 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t28, 0, (-t12 * t25 + t14 * t23) * pkin(3), 0, 0, 0, 0, 0, -t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t23 ^ 2 + t25 ^ 2) * pkin(3) ^ 2, t21, 0.2e1 * t36, 0, 0, 0, -0.2e1 * t19 * t29, 0.2e1 * t19 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, t9, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t37, t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, -t27 * t17, -t29 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
