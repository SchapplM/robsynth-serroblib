% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:03
% EndTime: 2021-01-15 11:14:04
% DurationCPUTime: 0.21s
% Computational Cost: add. (218->39), mult. (390->73), div. (0->0), fcn. (422->6), ass. (0->28)
t24 = cos(pkin(7));
t20 = -t24 * pkin(1) - pkin(2);
t36 = cos(qJ(3));
t15 = -t36 * pkin(3) + t20;
t40 = 0.2e1 * t15;
t25 = sin(qJ(3));
t39 = 0.2e1 * t25;
t21 = sin(pkin(8));
t38 = t21 * pkin(3);
t23 = cos(pkin(8));
t37 = t23 * pkin(3);
t22 = sin(pkin(7));
t33 = t22 * pkin(1) + pkin(6);
t28 = t36 * t33;
t10 = t36 * qJ(4) + t28;
t29 = (-qJ(4) - t33) * t25;
t5 = t21 * t10 - t23 * t29;
t7 = t23 * t10 + t21 * t29;
t35 = t5 ^ 2 + t7 ^ 2;
t12 = t21 * t25 - t23 * t36;
t14 = t21 * t36 + t23 * t25;
t34 = t12 ^ 2 + t14 ^ 2;
t32 = t5 * t12 + t7 * t14;
t30 = -0.2e1 * t7 * t12 + 0.2e1 * t5 * t14;
t18 = pkin(4) + t37;
t16 = qJ(5) + t38;
t3 = t12 * pkin(4) - t14 * qJ(5) + t15;
t1 = [1, 0, 0, (t22 ^ 2 + t24 ^ 2) * pkin(1) ^ 2, t25 ^ 2, t36 * t39, 0, 0, 0, -0.2e1 * t20 * t36, t20 * t39, t12 * t40, t14 * t40, t30, t15 ^ 2 + t35, 0.2e1 * t3 * t12, t30, -0.2e1 * t3 * t14, t3 ^ 2 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t32; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, t25, t36, 0, -t25 * t33, -t28, -t5, -t7, (-t12 * t21 - t14 * t23) * pkin(3), (t21 * t7 - t23 * t5) * pkin(3), -t5, -t16 * t12 - t18 * t14, t7, t7 * t16 - t5 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t25, -t12, -t14, 0, (-t12 * t23 + t14 * t21) * pkin(3), -t12, 0, t14, -t12 * t18 + t14 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t37, -0.2e1 * t38, 0, (t21 ^ 2 + t23 ^ 2) * pkin(3) ^ 2, 0.2e1 * t18, 0, 0.2e1 * t16, t16 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14, 0, t15, t12, 0, -t14, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
