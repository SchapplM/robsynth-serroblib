% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:01
% EndTime: 2021-01-15 11:24:03
% DurationCPUTime: 0.23s
% Computational Cost: add. (209->40), mult. (334->62), div. (0->0), fcn. (366->4), ass. (0->27)
t21 = sin(pkin(7));
t22 = cos(pkin(7));
t23 = sin(qJ(3));
t24 = cos(qJ(3));
t11 = -t21 * t23 + t22 * t24;
t12 = t21 * t24 + t22 * t23;
t42 = (t11 * t22 + t12 * t21) * pkin(3);
t32 = t11 ^ 2 + t12 ^ 2;
t38 = -pkin(1) - pkin(6);
t31 = -qJ(4) + t38;
t14 = t31 * t23;
t29 = t31 * t24;
t6 = t21 * t14 - t22 * t29;
t8 = t22 * t14 + t21 * t29;
t30 = -t6 * t11 + t8 * t12;
t17 = t23 * pkin(3) + qJ(2);
t40 = 0.2e1 * t17;
t39 = 0.2e1 * qJ(2);
t37 = t21 * pkin(3);
t36 = t22 * pkin(3);
t33 = t6 ^ 2 + t8 ^ 2;
t15 = qJ(5) + t37;
t18 = pkin(4) + t36;
t28 = t18 * t11 + t12 * t15;
t26 = 0.2e1 * t30;
t4 = t12 * pkin(4) - t11 * qJ(5) + t17;
t1 = [1, 0, 0, -2 * pkin(1), t39, (pkin(1) ^ 2) + qJ(2) ^ 2, t24 ^ 2, -0.2e1 * t24 * t23, 0, 0, 0, t23 * t39, t24 * t39, t12 * t40, t11 * t40, -t26, t17 ^ 2 + t33, 0.2e1 * t4 * t12, -t26, -0.2e1 * t4 * t11, t4 ^ 2 + t33; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t30, 0, -t32, 0, t30; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t24 * t38, -t23 * t38, -t6, -t8, -t42, (t21 * t8 - t22 * t6) * pkin(3), -t6, -t28, t8, t8 * t15 - t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, t11, -t12, 0, t42, t11, 0, t12, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t36, -0.2e1 * t37, 0, (t21 ^ 2 + t22 ^ 2) * pkin(3) ^ 2, 0.2e1 * t18, 0, 0.2e1 * t15, t15 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t17, t12, 0, -t11, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
