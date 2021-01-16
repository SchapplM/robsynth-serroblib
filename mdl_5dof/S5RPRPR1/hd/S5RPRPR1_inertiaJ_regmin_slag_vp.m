% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:52
% EndTime: 2021-01-15 11:33:54
% DurationCPUTime: 0.26s
% Computational Cost: add. (210->38), mult. (361->71), div. (0->0), fcn. (429->6), ass. (0->33)
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t17 = -t26 * t29 + t27 * t31;
t18 = t26 * t31 + t27 * t29;
t43 = (t17 * t27 + t18 * t26) * pkin(3);
t22 = t29 * pkin(3) + qJ(2);
t42 = 0.2e1 * t18 * pkin(4) + 0.2e1 * t22;
t41 = 0.2e1 * t22;
t40 = 0.2e1 * qJ(2);
t39 = t26 * pkin(3);
t38 = t27 * pkin(3);
t37 = t17 ^ 2 + t18 ^ 2;
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t36 = t30 * t17 - t28 * t18;
t32 = -pkin(1) - pkin(6);
t19 = (-qJ(4) + t32) * t29;
t24 = t31 * t32;
t20 = -t31 * qJ(4) + t24;
t9 = -t26 * t19 + t27 * t20;
t10 = t27 * t19 + t26 * t20;
t35 = t10 * t18 + t9 * t17;
t5 = t28 * t17 + t30 * t18;
t23 = pkin(4) + t38;
t13 = -t28 * t23 - t30 * t39;
t12 = t30 * t23 - t28 * t39;
t4 = -t18 * pkin(7) + t10;
t3 = -t17 * pkin(7) + t9;
t2 = -t28 * t3 - t30 * t4;
t1 = -t28 * t4 + t30 * t3;
t6 = [1, 0, 0, -2 * pkin(1), t40, (pkin(1) ^ 2) + qJ(2) ^ 2, t31 ^ 2, -0.2e1 * t31 * t29, 0, 0, 0, t29 * t40, t31 * t40, t18 * t41, t17 * t41, -0.2e1 * t35, t10 ^ 2 + t22 ^ 2 + t9 ^ 2, t36 ^ 2, -0.2e1 * t36 * t5, 0, 0, 0, t5 * t42, t36 * t42; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, t24, -t29 * t32, t9, -t10, -t43, (t10 * t26 + t27 * t9) * pkin(3), 0, 0, t36, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, t17, -t18, 0, t43, 0, 0, 0, 0, 0, t36, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t38, -0.2e1 * t39, 0, (t26 ^ 2 + t27 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t12, 0.2e1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, t22, 0, 0, 0, 0, 0, t5, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
