% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:27:19
% EndTime: 2021-01-15 12:27:22
% DurationCPUTime: 0.25s
% Computational Cost: add. (188->38), mult. (308->59), div. (0->0), fcn. (352->4), ass. (0->30)
t19 = sin(qJ(4));
t20 = sin(qJ(3));
t21 = cos(qJ(4));
t22 = cos(qJ(3));
t9 = -t19 * t20 + t21 * t22;
t7 = t9 ^ 2;
t8 = t19 * t22 + t21 * t20;
t33 = t8 ^ 2 + t7;
t12 = t20 * pkin(3) + qJ(2);
t5 = t8 * pkin(4) + t12;
t31 = 0.2e1 * t5;
t30 = 0.2e1 * t12;
t29 = 0.2e1 * qJ(2);
t28 = t9 * pkin(4);
t27 = t19 * pkin(3);
t23 = -pkin(1) - pkin(6);
t10 = (-pkin(7) + t23) * t20;
t15 = t22 * t23;
t11 = -t22 * pkin(7) + t15;
t3 = -t19 * t10 + t21 * t11;
t1 = -t9 * qJ(5) + t3;
t4 = -t21 * t10 - t19 * t11;
t2 = -t8 * qJ(5) - t4;
t26 = t1 * t9 + t2 * t8;
t18 = t21 * pkin(3);
t14 = t18 + pkin(4);
t25 = t9 * t14 + t8 * t27;
t24 = 0.2e1 * pkin(4);
t16 = -0.2e1 * t27;
t6 = [1, 0, 0, -2 * pkin(1), t29, (pkin(1) ^ 2) + qJ(2) ^ 2, t22 ^ 2, -0.2e1 * t22 * t20, 0, 0, 0, t20 * t29, t22 * t29, t7, -0.2e1 * t8 * t9, 0, 0, 0, t8 * t30, t9 * t30, t8 * t31, t9 * t31, -0.2e1 * t26, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t26; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, t15, -t20 * t23, 0, 0, t9, -t8, 0, t3, t4, t1, -t2, -t25, t1 * t14 + t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, 0, 0, 0, 0, t9, -t8, t9, -t8, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, t16, 0.2e1 * t14, t16, 0, t19 ^ 2 * pkin(3) ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t3, t4, t1, -t2, -t28, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, t9, -t8, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t18, -t27, t24 + t18, -t27, 0, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
