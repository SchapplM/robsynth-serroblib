% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:00
% EndTime: 2021-01-15 15:14:02
% DurationCPUTime: 0.14s
% Computational Cost: add. (85->33), mult. (158->50), div. (0->0), fcn. (181->6), ass. (0->27)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t5 = t14 * t17 - t15 * t19;
t30 = t5 ^ 2;
t16 = sin(qJ(4));
t29 = 0.2e1 * t16;
t18 = cos(qJ(4));
t28 = -0.2e1 * t18;
t7 = t14 * t19 + t15 * t17;
t27 = t16 * t7;
t26 = t18 * pkin(4);
t25 = t18 * t7;
t24 = t5 * t18;
t12 = t16 ^ 2;
t23 = t18 ^ 2 + t12;
t10 = t14 * pkin(2) + pkin(6);
t22 = qJ(5) + t10;
t11 = -t15 * pkin(2) - pkin(3);
t2 = t22 * t16;
t3 = t22 * t18;
t21 = t2 * t16 + t3 * t18;
t8 = t11 - t26;
t4 = t7 ^ 2;
t1 = t5 * t16;
t6 = [1, 0, 0, 0, t4 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t4 + t30; 0, 0, t19, -t17, (t14 * t7 - t15 * t5) * pkin(2), 0, 0, 0, 0, 0, -t24, t1, -t24, t1, t23 * t7, t21 * t7 + t5 * t8; 0, 1, 0, 0, (t14 ^ 2 + t15 ^ 2) * pkin(2) ^ 2, t12, t18 * t29, 0, 0, 0, t11 * t28, t11 * t29, t8 * t28, t8 * t29, 0.2e1 * t21, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t16 - t2 * t18; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25, -t27, -t25, 0, -pkin(4) * t27; 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t16 * t10, -t18 * t10, -t2, -t3, -t16 * pkin(4), -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t16, t18, -t16, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t16, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
