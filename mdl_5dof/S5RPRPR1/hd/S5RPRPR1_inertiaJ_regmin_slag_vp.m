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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:47:42
% EndTime: 2019-12-05 17:47:43
% DurationCPUTime: 0.23s
% Computational Cost: add. (192->37), mult. (327->66), div. (0->0), fcn. (389->6), ass. (0->31)
t31 = sin(pkin(8));
t32 = cos(pkin(8));
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t19 = -t31 * t34 + t32 * t36;
t20 = -t31 * t36 - t32 * t34;
t46 = (t19 * t32 - t20 * t31) * pkin(3);
t28 = t34 * pkin(3) + qJ(2);
t45 = -0.2e1 * t20 * pkin(4) + 0.2e1 * t28;
t44 = 0.2e1 * qJ(2);
t43 = pkin(3) * t31;
t37 = -pkin(1) - pkin(6);
t22 = (-qJ(4) + t37) * t34;
t29 = t36 * t37;
t23 = -t36 * qJ(4) + t29;
t10 = t32 * t22 + t31 * t23;
t42 = t19 ^ 2 + t20 ^ 2;
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t41 = t35 * t19 + t33 * t20;
t9 = -t31 * t22 + t32 * t23;
t40 = t10 * t20 - t9 * t19;
t5 = t33 * t19 - t35 * t20;
t27 = t32 * pkin(3) + pkin(4);
t13 = -t33 * t27 - t35 * t43;
t12 = t35 * t27 - t33 * t43;
t4 = t20 * pkin(7) + t10;
t3 = -t19 * pkin(7) + t9;
t2 = -t33 * t3 - t35 * t4;
t1 = t35 * t3 - t33 * t4;
t6 = [1, 0, 0, -2 * pkin(1), t44, (pkin(1) ^ 2) + qJ(2) ^ 2, t36 ^ 2, -0.2e1 * t36 * t34, 0, 0, 0, t34 * t44, t36 * t44, 0.2e1 * t40, t10 ^ 2 + t28 ^ 2 + t9 ^ 2, t41 ^ 2, -0.2e1 * t41 * t5, 0, 0, 0, t5 * t45, t41 * t45; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t42, -t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, t29, -t34 * t37, -t46, (t10 * t31 + t32 * t9) * pkin(3), 0, 0, t41, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, t46, 0, 0, 0, 0, 0, t41, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t31 ^ 2 + t32 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t12, 0.2e1 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t5, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
