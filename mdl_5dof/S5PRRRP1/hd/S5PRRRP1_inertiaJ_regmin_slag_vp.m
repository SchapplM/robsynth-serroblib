% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:40
% EndTime: 2021-01-15 16:14:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (95->41), mult. (169->59), div. (0->0), fcn. (154->4), ass. (0->26)
t13 = sin(qJ(4));
t27 = 0.2e1 * t13;
t15 = cos(qJ(4));
t26 = -0.2e1 * t15;
t25 = 0.2e1 * t15;
t24 = t13 * pkin(4);
t23 = sin(qJ(3)) * pkin(2);
t22 = t15 * pkin(4);
t21 = cos(qJ(3)) * pkin(2);
t10 = -pkin(3) - t21;
t20 = pkin(3) - t10;
t11 = -pkin(3) - t22;
t5 = t11 - t21;
t19 = t11 + t5;
t18 = qJ(5) + pkin(7);
t9 = pkin(7) + t23;
t17 = qJ(5) + t9;
t12 = t13 ^ 2;
t8 = t13 * t25;
t7 = t18 * t15;
t6 = t18 * t13;
t4 = t7 * t15;
t3 = t17 * t15;
t2 = t17 * t13;
t1 = t3 * t15;
t14 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t3 - t15 * t2; 0, 1, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t23, t12, t8, 0, 0, 0, t10 * t26, t10 * t27, t5 * t26, t5 * t27, 0.2e1 * t2 * t13 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t7 - t15 * t6; 0, 0, 0, 0, 1, t21, -t23, t12, t8, 0, 0, 0, t20 * t15, -t20 * t13, -t19 * t15, t19 * t13, t1 + t4 + (t2 + t6) * t13, t5 * t11 + t2 * t6 + t3 * t7; 0, 0, 0, 0, 1, 0, 0, t12, t8, 0, 0, 0, pkin(3) * t25, -0.2e1 * pkin(3) * t13, t11 * t26, t11 * t27, 0.2e1 * t6 * t13 + 0.2e1 * t4, t11 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, t15, -t13, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * t9, -t15 * t9, -t2, -t3, -t24, -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * pkin(7), -t15 * pkin(7), -t6, -t7, -t24, -t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t14;
