% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:50
% EndTime: 2019-12-05 16:15:51
% DurationCPUTime: 0.18s
% Computational Cost: add. (90->30), mult. (191->56), div. (0->0), fcn. (195->6), ass. (0->27)
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t29 = t22 ^ 2 + t23 ^ 2;
t30 = t29 * qJ(4);
t15 = -t23 * pkin(4) - pkin(3);
t34 = cos(qJ(3)) * pkin(2);
t7 = t15 - t34;
t37 = 0.2e1 * t7;
t36 = 0.2e1 * t15;
t35 = sin(qJ(3)) * pkin(2);
t16 = -pkin(3) - t34;
t33 = pkin(3) - t16;
t32 = t15 + t7;
t14 = qJ(4) + t35;
t31 = t29 * t14;
t26 = cos(qJ(5));
t24 = sin(qJ(5));
t19 = t23 * pkin(7);
t9 = t23 * qJ(4) + t19;
t8 = (-pkin(7) - qJ(4)) * t22;
t6 = t26 * t22 + t24 * t23;
t5 = t24 * t22 - t26 * t23;
t4 = t6 ^ 2;
t3 = t23 * t14 + t19;
t2 = (-pkin(7) - t14) * t22;
t1 = -0.2e1 * t6 * t5;
t10 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35, -0.2e1 * t16 * t23, 0.2e1 * t16 * t22, 0.2e1 * t31, t29 * t14 ^ 2 + t16 ^ 2, t4, t1, 0, 0, 0, t5 * t37, t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t34, -t35, t33 * t23, -t33 * t22, t30 + t31, -t16 * pkin(3) + t14 * t30, t4, t1, 0, 0, 0, t32 * t5, t32 * t6; 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t23, -0.2e1 * pkin(3) * t22, 0.2e1 * t30, t29 * qJ(4) ^ 2 + pkin(3) ^ 2, t4, t1, 0, 0, 0, t5 * t36, t6 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, t16, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, -t23, t22, 0, -pkin(3), 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t26 * t2 - t24 * t3, -t24 * t2 - t26 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, -t24 * t9 + t26 * t8, -t24 * t8 - t26 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
