% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR3
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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:42:09
% EndTime: 2021-01-15 15:42:11
% DurationCPUTime: 0.22s
% Computational Cost: add. (170->37), mult. (357->74), div. (0->0), fcn. (429->6), ass. (0->29)
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t13 = t22 * t25 - t23 * t27;
t21 = -t27 * pkin(3) - pkin(2);
t34 = 0.2e1 * t13 * pkin(4) + 0.2e1 * t21;
t33 = 0.2e1 * t21;
t32 = 0.2e1 * t27;
t31 = t22 * pkin(3);
t30 = t23 * pkin(3);
t29 = -qJ(4) - pkin(6);
t17 = t29 * t25;
t18 = t29 * t27;
t7 = t23 * t17 + t22 * t18;
t8 = t22 * t17 - t23 * t18;
t26 = cos(qJ(5));
t24 = sin(qJ(5));
t20 = pkin(4) + t30;
t15 = t22 * t27 + t23 * t25;
t11 = -t24 * t20 - t26 * t31;
t10 = t26 * t20 - t24 * t31;
t6 = -t24 * t13 + t26 * t15;
t5 = t26 * t13 + t24 * t15;
t4 = -t13 * pkin(7) + t8;
t3 = -t15 * pkin(7) + t7;
t2 = -t24 * t3 - t26 * t4;
t1 = -t24 * t4 + t26 * t3;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t7 + t15 * t8, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, t25 ^ 2, t25 * t32, 0, 0, 0, pkin(2) * t32, -0.2e1 * pkin(2) * t25, t13 * t33, t15 * t33, -0.2e1 * t8 * t13 - 0.2e1 * t7 * t15, t21 ^ 2 + t7 ^ 2 + t8 ^ 2, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t34, t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, -t13, -t15, 0, (-t13 * t23 + t15 * t22) * pkin(3), 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, t25, t27, 0, -t25 * pkin(6), -t27 * pkin(6), t7, -t8, (-t13 * t22 - t15 * t23) * pkin(3), (t22 * t8 + t23 * t7) * pkin(3), 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t30, -0.2e1 * t31, 0, (t22 ^ 2 + t23 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t10, 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, t21, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
