% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:39
% EndTime: 2021-01-15 15:22:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (163->36), mult. (330->67), div. (0->0), fcn. (366->4), ass. (0->24)
t32 = cos(qJ(3));
t19 = -t32 * pkin(3) - pkin(2);
t35 = 0.2e1 * t19;
t20 = sin(pkin(8));
t34 = t20 * pkin(3);
t21 = cos(pkin(8));
t33 = t21 * pkin(3);
t28 = t32 * pkin(6);
t14 = t32 * qJ(4) + t28;
t22 = sin(qJ(3));
t25 = (-qJ(4) - pkin(6)) * t22;
t6 = t20 * t14 - t21 * t25;
t8 = t21 * t14 + t20 * t25;
t31 = t6 ^ 2 + t8 ^ 2;
t30 = 0.2e1 * t32;
t10 = t20 * t22 - t21 * t32;
t12 = t20 * t32 + t21 * t22;
t29 = t10 ^ 2 + t12 ^ 2;
t27 = t10 * t6 + t12 * t8;
t24 = -0.2e1 * t8 * t10 + 0.2e1 * t6 * t12;
t17 = pkin(4) + t33;
t15 = qJ(5) + t34;
t3 = t10 * pkin(4) - t12 * qJ(5) + t19;
t1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, t27; 0, 1, 0, 0, t22 ^ 2, t22 * t30, 0, 0, 0, pkin(2) * t30, -0.2e1 * pkin(2) * t22, t10 * t35, t12 * t35, t24, t19 ^ 2 + t31, 0.2e1 * t3 * t10, t24, -0.2e1 * t3 * t12, t3 ^ 2 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t22, -t10, -t12, 0, (-t10 * t21 + t12 * t20) * pkin(3), -t10, 0, t12, -t10 * t17 + t12 * t15; 0, 0, 0, 0, 0, 0, t22, t32, 0, -t22 * pkin(6), -t28, -t6, -t8, (-t10 * t20 - t12 * t21) * pkin(3), (t20 * t8 - t21 * t6) * pkin(3), -t6, -t15 * t10 - t17 * t12, t8, t8 * t15 - t6 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t33, -0.2e1 * t34, 0, (t20 ^ 2 + t21 ^ 2) * pkin(3) ^ 2, 0.2e1 * t17, 0, 0.2e1 * t15, t15 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, t19, t10, 0, -t12, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
