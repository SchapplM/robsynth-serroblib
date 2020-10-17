% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR10
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
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (128->38), mult. (169->58), div. (0->0), fcn. (177->6), ass. (0->25)
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t25 = -pkin(1) - pkin(2);
t11 = t22 * qJ(2) - t24 * t25;
t10 = -pkin(3) - t11;
t12 = t24 * qJ(2) + t22 * t25;
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t4 = t19 * t10 + t20 * t12;
t21 = sin(qJ(5));
t7 = t19 * t22 - t20 * t24;
t30 = t7 * t21;
t23 = cos(qJ(5));
t29 = t7 * t23;
t3 = t20 * t10 - t19 * t12;
t1 = pkin(4) - t3;
t16 = -t20 * pkin(3) - pkin(4);
t28 = t1 - t16;
t27 = t21 * t23;
t18 = t21 ^ 2;
t15 = t19 * pkin(3) + pkin(7);
t13 = 0.2e1 * t27;
t9 = t19 * t24 + t20 * t22;
t2 = -pkin(7) + t4;
t5 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 1, 0.2e1 * t11, 0.2e1 * t12, t3 ^ 2 + t4 ^ 2, t18, t13, 0, 0, 0, 0.2e1 * t1 * t23, -0.2e1 * t1 * t21; 0, 0, 0, -1, 0, -pkin(1), 0, -t24, t22, -t3 * t7 + t4 * t9, 0, 0, 0, 0, 0, t29, -t30; 0, 0, 0, 0, 0, 1, 0, 0, 0, t7 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -1, -t11, -t12, (t19 * t4 + t20 * t3) * pkin(3), -t18, -0.2e1 * t27, 0, 0, 0, -t28 * t23, t28 * t21; 0, 0, 0, 0, 0, 0, 0, t24, -t22, (t19 * t9 - t20 * t7) * pkin(3), 0, 0, 0, 0, 0, -t29, t30; 0, 0, 0, 0, 0, 0, 1, 0, 0, (t19 ^ 2 + t20 ^ 2) * pkin(3) ^ 2, t18, t13, 0, 0, 0, -0.2e1 * t16 * t23, 0.2e1 * t16 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0, -t21 * t2, -t23 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t9, -t23 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, -t21 * t15, -t23 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
