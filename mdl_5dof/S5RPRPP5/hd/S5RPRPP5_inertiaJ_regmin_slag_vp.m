% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:41
% EndTime: 2019-12-31 18:16:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (90->32), mult. (121->47), div. (0->0), fcn. (96->2), ass. (0->23)
t16 = cos(qJ(3));
t26 = -0.2e1 * t16;
t25 = 2 * qJ(2);
t17 = pkin(3) + pkin(4);
t19 = -pkin(1) - pkin(6);
t24 = t16 * t19;
t13 = t16 ^ 2;
t15 = sin(qJ(3));
t10 = t15 ^ 2 + t13;
t23 = t15 * qJ(4);
t22 = t16 * qJ(4) - qJ(2);
t9 = t16 * pkin(3) + t23;
t21 = qJ(4) ^ 2;
t20 = 0.2e1 * qJ(4);
t11 = t15 * t19;
t8 = t15 * pkin(3) - t22;
t6 = (-qJ(5) - t19) * t16;
t5 = t17 * t16 + t23;
t4 = t15 * qJ(5) + t11;
t3 = t10 * t19;
t2 = -t17 * t15 + t22;
t1 = t4 * t15 - t6 * t16;
t7 = [1, 0, 0, -2 * pkin(1), t25, pkin(1) ^ 2 + qJ(2) ^ 2, t13, t15 * t26, 0, 0, 0, t15 * t25, t16 * t25, 0.2e1 * t8 * t15, -0.2e1 * t3, t8 * t26, t10 * t19 ^ 2 + t8 ^ 2, -0.2e1 * t2 * t15, 0.2e1 * t2 * t16, 0.2e1 * t1, t2 ^ 2 + t4 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t3, 0, 0, t10, t1; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t24, -t11, t24, -t9, t11, t9 * t19, -t6, t4, t5, t4 * qJ(4) - t6 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t16, 0, t15, t9, t16, t15, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t20, pkin(3) ^ 2 + t21, 0.2e1 * t17, t20, 0, t17 ^ 2 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t24, 0, 0, -t16, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
