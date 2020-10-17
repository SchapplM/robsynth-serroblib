% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (61->18), mult. (107->36), div. (0->0), fcn. (105->6), ass. (0->20)
t16 = sin(pkin(7));
t10 = t16 * pkin(1) + qJ(3);
t15 = sin(pkin(8));
t26 = 0.2e1 * t15 * pkin(4) + 0.2e1 * t10;
t25 = t10 ^ 2;
t24 = 0.2e1 * t10;
t18 = cos(pkin(7));
t12 = -t18 * pkin(1) - pkin(2);
t9 = -qJ(4) + t12;
t23 = -pkin(6) + t9;
t17 = cos(pkin(8));
t22 = t15 ^ 2 + t17 ^ 2;
t20 = cos(qJ(5));
t19 = sin(qJ(5));
t5 = -t19 * t15 + t20 * t17;
t4 = t20 * t15 + t19 * t17;
t3 = t23 * t17;
t2 = t23 * t15;
t1 = t22 * t9;
t6 = [1, 0, 0, (t16 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, 0.2e1 * t12, t24, t12 ^ 2 + t25, t15 * t24, t17 * t24, -0.2e1 * t1, t22 * t9 ^ 2 + t25, t5 ^ 2, -0.2e1 * t5 * t4, 0, 0, 0, t4 * t26, t5 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t12, 0, 0, -t22, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, t10, 0, 0, 0, 0, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, 0, -t19 * t2 + t20 * t3, -t19 * t3 - t20 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
