% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:49
% DurationCPUTime: 0.13s
% Computational Cost: add. (71->23), mult. (127->34), div. (0->0), fcn. (107->6), ass. (0->21)
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t11 = sin(pkin(8));
t23 = pkin(1) * t11;
t12 = cos(pkin(8));
t9 = t12 * pkin(1) + pkin(2);
t21 = t14 * t9 + t16 * t23;
t2 = qJ(4) + t21;
t24 = 0.2e1 * t2;
t18 = 0.2e1 * qJ(4);
t22 = qJ(4) + t2;
t5 = -t14 * t23 + t16 * t9;
t4 = -pkin(3) - t5;
t19 = -0.2e1 * pkin(3);
t17 = -pkin(3) - pkin(7);
t15 = cos(qJ(5));
t13 = sin(qJ(5));
t10 = t15 ^ 2;
t8 = -0.2e1 * t15 * t13;
t1 = -pkin(7) + t4;
t3 = [1, 0, 0, (t11 ^ 2 + t12 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t5, -0.2e1 * t21, 0.2e1 * t4, t24, t2 ^ 2 + t4 ^ 2, t10, t8, 0, 0, 0, t13 * t24, t15 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t5, -t21, t19 - t5, t18 + t21, -t4 * pkin(3) + t2 * qJ(4), t10, t8, 0, 0, 0, t22 * t13, t22 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, t19, t18, pkin(3) ^ 2 + qJ(4) ^ 2, t10, t8, 0, 0, 0, t13 * t18, t15 * t18; 0, 0, 0, 0, 0, 0, 0, 1, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, t15 * t1, -t13 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, t15 * t17, -t13 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
