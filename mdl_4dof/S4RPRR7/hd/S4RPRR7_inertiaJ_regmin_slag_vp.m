% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (106->29), mult. (245->70), div. (0->0), fcn. (284->6), ass. (0->31)
t19 = cos(pkin(7));
t13 = -t19 * pkin(2) - pkin(1);
t34 = 0.2e1 * t13;
t20 = sin(qJ(4));
t18 = sin(pkin(7));
t21 = sin(qJ(3));
t29 = cos(qJ(3));
t8 = t21 * t18 - t29 * t19;
t33 = t20 * t8;
t22 = cos(qJ(4));
t9 = t29 * t18 + t21 * t19;
t32 = t22 * t9;
t27 = pkin(5) + qJ(2);
t10 = t27 * t18;
t11 = t27 * t19;
t4 = t29 * t10 + t21 * t11;
t31 = t4 * t20;
t30 = t4 * t22;
t28 = t20 * t22;
t26 = t18 ^ 2 + t19 ^ 2;
t25 = -0.2e1 * t9 * t8;
t24 = -pkin(3) * t9 - pkin(6) * t8;
t17 = t22 ^ 2;
t16 = t20 ^ 2;
t7 = t9 ^ 2;
t6 = t22 * t8;
t5 = -t21 * t10 + t29 * t11;
t3 = t8 * pkin(3) - t9 * pkin(6) + t13;
t2 = t20 * t3 + t22 * t5;
t1 = -t20 * t5 + t22 * t3;
t12 = [1, 0, 0, 0.2e1 * pkin(1) * t19, -0.2e1 * pkin(1) * t18, 0.2e1 * t26 * qJ(2), t26 * qJ(2) ^ 2 + pkin(1) ^ 2, t7, t25, 0, 0, 0, t8 * t34, t9 * t34, t17 * t7, -0.2e1 * t7 * t28, 0.2e1 * t8 * t32, t20 * t25, t8 ^ 2, 0.2e1 * t1 * t8 + 0.2e1 * t9 * t31, -0.2e1 * t2 * t8 + 0.2e1 * t9 * t30; 0, 0, 0, -t19, t18, 0, -pkin(1), 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t6, -t33; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, -t4, -t5, t9 * t28, (-t16 + t17) * t9, t33, t6, 0, t24 * t20 - t30, t24 * t22 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t16, 0.2e1 * t28, 0, 0, 0, 0.2e1 * pkin(3) * t22, -0.2e1 * pkin(3) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t20 * t9, t8, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t20 * pkin(6), -t22 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
