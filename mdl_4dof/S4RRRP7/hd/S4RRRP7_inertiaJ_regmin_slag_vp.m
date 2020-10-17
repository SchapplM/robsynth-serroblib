% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:04
% EndTime: 2019-12-31 17:21:05
% DurationCPUTime: 0.22s
% Computational Cost: add. (133->47), mult. (293->101), div. (0->0), fcn. (261->4), ass. (0->37)
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t22 = -t18 * pkin(3) - t16 * qJ(4);
t8 = -pkin(2) + t22;
t39 = -0.2e1 * t8;
t17 = sin(qJ(2));
t38 = -0.2e1 * t17;
t19 = cos(qJ(2));
t37 = 0.2e1 * t19;
t36 = pkin(2) * t16;
t35 = pkin(2) * t18;
t34 = pkin(5) * t16;
t33 = pkin(5) * t18;
t32 = t16 * pkin(6);
t31 = t18 * pkin(6);
t27 = t18 * t19;
t9 = -t19 * pkin(2) - t17 * pkin(6) - pkin(1);
t4 = pkin(5) * t27 + t16 * t9;
t30 = t16 * t18;
t29 = t16 * t19;
t28 = t17 * t16;
t12 = t18 * t17;
t13 = t16 ^ 2;
t15 = t18 ^ 2;
t26 = t13 + t15;
t25 = t19 * qJ(4);
t24 = t17 * t37;
t1 = -t25 + t4;
t7 = t18 * t9;
t2 = -t7 + (pkin(3) + t34) * t19;
t23 = t1 * t18 + t2 * t16;
t21 = -pkin(3) * t16 + t18 * qJ(4);
t14 = t17 ^ 2;
t10 = pkin(6) * t29;
t5 = (pkin(5) - t21) * t17;
t3 = -pkin(5) * t29 + t7;
t6 = [1, 0, 0, t14, t24, 0, 0, 0, pkin(1) * t37, pkin(1) * t38, t15 * t14, -0.2e1 * t14 * t30, t27 * t38, t16 * t24, t19 ^ 2, 0.2e1 * t14 * t34 - 0.2e1 * t3 * t19, 0.2e1 * t14 * t33 + 0.2e1 * t4 * t19, 0.2e1 * t2 * t19 + 0.2e1 * t5 * t28, 0.2e1 * (-t1 * t16 + t18 * t2) * t17, -0.2e1 * t1 * t19 - 0.2e1 * t5 * t12, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t17, t19, 0, -t17 * pkin(5), -t19 * pkin(5), t16 * t12, (-t13 + t15) * t17, -t29, -t27, 0, t10 + (-t33 - t36) * t17, pkin(6) * t27 + (t34 - t35) * t17, -t5 * t18 + t8 * t28 + t10, t23, -t5 * t16 + (-pkin(6) * t19 - t17 * t8) * t18, t23 * pkin(6) + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t13, 0.2e1 * t30, 0, 0, 0, 0.2e1 * t35, -0.2e1 * t36, t18 * t39, 0.2e1 * t26 * pkin(6), t16 * t39, t26 * pkin(6) ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t28, -t19, t3, -t4, t7 + (-0.2e1 * pkin(3) - t34) * t19, t22 * t17, -0.2e1 * t25 + t4, -t2 * pkin(3) + t1 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t32, -t31, -t32, t21, t31, t21 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t12, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
