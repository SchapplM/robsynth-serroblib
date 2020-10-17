% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (171->47), mult. (246->70), div. (0->0), fcn. (212->4), ass. (0->32)
t17 = sin(qJ(4));
t14 = t17 ^ 2;
t19 = cos(qJ(4));
t24 = t19 ^ 2 + t14;
t36 = t24 * pkin(7);
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t21 = -pkin(1) - pkin(2);
t7 = t20 * qJ(2) + t18 * t21;
t5 = -pkin(7) + t7;
t23 = t24 * t5;
t35 = -0.2e1 * t17;
t34 = 0.2e1 * t19;
t6 = t18 * qJ(2) - t20 * t21;
t4 = pkin(3) + t6;
t33 = pkin(3) + t4;
t8 = -t19 * pkin(4) - t17 * qJ(5) - pkin(3);
t1 = t6 - t8;
t32 = -t1 + t8;
t31 = t17 * pkin(7);
t30 = t17 * t5;
t29 = t19 * pkin(7);
t28 = t19 * t5;
t27 = t17 * t18;
t26 = t17 * t19;
t25 = t19 * t18;
t2 = t24 * t18;
t9 = -t17 * pkin(4) + t19 * qJ(5);
t12 = t20 * t19;
t11 = t20 * t17;
t10 = 0.2e1 * t26;
t3 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 1, 0.2e1 * t6, 0.2e1 * t7, t14, t10, 0, 0, 0, t4 * t34, t4 * t35, t1 * t34, -0.2e1 * t23, 0.2e1 * t1 * t17, t24 * t5 ^ 2 + t1 ^ 2; 0, 0, 0, -1, 0, -pkin(1), 0, -t20, t18, 0, 0, 0, 0, 0, -t12, t11, -t12, -t2, -t11, -t1 * t20 + t5 * t2; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t18 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, -1, -t6, -t7, -t14, -0.2e1 * t26, 0, 0, 0, -t33 * t19, t33 * t17, t32 * t19, t23 - t36, t32 * t17, pkin(7) * t23 + t1 * t8; 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, 0, 0, 0, 0, t12, -t11, t12, t2, t11, pkin(7) * t2 - t20 * t8; 0, 0, 0, 0, 0, 0, 1, 0, 0, t14, t10, 0, 0, 0, pkin(3) * t34, pkin(3) * t35, -0.2e1 * t8 * t19, 0.2e1 * t36, t8 * t35, t24 * pkin(7) ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0, -t30, -t28, -t30, -t9, t28, t9 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25, -t27, 0, t25, t9 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, 0, -t31, -t29, -t31, t9, t29, t9 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
