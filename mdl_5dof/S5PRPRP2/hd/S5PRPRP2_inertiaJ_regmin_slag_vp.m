% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:51
% EndTime: 2021-01-15 15:04:52
% DurationCPUTime: 0.20s
% Computational Cost: add. (109->44), mult. (236->73), div. (0->0), fcn. (222->4), ass. (0->30)
t15 = sin(pkin(8));
t30 = -0.2e1 * t15;
t17 = sin(qJ(4));
t29 = t17 * t15;
t16 = cos(pkin(8));
t8 = t17 * t16;
t18 = cos(qJ(4));
t10 = t18 * t15;
t28 = t18 * t16;
t11 = t15 ^ 2;
t12 = t16 ^ 2;
t27 = t11 + t12;
t13 = t17 ^ 2;
t14 = t18 ^ 2;
t26 = t13 + t14;
t25 = qJ(3) * t17;
t24 = qJ(5) * t15;
t23 = t11 * qJ(3);
t22 = qJ(3) * t28;
t7 = -t16 * pkin(3) - t15 * pkin(6) - pkin(2);
t5 = t18 * t7;
t21 = -t18 * t24 + t5;
t1 = (-pkin(4) - t25) * t16 + t21;
t2 = t22 + (t7 - t24) * t17;
t20 = t1 * t18 + t2 * t17;
t9 = t14 * t11;
t6 = (pkin(4) * t17 + qJ(3)) * t15;
t4 = t17 * t7 + t22;
t3 = -t16 * t25 + t5;
t19 = [1, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t11 + t12 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t6 + (-t1 * t17 + t18 * t2) * t15; 0, 1, 0, 0, 0.2e1 * pkin(2) * t16, 0.2e1 * t27 * qJ(3), t27 * qJ(3) ^ 2 + pkin(2) ^ 2, t9, -0.2e1 * t18 * t11 * t17, t28 * t30, 0.2e1 * t15 * t8, t12, -0.2e1 * t3 * t16 + 0.2e1 * t17 * t23, 0.2e1 * t4 * t16 + 0.2e1 * t18 * t23, -0.2e1 * t1 * t16 + 0.2e1 * t6 * t29, 0.2e1 * t6 * t10 + 0.2e1 * t2 * t16, t20 * t30, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t16, 0, -pkin(2), 0, 0, 0, 0, 0, -t28, t8, -t28, t8, -t26 * t15, t20; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t10, -t29, -t10, 0, -pkin(4) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t29, -t16, t3, -t4, (-0.2e1 * pkin(4) - t25) * t16 + t21, -t2, -pkin(4) * t10, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t18, -t17, 0, t18 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t10, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t19;
