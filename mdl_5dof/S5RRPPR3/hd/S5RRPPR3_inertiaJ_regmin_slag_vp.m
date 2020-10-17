% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (88->29), mult. (152->47), div. (0->0), fcn. (128->6), ass. (0->24)
t18 = sin(qJ(5));
t28 = 0.2e1 * t18;
t20 = cos(qJ(5));
t27 = 0.2e1 * t20;
t26 = sin(qJ(2)) * pkin(1);
t14 = cos(qJ(2)) * pkin(1);
t13 = t14 + pkin(2);
t25 = pkin(2) + t13;
t16 = sin(pkin(8));
t10 = t16 * pkin(2) + qJ(4);
t17 = cos(pkin(8));
t23 = t17 * t26;
t6 = t16 * t13 + t23;
t2 = qJ(4) + t6;
t24 = t10 + t2;
t12 = -t17 * pkin(2) - pkin(3);
t7 = t16 * t26;
t5 = t17 * t13 - t7;
t4 = -pkin(3) - t5;
t15 = t20 ^ 2;
t9 = -0.2e1 * t20 * t18;
t8 = -pkin(7) + t12;
t1 = -pkin(7) + t4;
t3 = [1, 0, 0, 1, 0.2e1 * t14, -0.2e1 * t26, t5 ^ 2 + t6 ^ 2, 0.2e1 * t4, 0.2e1 * t2, t2 ^ 2 + t4 ^ 2, t15, t9, 0, 0, 0, t2 * t28, t2 * t27; 0, 0, 0, 1, t14, -t26, (t16 * t6 + t17 * t5) * pkin(2), -t25 * t17 - 0.2e1 * pkin(3) + t7, t25 * t16 + 0.2e1 * qJ(4) + t23, t2 * t10 + t4 * t12, t15, t9, 0, 0, 0, t24 * t18, t24 * t20; 0, 0, 0, 1, 0, 0, (t16 ^ 2 + t17 ^ 2) * pkin(2) ^ 2, 0.2e1 * t12, 0.2e1 * t10, t10 ^ 2 + t12 ^ 2, t15, t9, 0, 0, 0, t10 * t28, t10 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, t20 * t1, -t18 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, t20 * t8, -t18 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
