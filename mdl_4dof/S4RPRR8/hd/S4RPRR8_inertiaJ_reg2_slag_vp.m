% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (121->32), mult. (211->54), div. (0->0), fcn. (221->4), ass. (0->28)
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t4 = t17 * t20 + t19 * t18;
t6 = -t17 * t18 + t19 * t20;
t32 = (t17 * t4 + t19 * t6) * pkin(3);
t3 = t6 ^ 2;
t30 = t4 ^ 2;
t31 = t3 + t30;
t11 = t18 * pkin(3) + qJ(2);
t29 = 0.2e1 * t11;
t28 = 0.2e1 * qJ(2);
t27 = t17 * pkin(3);
t26 = t19 * pkin(3);
t14 = t18 ^ 2;
t15 = t20 ^ 2;
t10 = t14 + t15;
t21 = -pkin(1) - pkin(5);
t8 = (-pkin(6) + t21) * t18;
t13 = t20 * t21;
t9 = -t20 * pkin(6) + t13;
t1 = -t17 * t8 + t19 * t9;
t2 = t17 * t9 + t19 * t8;
t25 = t1 * t6 + t2 * t4;
t22 = qJ(2) ^ 2;
t7 = t10 * t21;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t28, (pkin(1) ^ 2) + t22, t15, -0.2e1 * t20 * t18, 0, t14, 0, 0, t18 * t28, t20 * t28, -0.2e1 * t7, t10 * t21 ^ 2 + t22, t3, -0.2e1 * t6 * t4, 0, t30, 0, 0, t4 * t29, t6 * t29, -0.2e1 * t25, t1 ^ 2 + t11 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t10, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, t13, -t18 * t21, 0, 0, 0, 0, t6, 0, -t4, 0, t1, -t2, -t32, (t1 * t19 + t17 * t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t18, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t4, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t27, 0, (t17 ^ 2 + t19 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t4, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
