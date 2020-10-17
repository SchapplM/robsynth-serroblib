% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:39:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (149->26), mult. (263->51), div. (0->0), fcn. (324->6), ass. (0->30)
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t14 = t30 * t24 + t28 * t25;
t19 = t24 * pkin(3) + qJ(2);
t40 = 0.2e1 * t14 * pkin(4) + 0.2e1 * t19;
t39 = 0.2e1 * t19;
t38 = 0.2e1 * qJ(2);
t27 = sin(qJ(5));
t37 = t27 * pkin(4);
t29 = cos(qJ(5));
t36 = t29 * pkin(4);
t26 = -pkin(1) - qJ(3);
t35 = -pkin(6) + t26;
t18 = t24 ^ 2 + t25 ^ 2;
t15 = -t28 * t24 + t30 * t25;
t34 = -t27 * t14 + t29 * t15;
t16 = t35 * t24;
t17 = t35 * t25;
t33 = -t28 * t16 + t30 * t17;
t6 = t29 * t14 + t27 * t15;
t32 = -t30 * t16 - t28 * t17;
t31 = qJ(2) ^ 2;
t13 = t18 * t26;
t4 = -t14 * pkin(7) - t32;
t3 = -t15 * pkin(7) + t33;
t2 = -t27 * t3 - t29 * t4;
t1 = -t27 * t4 + t29 * t3;
t5 = [1, 0, 0, -2 * pkin(1), t38, (pkin(1) ^ 2) + t31, t24 * t38, t25 * t38, -0.2e1 * t13, t18 * t26 ^ 2 + t31, t15 ^ 2, -0.2e1 * t15 * t14, 0, 0, 0, t14 * t39, t15 * t39, t34 ^ 2, -0.2e1 * t34 * t6, 0, 0, 0, t6 * t40, t34 * t40; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t18, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t24, t25, 0, qJ(2), 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t6, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t33, t32, 0, 0, t34, -t6, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, 0, 0, 0, 0, t34, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t6, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
