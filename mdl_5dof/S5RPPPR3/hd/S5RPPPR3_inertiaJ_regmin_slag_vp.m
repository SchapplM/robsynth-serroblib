% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:01
% EndTime: 2019-12-31 17:44:02
% DurationCPUTime: 0.16s
% Computational Cost: add. (68->26), mult. (123->40), div. (0->0), fcn. (120->6), ass. (0->21)
t18 = sin(pkin(8));
t20 = cos(pkin(8));
t11 = t18 ^ 2 + t20 ^ 2;
t21 = cos(pkin(7));
t15 = -t21 * pkin(1) - pkin(2);
t25 = t18 * qJ(4) - t15;
t28 = 0.2e1 * (pkin(3) + pkin(4)) * t20 + 0.2e1 * t25;
t27 = -0.2e1 * t20;
t19 = sin(pkin(7));
t13 = t19 * pkin(1) + qJ(3);
t26 = t11 * t13 ^ 2;
t23 = cos(qJ(5));
t22 = sin(qJ(5));
t9 = t18 * t13;
t7 = t18 * t23 - t20 * t22;
t6 = t18 * t22 + t20 * t23;
t5 = (-pkin(6) + t13) * t20;
t4 = -t18 * pkin(6) + t9;
t3 = -t20 * pkin(3) - t25;
t1 = 0.2e1 * t11 * t13;
t2 = [1, 0, 0, (t19 ^ 2 + t21 ^ 2) * pkin(1) ^ 2, t15 * t27, 0.2e1 * t15 * t18, t1, t15 ^ 2 + t26, t3 * t27, t1, -0.2e1 * t3 * t18, t3 ^ 2 + t26, t7 ^ 2, -0.2e1 * t7 * t6, 0, 0, 0, t6 * t28, t7 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t11, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t20, t18, 0, t15, -t20, 0, -t18, t3, 0, 0, 0, 0, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t22 * t5 + t23 * t4, -t22 * t4 - t23 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t2;
