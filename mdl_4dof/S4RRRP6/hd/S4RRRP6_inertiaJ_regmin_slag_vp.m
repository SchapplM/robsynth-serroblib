% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP6
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
% MM_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (96->41), mult. (215->86), div. (0->0), fcn. (197->4), ass. (0->31)
t15 = sin(qJ(2));
t30 = -0.2e1 * t15;
t17 = cos(qJ(2));
t29 = 0.2e1 * t17;
t16 = cos(qJ(3));
t28 = pkin(2) * t16;
t14 = sin(qJ(3));
t27 = pkin(3) * t14;
t26 = pkin(5) * t14;
t25 = t14 * t16;
t24 = t14 * t17;
t23 = t16 * t15;
t22 = t16 * t17;
t21 = -qJ(4) - pkin(6);
t20 = qJ(4) * t15;
t19 = t15 * t29;
t18 = pkin(5) * t22;
t13 = t16 ^ 2;
t12 = t15 ^ 2;
t11 = t14 ^ 2;
t10 = -t16 * pkin(3) - pkin(2);
t9 = t21 * t16;
t8 = t21 * t14;
t7 = -t17 * pkin(2) - t15 * pkin(6) - pkin(1);
t6 = (pkin(5) + t27) * t15;
t5 = t16 * t7;
t4 = t14 * t7 + t18;
t3 = -pkin(5) * t24 + t5;
t2 = t18 + (t7 - t20) * t14;
t1 = -t16 * t20 + t5 + (-pkin(3) - t26) * t17;
t31 = [1, 0, 0, t12, t19, 0, 0, 0, pkin(1) * t29, pkin(1) * t30, t13 * t12, -0.2e1 * t12 * t25, t22 * t30, t14 * t19, t17 ^ 2, 0.2e1 * t12 * t26 - 0.2e1 * t3 * t17, 0.2e1 * t12 * pkin(5) * t16 + 0.2e1 * t4 * t17, 0.2e1 * (-t1 * t16 - t14 * t2) * t15, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t15, t17, 0, -t15 * pkin(5), -t17 * pkin(5), t14 * t23, (-t11 + t13) * t15, -t24, -t22, 0, -pkin(5) * t23 + (-pkin(2) * t15 + pkin(6) * t17) * t14, pkin(6) * t22 + (t26 - t28) * t15, (-t15 * t8 + t2) * t16 + (t15 * t9 - t1) * t14, t1 * t8 + t6 * t10 - t2 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t11, 0.2e1 * t25, 0, 0, 0, 0.2e1 * t28, -0.2e1 * pkin(2) * t14, -0.2e1 * t8 * t14 - 0.2e1 * t9 * t16, t10 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t14 * t15, -t17, t3, -t4, -pkin(3) * t23, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * pkin(6), -t16 * pkin(6), -t27, t8 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t31;
