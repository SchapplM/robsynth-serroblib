% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:37
% DurationCPUTime: 0.19s
% Computational Cost: add. (51->25), mult. (116->44), div. (0->0), fcn. (83->2), ass. (0->24)
t16 = sin(qJ(2));
t12 = t16 ^ 2;
t17 = cos(qJ(2));
t13 = t17 ^ 2;
t28 = t12 + t13;
t7 = t17 * qJ(3);
t27 = -t16 * pkin(2) + t7;
t26 = -0.2e1 * t16;
t25 = 0.2e1 * t17;
t23 = t28 * pkin(5) ^ 2;
t22 = t16 * t17;
t14 = pkin(2) + qJ(4);
t21 = -t16 * qJ(3) - pkin(1);
t19 = qJ(3) ^ 2;
t18 = 0.2e1 * qJ(3);
t11 = t17 * pkin(5);
t9 = t16 * pkin(5);
t6 = 0.2e1 * t22;
t5 = t17 * pkin(3) + t11;
t4 = t16 * pkin(3) + t9;
t3 = -t17 * pkin(2) + t21;
t2 = 0.2e1 * t28 * pkin(5);
t1 = -t14 * t17 + t21;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t12, t6, 0, t13, 0, 0, pkin(1) * t25, pkin(1) * t26, t2, pkin(1) ^ 2 + t23, 0, 0, 0, t12, t6, t13, t2, t3 * t25, t3 * t26, t3 ^ 2 + t23, 0, 0, 0, t13, -0.2e1 * t22, t12, 0.2e1 * t4 * t16 + 0.2e1 * t5 * t17, t1 * t26, -0.2e1 * t1 * t17, t1 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t17, 0, -t9, -t11, 0, 0, 0, -t16, -t17, 0, 0, 0, t27, t9, t11, t27 * pkin(5), 0, -t17, t16, 0, 0, 0, -t14 * t16 + t7, t5, -t4, t5 * qJ(3) - t4 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t18, pkin(2) ^ 2 + t19, 1, 0, 0, 0, 0, 0, 0, t18, 0.2e1 * t14, t14 ^ 2 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, t9, 0, 0, 0, 0, 0, 0, t16, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
