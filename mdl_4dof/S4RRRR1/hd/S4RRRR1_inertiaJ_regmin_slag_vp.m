% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (60->25), mult. (130->36), div. (0->0), fcn. (122->6), ass. (0->25)
t16 = sin(qJ(4));
t28 = pkin(3) * t16;
t17 = sin(qJ(3));
t27 = t17 * pkin(2);
t26 = sin(qJ(2)) * pkin(1);
t19 = cos(qJ(4));
t13 = cos(qJ(2)) * pkin(1);
t11 = t13 + pkin(2);
t20 = cos(qJ(3));
t22 = -t20 * t11 + t17 * t26;
t2 = -pkin(3) + t22;
t25 = t2 * t19;
t12 = t20 * pkin(2);
t10 = -t12 - pkin(3);
t24 = t10 * t19;
t23 = t20 * t26;
t5 = -t17 * t11 - t23;
t15 = t16 ^ 2;
t14 = pkin(3) * t19;
t9 = pkin(7) + t27;
t8 = 0.2e1 * t16 * t19;
t7 = t10 * t16;
t3 = pkin(7) - t5;
t1 = t2 * t16;
t4 = [1, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t26, 1, -0.2e1 * t22, 0.2e1 * t5, t15, t8, 0, 0, 0, -0.2e1 * t25, 0.2e1 * t1; 0, 0, 0, 1, t13, -t26, 1, t12 - t22, -t23 + (-pkin(2) - t11) * t17, t15, t8, 0, 0, 0, (-t10 - t2) * t19, t7 + t1; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t12, -0.2e1 * t27, t15, t8, 0, 0, 0, -0.2e1 * t24, 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 1, -t22, t5, t15, t8, 0, 0, 0, t14 - t25, t1 - t28; 0, 0, 0, 0, 0, 0, 1, t12, -t27, t15, t8, 0, 0, 0, t14 - t24, t7 - t28; 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, t8, 0, 0, 0, 0.2e1 * t14, -0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * t3, -t19 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * t9, -t19 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t19, 0, -t16 * pkin(7), -t19 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
