% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.17s
% Computational Cost: add. (147->47), mult. (176->57), div. (0->0), fcn. (158->6), ass. (0->30)
t27 = sin(qJ(5));
t38 = -0.2e1 * t27;
t29 = cos(qJ(5));
t37 = 0.2e1 * t29;
t35 = cos(qJ(2)) * pkin(1);
t20 = pkin(2) + t35;
t14 = -pkin(3) - t20;
t26 = cos(pkin(8));
t12 = t26 * t14;
t23 = sin(qJ(2)) * pkin(1);
t17 = t23 + qJ(3);
t25 = sin(pkin(8));
t3 = t25 * t17 - t12;
t1 = pkin(4) + t3;
t31 = -pkin(2) - pkin(3);
t19 = t26 * t31;
t8 = t25 * qJ(3) - t19;
t6 = pkin(4) + t8;
t36 = t1 + t6;
t34 = t26 * t29;
t5 = t25 * t14 + t26 * t17;
t10 = t26 * qJ(3) + t25 * t31;
t33 = 0.2e1 * pkin(2);
t32 = 0.2e1 * qJ(3);
t24 = t27 ^ 2;
t16 = t26 * t27;
t15 = t27 * t37;
t7 = -pkin(7) + t10;
t2 = -pkin(7) + t5;
t4 = [1, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t23, 0.2e1 * t20, 0.2e1 * t17, t17 ^ 2 + t20 ^ 2, 0.2e1 * t3, 0.2e1 * t5, t3 ^ 2 + t5 ^ 2, t24, t15, 0, 0, 0, t1 * t37, t1 * t38; 0, 0, 0, 1, t35, -t23, t33 + t35, t32 + t23, t20 * pkin(2) + t17 * qJ(3), -t12 - t19 + (qJ(3) + t17) * t25, t10 + t5, t5 * t10 + t3 * t8, t24, t15, 0, 0, 0, t36 * t29, -t36 * t27; 0, 0, 0, 1, 0, 0, t33, t32, pkin(2) ^ 2 + qJ(3) ^ 2, 0.2e1 * t8, 0.2e1 * t10, t10 ^ 2 + t8 ^ 2, t24, t15, 0, 0, 0, t6 * t37, t6 * t38; 0, 0, 0, 0, 0, 0, -1, 0, -t20, -t26, t25, t5 * t25 - t3 * t26, 0, 0, 0, 0, 0, -t34, t16; 0, 0, 0, 0, 0, 0, -1, 0, -pkin(2), -t26, t25, t10 * t25 - t8 * t26, 0, 0, 0, 0, 0, -t34, t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t29, 0, -t27 * t2, -t29 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t29, 0, -t27 * t7, -t29 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t25, -t29 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
