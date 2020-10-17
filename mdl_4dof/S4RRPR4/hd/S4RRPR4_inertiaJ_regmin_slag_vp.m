% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:33
% DurationCPUTime: 0.16s
% Computational Cost: add. (87->28), mult. (185->56), div. (0->0), fcn. (185->6), ass. (0->27)
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t29 = t22 ^ 2 + t23 ^ 2;
t30 = t29 * qJ(3);
t15 = -t23 * pkin(3) - pkin(2);
t34 = cos(qJ(2)) * pkin(1);
t7 = t15 - t34;
t37 = 0.2e1 * t7;
t36 = 0.2e1 * t15;
t35 = sin(qJ(2)) * pkin(1);
t16 = -pkin(2) - t34;
t33 = pkin(2) - t16;
t32 = t15 + t7;
t14 = qJ(3) + t35;
t31 = t29 * t14;
t26 = cos(qJ(4));
t24 = sin(qJ(4));
t19 = t23 * pkin(6);
t9 = t23 * qJ(3) + t19;
t8 = (-pkin(6) - qJ(3)) * t22;
t6 = t26 * t22 + t24 * t23;
t5 = t24 * t22 - t26 * t23;
t4 = t6 ^ 2;
t3 = t23 * t14 + t19;
t2 = (-pkin(6) - t14) * t22;
t1 = -0.2e1 * t6 * t5;
t10 = [1, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35, -0.2e1 * t16 * t23, 0.2e1 * t16 * t22, 0.2e1 * t31, t29 * t14 ^ 2 + t16 ^ 2, t4, t1, 0, 0, 0, t5 * t37, t6 * t37; 0, 0, 0, 1, t34, -t35, t33 * t23, -t33 * t22, t30 + t31, -t16 * pkin(2) + t14 * t30, t4, t1, 0, 0, 0, t32 * t5, t32 * t6; 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t23, -0.2e1 * pkin(2) * t22, 0.2e1 * t30, t29 * qJ(3) ^ 2 + pkin(2) ^ 2, t4, t1, 0, 0, 0, t5 * t36, t6 * t36; 0, 0, 0, 0, 0, 0, -t23, t22, 0, t16, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, -t23, t22, 0, -pkin(2), 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t26 * t2 - t24 * t3, -t24 * t2 - t26 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, -t24 * t9 + t26 * t8, -t24 * t8 - t26 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
