% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (39->15), mult. (82->35), div. (0->0), fcn. (84->6), ass. (0->16)
t14 = cos(pkin(7));
t15 = cos(pkin(6));
t9 = -t15 * pkin(1) - pkin(2);
t21 = -0.2e1 * t14 * pkin(3) + 0.2e1 * t9;
t13 = sin(pkin(6));
t7 = t13 * pkin(1) + qJ(3);
t20 = pkin(5) + t7;
t12 = sin(pkin(7));
t19 = t12 ^ 2 + t14 ^ 2;
t17 = cos(qJ(4));
t16 = sin(qJ(4));
t4 = t17 * t12 + t16 * t14;
t3 = t16 * t12 - t17 * t14;
t2 = t20 * t14;
t1 = t20 * t12;
t5 = [1, 0, 0, (t13 ^ 2 + t15 ^ 2) * pkin(1) ^ 2, -0.2e1 * t9 * t14, 0.2e1 * t9 * t12, 0.2e1 * t19 * t7, t19 * t7 ^ 2 + t9 ^ 2, t4 ^ 2, -0.2e1 * t4 * t3, 0, 0, 0, t3 * t21, t4 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t14, t12, 0, t9, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, -t17 * t1 - t16 * t2, t16 * t1 - t17 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
