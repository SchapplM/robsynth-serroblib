% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (103->25), mult. (206->53), div. (0->0), fcn. (204->4), ass. (0->19)
t19 = cos(qJ(2));
t25 = 0.2e1 * t19;
t24 = -qJ(3) - pkin(5);
t10 = t24 * t19;
t16 = sin(pkin(6));
t17 = cos(pkin(6));
t18 = sin(qJ(2));
t22 = t24 * t18;
t4 = -t16 * t10 - t17 * t22;
t6 = -t17 * t10 + t16 * t22;
t23 = t4 ^ 2 + t6 ^ 2;
t15 = -t19 * pkin(2) - pkin(1);
t7 = t16 * t18 - t17 * t19;
t8 = t16 * t19 + t17 * t18;
t21 = 0.2e1 * t4 * t8 - 0.2e1 * t6 * t7;
t13 = t17 * pkin(2) + pkin(3);
t11 = t16 * pkin(2) + qJ(4);
t2 = t7 * pkin(3) - t8 * qJ(4) + t15;
t1 = [1, 0, 0, t18 ^ 2, t18 * t25, 0, 0, 0, pkin(1) * t25, -0.2e1 * pkin(1) * t18, t21, t15 ^ 2 + t23, 0.2e1 * t2 * t7, t21, -0.2e1 * t2 * t8, t2 ^ 2 + t23; 0, 0, 0, 0, 0, t18, t19, 0, -t18 * pkin(5), -t19 * pkin(5), (-t16 * t7 - t17 * t8) * pkin(2), (t16 * t6 - t17 * t4) * pkin(2), -t4, -t11 * t7 - t13 * t8, t6, t6 * t11 - t4 * t13; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t16 ^ 2 + t17 ^ 2) * pkin(2) ^ 2, 0.2e1 * t13, 0, 0.2e1 * t11, t11 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t7, 0, -t8, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
