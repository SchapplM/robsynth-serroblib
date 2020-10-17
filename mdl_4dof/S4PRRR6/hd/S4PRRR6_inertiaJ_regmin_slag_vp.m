% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:07
% DurationCPUTime: 0.10s
% Computational Cost: add. (38->19), mult. (92->38), div. (0->0), fcn. (119->6), ass. (0->20)
t15 = cos(qJ(3));
t21 = -0.2e1 * t15 * pkin(3) - (2 * pkin(2));
t20 = 0.2e1 * t15;
t19 = pkin(5) + pkin(6);
t11 = sin(qJ(4));
t18 = t11 * pkin(3);
t14 = cos(qJ(4));
t17 = t14 * pkin(3);
t12 = sin(qJ(3));
t6 = t11 * t15 + t14 * t12;
t5 = t11 * t12 - t14 * t15;
t16 = cos(qJ(2));
t13 = sin(qJ(2));
t8 = t19 * t15;
t7 = t19 * t12;
t4 = t5 * t13;
t3 = t6 * t13;
t2 = t11 * t7 - t14 * t8;
t1 = -t11 * t8 - t14 * t7;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t16, -t13, 0, 0, 0, 0, 0, t16 * t15, -t16 * t12, 0, 0, 0, 0, 0, -t16 * t5, -t16 * t6; 0, 1, 0, 0, t12 ^ 2, t12 * t20, 0, 0, 0, pkin(2) * t20, -0.2e1 * pkin(2) * t12, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t21, t6 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t13, -t15 * t13, 0, 0, 0, 0, 0, -t3, t4; 0, 0, 0, 0, 0, 0, t12, t15, 0, -t12 * pkin(5), -t15 * pkin(5), 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
