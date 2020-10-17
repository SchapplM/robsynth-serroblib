% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (48->27), mult. (112->59), div. (0->0), fcn. (112->6), ass. (0->24)
t12 = sin(qJ(3));
t24 = 0.2e1 * t12;
t13 = cos(qJ(4));
t23 = pkin(3) * t13;
t11 = sin(qJ(4));
t9 = sin(pkin(7));
t4 = t9 * pkin(1) + pkin(5);
t22 = t11 * t4;
t21 = t11 * t12;
t20 = t11 * t13;
t14 = cos(qJ(3));
t19 = t11 * t14;
t18 = t13 * t12;
t17 = t13 * t14;
t16 = t14 * t24;
t10 = cos(pkin(7));
t5 = -t10 * pkin(1) - pkin(2);
t8 = t13 ^ 2;
t7 = t12 ^ 2;
t6 = t11 ^ 2;
t3 = -t14 * pkin(3) - t12 * pkin(6) + t5;
t2 = t11 * t3 + t4 * t17;
t1 = t13 * t3 - t4 * t19;
t15 = [1, 0, 0, (t10 ^ 2 + t9 ^ 2) * pkin(1) ^ 2, t7, t16, 0, 0, 0, -0.2e1 * t5 * t14, t5 * t24, t8 * t7, -0.2e1 * t7 * t20, -0.2e1 * t12 * t17, t11 * t16, t14 ^ 2, -0.2e1 * t1 * t14 + 0.2e1 * t7 * t22, 0.2e1 * t7 * t4 * t13 + 0.2e1 * t2 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t12, t14, 0, -t12 * t4, -t14 * t4, t11 * t18, (-t6 + t8) * t12, -t19, -t17, 0, -t4 * t18 + (-pkin(3) * t12 + pkin(6) * t14) * t11, pkin(6) * t17 + (t22 - t23) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, 0, 0, 0, 0, 0, t17, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t6, 0.2e1 * t20, 0, 0, 0, 0.2e1 * t23, -0.2e1 * pkin(3) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t21, -t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, 0, -t11 * pkin(6), -t13 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
