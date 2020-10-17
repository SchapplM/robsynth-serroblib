% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:24
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->31), mult. (273->56), div. (0->0), fcn. (300->6), ass. (0->24)
t18 = cos(pkin(8));
t20 = sin(qJ(4));
t16 = sin(pkin(8));
t26 = cos(qJ(4));
t23 = t26 * t16;
t8 = t20 * t18 + t23;
t29 = -0.2e1 * t8;
t19 = cos(pkin(7));
t13 = -t19 * pkin(1) - pkin(2);
t9 = -t18 * pkin(3) + t13;
t28 = 0.2e1 * t9;
t17 = sin(pkin(7));
t11 = t17 * pkin(1) + qJ(3);
t27 = pkin(6) + t11;
t25 = t20 * t16;
t24 = t16 ^ 2 + t18 ^ 2;
t7 = -t26 * t18 + t25;
t22 = -t7 * pkin(4) + t8 * qJ(5);
t6 = t8 ^ 2;
t5 = t27 * t18;
t3 = -t27 * t25 + t26 * t5;
t2 = t20 * t5 + t27 * t23;
t1 = -t22 + t9;
t4 = [1, 0, 0, (t17 ^ 2 + t19 ^ 2) * pkin(1) ^ 2, -0.2e1 * t13 * t18, 0.2e1 * t13 * t16, 0.2e1 * t24 * t11, t24 * t11 ^ 2 + t13 ^ 2, t6, t7 * t29, 0, 0, 0, t7 * t28, t8 * t28, 0.2e1 * t1 * t7, 0.2e1 * t2 * t8 - 0.2e1 * t3 * t7, t1 * t29, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t7 + t3 * t8; 0, 0, 0, 1, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t6; 0, 0, 0, 0, -t18, t16, 0, t13, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, -t2, -t3, -t2, -pkin(4) * t8 - t7 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
