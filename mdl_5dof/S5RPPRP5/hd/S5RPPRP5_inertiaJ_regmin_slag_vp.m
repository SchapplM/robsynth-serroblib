% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.20s
% Computational Cost: add. (124->40), mult. (237->61), div. (0->0), fcn. (247->4), ass. (0->22)
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t28 = t17 ^ 2 + t18 ^ 2;
t19 = sin(qJ(4));
t20 = cos(qJ(4));
t6 = t17 * t19 + t18 * t20;
t27 = 0.2e1 * t6;
t7 = t17 * t20 - t18 * t19;
t26 = -0.2e1 * t7;
t25 = -0.2e1 * t17;
t24 = t28 * qJ(2) ^ 2;
t13 = t17 * qJ(2);
t23 = -t17 * pkin(6) + t13;
t22 = t17 * qJ(3) + pkin(1);
t4 = (pkin(2) + pkin(3)) * t18 + t22;
t10 = (-pkin(6) + qJ(2)) * t18;
t9 = -t18 * pkin(2) - t22;
t8 = 0.2e1 * t28 * qJ(2);
t3 = t20 * t10 + t19 * t23;
t2 = t19 * t10 - t20 * t23;
t1 = t6 * pkin(4) - t7 * qJ(5) + t4;
t5 = [1, 0, 0, 0.2e1 * pkin(1) * t18, pkin(1) * t25, t8, pkin(1) ^ 2 + t24, -0.2e1 * t9 * t18, t8, t9 * t25, t9 ^ 2 + t24, t7 ^ 2, t6 * t26, 0, 0, 0, t4 * t27, 0.2e1 * t4 * t7, t1 * t27, 0.2e1 * t2 * t7 - 0.2e1 * t3 * t6, t1 * t26, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -t18, t17, 0, -pkin(1), -t18, 0, -t17, t9, 0, 0, 0, 0, 0, -t6, -t7, -t6, 0, t7, -t1; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t6 - t20 * t7, 0, t3 * t19 - t2 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -t7 * pkin(4) - t6 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t20, 0, t19, t20 * pkin(4) + t19 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
