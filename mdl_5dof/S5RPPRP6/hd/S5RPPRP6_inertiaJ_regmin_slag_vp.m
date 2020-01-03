% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t17 = sin(pkin(7));
t20 = sin(qJ(4));
t18 = cos(pkin(7));
t27 = cos(qJ(4));
t24 = t27 * t18;
t6 = t20 * t17 - t24;
t31 = t6 ^ 2;
t26 = t20 * t18;
t7 = t27 * t17 + t26;
t30 = 0.2e1 * t7;
t29 = 2 * qJ(2);
t19 = -pkin(1) - qJ(3);
t28 = -pkin(6) + t19;
t10 = t17 ^ 2 + t18 ^ 2;
t11 = t17 * pkin(3) + qJ(2);
t25 = t7 ^ 2 + t31;
t9 = t28 * t17;
t2 = t20 * t9 - t28 * t24;
t3 = t28 * t26 + t27 * t9;
t23 = t2 * t6 + t3 * t7;
t22 = t6 * pkin(4) - t7 * qJ(5);
t21 = qJ(2) ^ 2;
t5 = t10 * t19;
t1 = t7 * pkin(4) + t6 * qJ(5) + t11;
t4 = [1, 0, 0, -2 * pkin(1), t29, pkin(1) ^ 2 + t21, t17 * t29, t18 * t29, -0.2e1 * t5, t10 * t19 ^ 2 + t21, t31, t6 * t30, 0, 0, 0, t11 * t30, -0.2e1 * t11 * t6, t1 * t30, -0.2e1 * t23, 0.2e1 * t1 * t6, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t10, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, t23; 0, 0, 0, 0, 0, 1, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, t17, t18, 0, qJ(2), 0, 0, 0, 0, 0, t7, -t6, t7, 0, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, -t2, -t3, -t2, t22, t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, -t6, 0, t7, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
