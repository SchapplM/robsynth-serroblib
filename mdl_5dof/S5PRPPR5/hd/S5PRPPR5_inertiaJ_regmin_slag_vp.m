% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t17 = cos(qJ(5));
t20 = 0.2e1 * t17;
t19 = -pkin(2) - pkin(3);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t7 = t14 * qJ(3) + t13 * t19;
t5 = t13 * qJ(3) - t14 * t19;
t18 = cos(qJ(2));
t16 = sin(qJ(2));
t15 = sin(qJ(5));
t4 = -pkin(6) + t7;
t3 = pkin(4) + t5;
t2 = -t18 * t13 + t16 * t14;
t1 = -t16 * t13 - t18 * t14;
t6 = [1, 0, 0, 0, 0, 0, t16 ^ 2 + t18 ^ 2, 0, 0, t1 ^ 2 + t2 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t18, -t16, t18, t16, t18 * pkin(2) + t16 * qJ(3), -t1, t2, -t1 * t5 + t2 * t7, 0, 0, 0, 0, 0, -t1 * t17, t1 * t15; 0, 1, 0, 0, 0.2e1 * pkin(2), 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0.2e1 * t5, 0.2e1 * t7, t5 ^ 2 + t7 ^ 2, t15 ^ 2, t15 * t20, 0, 0, 0, t3 * t20, -0.2e1 * t3 * t15; 0, 0, 0, 0, 0, 0, -t18, 0, 0, t1 * t14 + t2 * t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, -pkin(2), -t14, t13, t7 * t13 - t5 * t14, 0, 0, 0, 0, 0, -t14 * t17, t15 * t14; 0, 0, 0, 0, 0, 0, 1, 0, 0, t13 ^ 2 + t14 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t2, -t17 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, -t15 * t4, -t17 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t13, -t17 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
