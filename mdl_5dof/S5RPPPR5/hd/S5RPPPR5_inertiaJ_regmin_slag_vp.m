% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t21 = sin(pkin(7));
t20 = sin(pkin(8));
t22 = cos(pkin(8));
t27 = t20 ^ 2 + t22 ^ 2;
t30 = t27 * t21;
t24 = sin(qJ(5));
t25 = cos(qJ(5));
t6 = -t25 * t20 - t24 * t22;
t29 = 0.2e1 * t6;
t23 = cos(pkin(7));
t26 = -pkin(1) - pkin(2);
t11 = t23 * qJ(2) + t21 * t26;
t7 = -qJ(4) + t11;
t28 = pkin(6) - t7;
t9 = t21 * qJ(2) - t23 * t26;
t8 = pkin(3) + t9;
t5 = t24 * t20 - t25 * t22;
t19 = t23 ^ 2;
t17 = t21 ^ 2;
t3 = t22 * pkin(4) + t8;
t2 = t28 * t22;
t1 = t28 * t20;
t4 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t9, 0.2e1 * t11, t11 ^ 2 + t9 ^ 2, 0.2e1 * t8 * t22, -0.2e1 * t8 * t20, -0.2e1 * t27 * t7, t27 * t7 ^ 2 + t8 ^ 2, t6 ^ 2, t5 * t29, 0, 0, 0, -0.2e1 * t3 * t5, t3 * t29; 0, 0, 0, -1, 0, -pkin(1), -t23, t21, t11 * t21 - t9 * t23, -t23 * t22, t23 * t20, -t30, -t8 * t23 + t7 * t30, 0, 0, 0, 0, 0, t23 * t5, -t23 * t6; 0, 0, 0, 0, 0, 1, 0, 0, t17 + t19, 0, 0, 0, t27 * t17 + t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, t8, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t25 * t1 + t24 * t2, -t24 * t1 + t25 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t21, t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
