% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t9 = sin(qJ(4));
t14 = 0.2e1 * t9;
t12 = cos(qJ(2));
t11 = cos(qJ(4));
t10 = sin(qJ(2));
t8 = cos(pkin(7));
t7 = sin(pkin(7));
t6 = -t8 * pkin(2) - pkin(3);
t5 = t7 * pkin(2) + pkin(5);
t3 = t8 * t10 + t7 * t12;
t1 = t7 * t10 - t8 * t12;
t2 = [1, 0, 0, 0, t1 ^ 2 + t3 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t12, -t10, (-t1 * t8 + t3 * t7) * pkin(2), 0, 0, 0, 0, 0, -t1 * t11, t1 * t9; 0, 1, 0, 0, (t7 ^ 2 + t8 ^ 2) * pkin(2) ^ 2, t9 ^ 2, t11 * t14, 0, 0, 0, -0.2e1 * t6 * t11, t6 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t3, -t11 * t3; 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, -t9 * t5, -t11 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t2;
