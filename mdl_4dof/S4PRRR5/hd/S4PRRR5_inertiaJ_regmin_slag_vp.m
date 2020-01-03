% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR5
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
% MM_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = cos(qJ(4));
t19 = 0.2e1 * t12;
t13 = cos(qJ(3));
t16 = t13 * pkin(2);
t7 = -pkin(3) - t16;
t18 = pkin(3) - t7;
t10 = sin(qJ(3));
t17 = t10 * pkin(2);
t11 = sin(qJ(2));
t14 = cos(qJ(2));
t2 = t10 * t11 - t13 * t14;
t15 = t2 * t12;
t9 = sin(qJ(4));
t8 = t9 ^ 2;
t6 = pkin(6) + t17;
t4 = t9 * t19;
t3 = t10 * t14 + t13 * t11;
t1 = t2 * t9;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t14, -t11, 0, -t2, -t3, 0, 0, 0, 0, 0, -t15, t1; 0, 1, 0, 0, 1, 0.2e1 * t16, -0.2e1 * t17, t8, t4, 0, 0, 0, -0.2e1 * t7 * t12, 0.2e1 * t7 * t9; 0, 0, 0, 0, 0, -t2, -t3, 0, 0, 0, 0, 0, -t15, t1; 0, 0, 0, 0, 1, t16, -t17, t8, t4, 0, 0, 0, t18 * t12, -t18 * t9; 0, 0, 0, 0, 1, 0, 0, t8, t4, 0, 0, 0, pkin(3) * t19, -0.2e1 * pkin(3) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t3, -t12 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t12, 0, -t9 * t6, -t12 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t12, 0, -t9 * pkin(6), -t12 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
