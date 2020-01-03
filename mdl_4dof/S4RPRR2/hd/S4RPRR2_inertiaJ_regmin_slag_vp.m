% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR2
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
% MM_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = cos(qJ(4));
t17 = 0.2e1 * t12;
t8 = sin(pkin(7));
t16 = pkin(1) * t8;
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t9 = cos(pkin(7));
t6 = t9 * pkin(1) + pkin(2);
t3 = -t11 * t16 + t13 * t6;
t1 = -pkin(3) - t3;
t15 = pkin(3) - t1;
t4 = -t11 * t6 - t13 * t16;
t10 = sin(qJ(4));
t7 = t10 ^ 2;
t5 = t10 * t17;
t2 = pkin(6) - t4;
t14 = [1, 0, 0, (t8 ^ 2 + t9 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t3, 0.2e1 * t4, t7, t5, 0, 0, 0, -0.2e1 * t1 * t12, 0.2e1 * t1 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t3, t4, t7, t5, 0, 0, 0, t15 * t12, -t15 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, t7, t5, 0, 0, 0, pkin(3) * t17, -0.2e1 * pkin(3) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, -t10 * t2, -t12 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, -t10 * pkin(6), -t12 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t14;
