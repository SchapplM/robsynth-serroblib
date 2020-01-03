% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t8 = sin(qJ(2));
t21 = -0.2e1 * t8;
t9 = cos(qJ(2));
t20 = 0.2e1 * t9;
t19 = t8 * pkin(5);
t10 = pkin(2) + pkin(3);
t6 = t8 ^ 2;
t18 = t9 ^ 2 + t6;
t17 = t9 * qJ(3);
t16 = t8 * qJ(3) + pkin(1);
t15 = -t8 * pkin(2) + t17;
t13 = qJ(3) ^ 2;
t12 = 0.2e1 * qJ(3);
t5 = t9 * pkin(5);
t4 = -t9 * qJ(4) + t5;
t3 = (pkin(5) - qJ(4)) * t8;
t2 = -t9 * pkin(2) - t16;
t1 = t10 * t9 + t16;
t7 = [1, 0, 0, t6, t8 * t20, 0, 0, 0, pkin(1) * t20, pkin(1) * t21, -0.2e1 * t2 * t9, 0.2e1 * t18 * pkin(5), t2 * t21, t18 * pkin(5) ^ 2 + t2 ^ 2, t1 * t20, 0.2e1 * t1 * t8, -0.2e1 * t3 * t8 - 0.2e1 * t4 * t9, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t8, t9, 0, -t19, -t5, -t19, t15, t5, t15 * pkin(5), -t3, t4, t10 * t8 - t17, t4 * qJ(3) - t3 * t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t12, pkin(2) ^ 2 + t13, 0.2e1 * t10, t12, 0, t10 ^ 2 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, t19, 0, 0, -t8, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), -1, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
