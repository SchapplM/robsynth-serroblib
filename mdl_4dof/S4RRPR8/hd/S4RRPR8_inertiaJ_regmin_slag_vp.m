% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t18 = cos(qJ(2));
t16 = sin(qJ(2));
t22 = t16 * qJ(3) + pkin(1);
t24 = pkin(2) + pkin(3);
t27 = 0.2e1 * t24 * t18 + 0.2e1 * t22;
t26 = -0.2e1 * t16;
t25 = 0.2e1 * t18;
t13 = t16 ^ 2;
t23 = t18 ^ 2 + t13;
t21 = -t16 * pkin(2) + t18 * qJ(3);
t17 = cos(qJ(4));
t15 = sin(qJ(4));
t12 = t18 * pkin(5);
t11 = t16 * pkin(5);
t10 = -t18 * pkin(6) + t12;
t9 = -t16 * pkin(6) + t11;
t8 = -t18 * pkin(2) - t22;
t7 = t17 * qJ(3) - t15 * t24;
t6 = t15 * qJ(3) + t17 * t24;
t5 = -t18 * t15 + t16 * t17;
t4 = t16 * t15 + t18 * t17;
t2 = t17 * t10 + t15 * t9;
t1 = t15 * t10 - t17 * t9;
t3 = [1, 0, 0, t13, t16 * t25, 0, 0, 0, pkin(1) * t25, pkin(1) * t26, -0.2e1 * t8 * t18, 0.2e1 * t23 * pkin(5), t8 * t26, t23 * pkin(5) ^ 2 + t8 ^ 2, t5 ^ 2, -0.2e1 * t5 * t4, 0, 0, 0, t4 * t27, t5 * t27; 0, 0, 0, 0, 0, t16, t18, 0, -t11, -t12, -t11, t21, t12, t21 * pkin(5), 0, 0, -t5, t4, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t6, 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t17, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
