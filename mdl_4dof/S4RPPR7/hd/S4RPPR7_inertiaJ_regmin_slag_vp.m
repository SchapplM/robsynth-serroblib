% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = sin(pkin(6));
t19 = 0.2e1 * t11 * pkin(3) + (2 * qJ(2));
t18 = 2 * qJ(2);
t12 = cos(pkin(6));
t6 = t11 ^ 2 + t12 ^ 2;
t13 = -pkin(1) - qJ(3);
t17 = -pkin(5) + t13;
t16 = qJ(2) ^ 2;
t15 = cos(qJ(4));
t14 = sin(qJ(4));
t5 = t17 * t12;
t4 = t17 * t11;
t3 = -t14 * t11 + t15 * t12;
t2 = t15 * t11 + t14 * t12;
t1 = t6 * t13;
t7 = [1, 0, 0, -2 * pkin(1), t18, pkin(1) ^ 2 + t16, t11 * t18, t12 * t18, -0.2e1 * t1, t6 * t13 ^ 2 + t16, t3 ^ 2, -0.2e1 * t3 * t2, 0, 0, 0, t2 * t19, t3 * t19; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t6, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t11, t12, 0, qJ(2), 0, 0, 0, 0, 0, t2, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, 0, -t14 * t4 + t15 * t5, -t14 * t5 - t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
