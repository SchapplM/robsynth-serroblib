% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t6 = sin(qJ(3));
t20 = -0.2e1 * t6;
t8 = cos(qJ(3));
t19 = 0.2e1 * t8;
t18 = t6 * pkin(5);
t7 = sin(qJ(2));
t17 = t6 * t7;
t16 = t8 * pkin(5);
t15 = t8 * t7;
t9 = cos(qJ(2));
t14 = t9 * t6;
t3 = t6 ^ 2;
t13 = t8 ^ 2 + t3;
t12 = t13 * t7;
t11 = -t6 * pkin(3) + t8 * qJ(4);
t2 = t9 * t8;
t1 = -t8 * pkin(3) - t6 * qJ(4) - pkin(2);
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t7 ^ 2 + t9 ^ 2; 0, 0, t9, -t7, 0, 0, 0, 0, 0, t2, -t14, t2, t12, t14, pkin(5) * t12 - t9 * t1; 0, 1, 0, 0, t3, t6 * t19, 0, 0, 0, pkin(2) * t19, pkin(2) * t20, -0.2e1 * t1 * t8, 0.2e1 * t13 * pkin(5), t1 * t20, t13 * pkin(5) ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t15, -t17, 0, t15, t11 * t7; 0, 0, 0, 0, 0, 0, t6, t8, 0, -t18, -t16, -t18, t11, t16, t11 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
