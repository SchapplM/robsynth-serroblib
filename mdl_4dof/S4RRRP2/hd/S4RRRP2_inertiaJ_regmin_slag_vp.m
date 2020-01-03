% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = cos(qJ(3));
t24 = 0.2e1 * t16;
t14 = sin(qJ(3));
t23 = t14 * pkin(3);
t22 = sin(qJ(2)) * pkin(1);
t21 = t16 * pkin(6);
t9 = pkin(6) + t22;
t20 = t16 * t9;
t19 = cos(qJ(2)) * pkin(1);
t10 = -pkin(2) - t19;
t18 = pkin(2) - t10;
t11 = -t16 * pkin(3) - pkin(2);
t13 = t14 ^ 2;
t12 = t16 * qJ(4);
t8 = t14 * t24;
t7 = t12 + t21;
t6 = (-qJ(4) - pkin(6)) * t14;
t5 = t11 - t19;
t4 = t7 * t16;
t3 = t12 + t20;
t2 = (-qJ(4) - t9) * t14;
t1 = t3 * t16;
t15 = [1, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t22, t13, t8, 0, 0, 0, -0.2e1 * t10 * t16, 0.2e1 * t10 * t14, -0.2e1 * t2 * t14 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 1, t19, -t22, t13, t8, 0, 0, 0, t18 * t16, -t18 * t14, t1 + t4 + (-t2 - t6) * t14, t5 * t11 + t2 * t6 + t3 * t7; 0, 0, 0, 1, 0, 0, t13, t8, 0, 0, 0, pkin(2) * t24, -0.2e1 * pkin(2) * t14, -0.2e1 * t6 * t14 + 0.2e1 * t4, t11 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * t9, -t20, -t23, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * pkin(6), -t21, -t23, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t15;
