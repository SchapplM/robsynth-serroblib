% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t13 = -t22 * t25 + t23 * t27;
t21 = -t27 * pkin(2) - pkin(1);
t32 = -0.2e1 * t13 * pkin(3) + 0.2e1 * t21;
t31 = 0.2e1 * t27;
t30 = pkin(2) * t22;
t29 = -qJ(3) - pkin(5);
t18 = t29 * t25;
t19 = t29 * t27;
t8 = t22 * t18 - t23 * t19;
t7 = t23 * t18 + t22 * t19;
t26 = cos(qJ(4));
t24 = sin(qJ(4));
t20 = t23 * pkin(2) + pkin(3);
t14 = t22 * t27 + t23 * t25;
t11 = -t24 * t20 - t26 * t30;
t10 = t26 * t20 - t24 * t30;
t6 = t24 * t13 + t26 * t14;
t5 = -t26 * t13 + t24 * t14;
t4 = t13 * pkin(6) + t8;
t3 = -t14 * pkin(6) + t7;
t2 = -t24 * t3 - t26 * t4;
t1 = -t24 * t4 + t26 * t3;
t9 = [1, 0, 0, t25 ^ 2, t25 * t31, 0, 0, 0, pkin(1) * t31, -0.2e1 * pkin(1) * t25, 0.2e1 * t8 * t13 - 0.2e1 * t7 * t14, t21 ^ 2 + t7 ^ 2 + t8 ^ 2, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t32, t6 * t32; 0, 0, 0, 0, 0, t25, t27, 0, -t25 * pkin(5), -t27 * pkin(5), (t13 * t22 - t14 * t23) * pkin(2), (t22 * t8 + t23 * t7) * pkin(2), 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t22 ^ 2 + t23 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t10, 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
