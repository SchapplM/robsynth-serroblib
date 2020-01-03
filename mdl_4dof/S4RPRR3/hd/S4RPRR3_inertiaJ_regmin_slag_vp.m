% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR3
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
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = cos(qJ(3));
t12 = cos(pkin(7));
t9 = -t12 * pkin(1) - pkin(2);
t22 = -0.2e1 * t16 * pkin(3) + 0.2e1 * t9;
t14 = sin(qJ(3));
t21 = 0.2e1 * t14;
t11 = sin(pkin(7));
t8 = t11 * pkin(1) + pkin(5);
t20 = pkin(6) + t8;
t13 = sin(qJ(4));
t19 = t13 * pkin(3);
t15 = cos(qJ(4));
t18 = t15 * pkin(3);
t6 = t13 * t16 + t15 * t14;
t5 = t13 * t14 - t15 * t16;
t4 = t20 * t16;
t3 = t20 * t14;
t2 = t13 * t3 - t15 * t4;
t1 = -t13 * t4 - t15 * t3;
t7 = [1, 0, 0, (t11 ^ 2 + t12 ^ 2) * pkin(1) ^ 2, t14 ^ 2, t16 * t21, 0, 0, 0, -0.2e1 * t9 * t16, t9 * t21, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t22, t6 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t14, t16, 0, -t14 * t8, -t16 * t8, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, -0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
