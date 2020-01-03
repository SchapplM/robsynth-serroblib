% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = sin(qJ(3));
t26 = sin(qJ(2));
t28 = cos(qJ(3));
t29 = cos(qJ(2));
t13 = t25 * t26 - t28 * t29;
t21 = -t29 * pkin(2) - pkin(1);
t36 = 0.2e1 * t13 * pkin(3) + 0.2e1 * t21;
t35 = 0.2e1 * t21;
t34 = 0.2e1 * t29;
t33 = pkin(5) + pkin(6);
t24 = sin(qJ(4));
t32 = t24 * pkin(3);
t31 = t25 * pkin(2);
t27 = cos(qJ(4));
t30 = t27 * t31;
t16 = t33 * t26;
t17 = t33 * t29;
t7 = -t28 * t16 - t25 * t17;
t23 = t28 * pkin(2);
t20 = t23 + pkin(3);
t10 = t27 * t20 - t24 * t31;
t8 = t25 * t16 - t28 * t17;
t22 = t27 * pkin(3);
t14 = t25 * t29 + t28 * t26;
t11 = -t24 * t20 - t30;
t6 = -t24 * t13 + t27 * t14;
t5 = t27 * t13 + t24 * t14;
t4 = -t13 * pkin(7) - t8;
t3 = -t14 * pkin(7) + t7;
t2 = -t24 * t3 - t27 * t4;
t1 = -t24 * t4 + t27 * t3;
t9 = [1, 0, 0, t26 ^ 2, t26 * t34, 0, 0, 0, pkin(1) * t34, -0.2e1 * pkin(1) * t26, t14 ^ 2, -0.2e1 * t14 * t13, 0, 0, 0, t13 * t35, t14 * t35, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t36, t6 * t36; 0, 0, 0, 0, 0, t26, t29, 0, -t26 * pkin(5), -t29 * pkin(5), 0, 0, t14, -t13, 0, t7, t8, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t23, -0.2e1 * t31, 0, 0, 0, 0, 1, 0.2e1 * t10, 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t7, t8, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t23, -t31, 0, 0, 0, 0, 1, t10 + t22, -t30 + (-pkin(3) - t20) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t22, -0.2e1 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t22, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
