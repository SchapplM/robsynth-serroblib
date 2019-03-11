% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(4));
t21 = pkin(1) + qJ(3);
t26 = cos(qJ(4));
t31 = -t26 * qJ(5) + t21;
t5 = t24 * pkin(4) + t31;
t40 = -0.2e1 * t5;
t39 = 0.2e1 * t21;
t38 = 0.2e1 * qJ(5);
t37 = pkin(4) + pkin(8);
t23 = sin(qJ(6));
t11 = t23 * t24;
t25 = cos(qJ(6));
t36 = t25 * t23;
t13 = t25 * t24;
t20 = -pkin(7) + qJ(2);
t35 = t26 * t20;
t34 = t26 * t24;
t17 = t24 ^ 2;
t19 = t26 ^ 2;
t9 = t17 + t19;
t33 = t24 * qJ(5);
t32 = 0.2e1 * t34;
t8 = pkin(4) * t26 + t33;
t30 = t26 * t37 + t33;
t29 = (qJ(2) ^ 2);
t28 = 2 * qJ(2);
t18 = t25 ^ 2;
t16 = t23 ^ 2;
t14 = t25 * t26;
t12 = t23 * t26;
t10 = t24 * t20;
t7 = (pkin(5) - t20) * t26;
t6 = -pkin(5) * t24 + t10;
t4 = t9 * t20;
t3 = t24 * t37 + t31;
t2 = t23 * t7 + t25 * t3;
t1 = -t23 * t3 + t25 * t7;
t15 = [1, 0, 0, -2 * pkin(1), t28, pkin(1) ^ 2 + t29, t28, t39, t21 ^ 2 + t29, t19, -0.2e1 * t34, 0, 0, 0, t24 * t39, t26 * t39, -0.2e1 * t4, t24 * t40, t26 * t40, t20 ^ 2 * t9 + t5 ^ 2, t16 * t17, 0.2e1 * t17 * t36, t23 * t32, t25 * t32, t19, 0.2e1 * t1 * t26 - 0.2e1 * t13 * t6, 0.2e1 * t11 * t6 - 0.2e1 * t2 * t26; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t21, 0, 0, 0, 0, 0, -t24, -t26, 0, t24, t26, -t5, 0, 0, 0, 0, 0, t12, t14; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, t4, 0, 0, 0, 0, 0, -t9 * t25, t9 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, t35, -t10, -t8, -t35, t10, t8 * t20, t23 * t13 (-t16 + t18) * t24, t14, -t12, 0, t6 * t23 - t25 * t30, t23 * t30 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, -t26, t24, t8, 0, 0, 0, 0, 0, t11, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t38, pkin(4) ^ 2 + qJ(5) ^ 2, t18, -0.2e1 * t36, 0, 0, 0, t23 * t38, t25 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, -t35, 0, 0, 0, 0, 0, t14, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t13, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23, 0, -t25 * t37, t23 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
