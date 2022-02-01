% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:46
% DurationCPUTime: 0.17s
% Computational Cost: add. (125->30), mult. (234->44), div. (0->0), fcn. (222->8), ass. (0->31)
t19 = sin(pkin(9));
t34 = pkin(1) * t19;
t21 = sin(qJ(5));
t33 = pkin(4) * t21;
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t20 = cos(pkin(9));
t13 = t20 * pkin(1) + pkin(2);
t23 = sin(qJ(3));
t26 = cos(qJ(3));
t8 = t26 * t13 - t23 * t34;
t7 = pkin(3) + t8;
t9 = t23 * t13 + t26 * t34;
t28 = t22 * t9 - t25 * t7;
t2 = -pkin(4) + t28;
t24 = cos(qJ(5));
t32 = t2 * t24;
t31 = t22 * pkin(3);
t30 = t25 * t9;
t16 = t25 * pkin(3);
t15 = -t16 - pkin(4);
t29 = t15 * t24;
t5 = -t22 * t7 - t30;
t18 = t21 ^ 2;
t17 = pkin(4) * t24;
t14 = pkin(8) + t31;
t12 = 0.2e1 * t21 * t24;
t11 = t15 * t21;
t3 = pkin(8) - t5;
t1 = t2 * t21;
t4 = [1, 0, 0, (t19 ^ 2 + t20 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t8, -0.2e1 * t9, 1, -0.2e1 * t28, 0.2e1 * t5, t18, t12, 0, 0, 0, -0.2e1 * t32, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t8, -t9, 1, t16 - t28, -t30 + (-pkin(3) - t7) * t22, t18, t12, 0, 0, 0, (-t15 - t2) * t24, t11 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t16, -0.2e1 * t31, t18, t12, 0, 0, 0, -0.2e1 * t29, 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 1, -t28, t5, t18, t12, 0, 0, 0, t17 - t32, t1 - t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t16, -t31, t18, t12, 0, 0, 0, t17 - t29, t11 - t33; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t18, t12, 0, 0, 0, 0.2e1 * t17, -0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t24, 0, -t21 * t3, -t24 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t24, 0, -t21 * t14, -t24 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t24, 0, -t21 * pkin(8), -t24 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
