% Calculate minimal parameter regressor of gravitation load for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_minpar_matlab.m
t24 = sin(pkin(4));
t36 = pkin(3) * t24;
t26 = sin(qJ(1));
t27 = cos(qJ(1));
t34 = qJ(2) * t24;
t35 = t27 * pkin(1) + t26 * t34;
t22 = pkin(4) - pkin(6);
t16 = -sin(t22) / 0.2e1;
t21 = pkin(4) + pkin(6);
t18 = sin(t21);
t33 = t18 / 0.2e1 + t16;
t32 = -t26 * pkin(1) + t27 * t34;
t17 = cos(t22) / 0.2e1;
t19 = cos(t21);
t13 = t17 + t19 / 0.2e1;
t23 = sin(pkin(6));
t6 = -t27 * t13 + t26 * t23;
t8 = t26 * t13 + t27 * t23;
t2 = g(1) * t6 - g(2) * t8;
t25 = cos(pkin(6));
t7 = t26 * t25 + t27 * t33;
t9 = t27 * t25 - t26 * t33;
t3 = g(1) * t7 - g(2) * t9;
t31 = g(1) * t27 + g(2) * t26;
t30 = g(1) * t26 - g(2) * t27;
t29 = t9 * pkin(2) + t8 * qJ(3) + t35;
t28 = -t7 * pkin(2) - t6 * qJ(3) + t32;
t12 = t31 * t24;
t10 = -g(3) * cos(pkin(4)) - t30 * t24;
t1 = -g(1) * t8 - g(2) * t6 - g(3) * (-t18 / 0.2e1 + t16);
t4 = [0, t30, t31, t3, -t2, -t12, -g(1) * t32 - g(2) * t35, -t12, -t3, t2, -g(1) * t28 - g(2) * t29, -t12, t2, t3, -g(1) * (-t7 * qJ(4) + t27 * t36 + t28) - g(2) * (t9 * qJ(4) + t26 * t36 + t29); 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t17 - t19 / 0.2e1);];
taug_reg  = t4;
