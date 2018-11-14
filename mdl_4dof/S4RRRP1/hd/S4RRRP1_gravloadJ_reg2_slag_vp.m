% Calculate inertial parameters regressor of gravitation load for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t13 = qJ(1) + qJ(2);
t9 = sin(t13);
t20 = pkin(2) * t9;
t11 = qJ(3) + t13;
t7 = sin(t11);
t19 = pkin(3) * t7;
t15 = cos(qJ(1));
t10 = cos(t13);
t6 = pkin(2) * t10;
t18 = t15 * pkin(1) + t6;
t8 = cos(t11);
t1 = g(1) * t7 - g(2) * t8;
t14 = sin(qJ(1));
t17 = -t14 * pkin(1) - t20;
t3 = g(1) * t9 - g(2) * t10;
t16 = g(1) * t14 - g(2) * t15;
t5 = pkin(3) * t8;
t4 = g(1) * t10 + g(2) * t9;
t2 = g(1) * t8 + g(2) * t7;
t12 = [0, 0, 0, 0, 0, 0, t16, g(1) * t15 + g(2) * t14, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t16 * pkin(1), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * t17 - g(2) * t18, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t17 - t19) - g(2) * (t5 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t19 - t20) - g(2) * (t5 + t6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
taug_reg  = t12;
