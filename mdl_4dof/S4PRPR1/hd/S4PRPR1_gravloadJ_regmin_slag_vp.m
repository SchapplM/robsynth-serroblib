% Calculate minimal parameter regressor of gravitation load for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_minpar_matlab.m
t14 = cos(qJ(4));
t13 = sin(qJ(4));
t10 = pkin(6) + qJ(2);
t8 = sin(t10);
t9 = cos(t10);
t1 = -t8 * t13 - t9 * t14;
t2 = t9 * t13 - t8 * t14;
t12 = g(1) * t2 - g(2) * t1;
t11 = g(1) * t1 + g(2) * t2;
t4 = g(1) * t9 + g(2) * t8;
t3 = g(1) * t8 - g(2) * t9;
t5 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0; 0, 0, t3, t4, t3, -t4, -g(1) * (-pkin(2) * t8 + qJ(3) * t9) - g(2) * (pkin(2) * t9 + qJ(3) * t8) 0, -t12, t11; 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11;];
taug_reg  = t5;
