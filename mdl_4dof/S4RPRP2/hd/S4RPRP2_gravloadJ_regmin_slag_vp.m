% Calculate minimal parameter regressor of gravitation load for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_minpar_matlab.m
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t21 = t16 * pkin(1) + t14 * qJ(2);
t13 = sin(qJ(3));
t20 = t13 * t16;
t19 = t14 * t13;
t15 = cos(qJ(3));
t1 = -t16 * t15 - t19;
t2 = -t14 * t15 + t20;
t18 = g(1) * t2 - g(2) * t1;
t17 = g(1) * t1 + g(2) * t2;
t10 = t16 * qJ(2);
t8 = pkin(3) * t15 + pkin(2);
t4 = g(1) * t16 + g(2) * t14;
t3 = g(1) * t14 - g(2) * t16;
t5 = [0, t3, t4, t3, -t4, -g(1) * (-t14 * pkin(1) + t10) - g(2) * t21, 0, -t18, t17, -g(1) * (pkin(3) * t20 + t10 + (-pkin(1) - t8) * t14) - g(2) * (pkin(3) * t19 + t16 * t8 + t21); 0, 0, 0, 0, 0, -t3, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, t18, -t17, t18 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
taug_reg  = t5;
