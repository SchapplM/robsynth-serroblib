% Calculate minimal parameter regressor of gravitation load for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t9 = qJ(1) + pkin(6);
t8 = qJ(3) + t9;
t6 = sin(t8);
t7 = cos(t8);
t14 = t7 * pkin(3) + t6 * qJ(4);
t13 = -t6 * pkin(3) + t7 * qJ(4);
t10 = sin(qJ(1));
t11 = cos(qJ(1));
t12 = g(1) * t10 - g(2) * t11;
t2 = g(1) * t7 + g(2) * t6;
t1 = g(1) * t6 - g(2) * t7;
t3 = [0, t12, g(1) * t11 + g(2) * t10, t12 * pkin(1), 0, t1, t2, t1, -t2, -g(1) * (-pkin(2) * sin(t9) - t10 * pkin(1) + t13) - g(2) * (pkin(2) * cos(t9) + t11 * pkin(1) + t14); 0, 0, 0, -g(3), 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, t1, t2, t1, -t2, -g(1) * t13 - g(2) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t3;
