% Calculate minimal parameter regressor of gravitation load for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% taug_reg [4x8]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = cos(qJ(3));
t11 = sin(qJ(3));
t10 = sin(pkin(5));
t8 = cos(pkin(5));
t2 = -t10 * t11 - t8 * t12;
t3 = -t10 * t12 + t8 * t11;
t1 = g(1) * t3 - g(2) * t2;
t9 = g(1) * t2 + g(2) * t3;
t4 = -g(1) * t10 + g(2) * t8;
t5 = [-g(3), -g(3), 0, 0, 0, 0, 0, -g(3); 0, t4, 0, 0, 0, 0, 0, t4; 0, 0, 0, t1, -t9, t1, t9, -g(1) * (-t3 * pkin(3) - t2 * qJ(4)) - g(2) * (t2 * pkin(3) - t3 * qJ(4)); 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
