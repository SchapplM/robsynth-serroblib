% Calculate minimal parameter regressor of gravitation load for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x8]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = cos(qJ(3));
t12 = sin(qJ(3));
t11 = cos(pkin(6));
t10 = sin(pkin(6));
t9 = qJ(3) + qJ(4);
t8 = cos(t9);
t7 = sin(t9);
t6 = -t10 * t13 + t11 * t12;
t5 = -t10 * t12 - t11 * t13;
t4 = -t10 * t8 + t11 * t7;
t3 = -t10 * t7 - t11 * t8;
t2 = g(1) * t4 - g(2) * t3;
t1 = -g(1) * t3 - g(2) * t4;
t14 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0; 0, -g(1) * t10 + g(2) * t11, 0, 0, 0, 0, 0, 0; 0, 0, 0, g(1) * t6 - g(2) * t5, -g(1) * t5 - g(2) * t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg  = t14;
