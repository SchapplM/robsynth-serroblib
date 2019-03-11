% Calculate minimal parameter regressor of gravitation load for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t11 = qJ(1) + qJ(2);
t8 = sin(t11);
t9 = cos(t11);
t18 = t9 * pkin(2) + t8 * qJ(3);
t12 = sin(qJ(1));
t17 = t12 * pkin(1);
t13 = cos(qJ(1));
t16 = t13 * pkin(1) + t18;
t4 = t9 * qJ(3);
t15 = -t8 * pkin(2) + t4;
t14 = t4 + (-pkin(2) - pkin(3)) * t8;
t5 = t9 * pkin(3);
t2 = g(1) * t9 + g(2) * t8;
t1 = g(1) * t8 - g(2) * t9;
t3 = [0, g(1) * t12 - g(2) * t13, g(1) * t13 + g(2) * t12, 0, t1, t2, t1, -t2, -g(1) * (t15 - t17) - g(2) * t16, t1, -t2, -g(1) * (t14 - t17) - g(2) * (t5 + t16); 0, 0, 0, 0, t1, t2, t1, -t2, -g(1) * t15 - g(2) * t18, t1, -t2, -g(1) * t14 - g(2) * (t5 + t18); 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3);];
taug_reg  = t3;
