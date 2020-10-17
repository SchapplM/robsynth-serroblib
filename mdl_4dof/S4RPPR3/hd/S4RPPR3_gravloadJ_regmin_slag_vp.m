% Calculate minimal parameter regressor of gravitation load for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:56
% DurationCPUTime: 0.08s
% Computational Cost: add. (50->19), mult. (46->25), div. (0->0), fcn. (42->8), ass. (0->12)
t6 = qJ(1) + pkin(6);
t2 = sin(t6);
t4 = cos(t6);
t13 = g(1) * t4 + g(2) * t2;
t12 = g(1) * t2 - g(2) * t4;
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t11 = g(1) * t9 - g(2) * t10;
t5 = pkin(7) + qJ(4);
t3 = cos(t5);
t1 = sin(t5);
t7 = [0, t11, g(1) * t10 + g(2) * t9, t11 * pkin(1), t12 * cos(pkin(7)), -t12 * sin(pkin(7)), -t13, -g(1) * (-t9 * pkin(1) - t2 * pkin(2) + t4 * qJ(3)) - g(2) * (t10 * pkin(1) + t4 * pkin(2) + t2 * qJ(3)), 0, 0, 0, 0, 0, t12 * t3, -t12 * t1; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t3 + t13 * t1, g(3) * t1 + t13 * t3;];
taug_reg = t7;
