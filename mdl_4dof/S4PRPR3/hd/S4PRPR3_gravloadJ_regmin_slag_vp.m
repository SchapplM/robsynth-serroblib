% Calculate minimal parameter regressor of gravitation load for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->15), mult. (40->18), div. (0->0), fcn. (38->6), ass. (0->9)
t8 = pkin(6) + qJ(2);
t4 = sin(t8);
t6 = cos(t8);
t2 = g(1) * t6 + g(2) * t4;
t1 = g(1) * t4 - g(2) * t6;
t7 = pkin(7) + qJ(4);
t5 = cos(t7);
t3 = sin(t7);
t9 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t1, t2, t1 * cos(pkin(7)), -t1 * sin(pkin(7)), -t2, -g(1) * (-t4 * pkin(2) + t6 * qJ(3)) - g(2) * (t6 * pkin(2) + t4 * qJ(3)), 0, 0, 0, 0, 0, t1 * t5, -t1 * t3; 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t5 + t2 * t3, g(3) * t3 + t2 * t5;];
taug_reg = t9;
