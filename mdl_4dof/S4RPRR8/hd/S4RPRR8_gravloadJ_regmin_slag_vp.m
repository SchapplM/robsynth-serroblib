% Calculate minimal parameter regressor of gravitation load for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% taug_reg [4x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (40->17), mult. (62->22), div. (0->0), fcn. (60->6), ass. (0->12)
t10 = sin(qJ(1));
t12 = cos(qJ(1));
t4 = g(1) * t12 + g(2) * t10;
t3 = g(1) * t10 - g(2) * t12;
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t8 = qJ(3) + qJ(4);
t6 = cos(t8);
t5 = sin(t8);
t2 = g(3) * t5 - t3 * t6;
t1 = g(3) * t6 + t3 * t5;
t7 = [0, t3, t4, -t3, -t4, -g(1) * (-t10 * pkin(1) + t12 * qJ(2)) - g(2) * (t12 * pkin(1) + t10 * qJ(2)), 0, 0, 0, 0, 0, -t4 * t9, -t4 * t11, 0, 0, 0, 0, 0, -t4 * t5, -t4 * t6; 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t9 - t3 * t11, g(3) * t11 + t3 * t9, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t7;
