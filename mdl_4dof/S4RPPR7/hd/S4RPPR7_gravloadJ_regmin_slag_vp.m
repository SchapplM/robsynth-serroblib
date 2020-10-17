% Calculate minimal parameter regressor of gravitation load for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% taug_reg [4x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->21), mult. (56->22), div. (0->0), fcn. (52->6), ass. (0->10)
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t14 = t13 * pkin(1) + t12 * qJ(2);
t2 = g(1) * t13 + g(2) * t12;
t1 = g(1) * t12 - g(2) * t13;
t9 = pkin(6) + qJ(4);
t6 = t13 * qJ(2);
t4 = cos(t9);
t3 = sin(t9);
t5 = [0, t1, t2, -t1, -t2, -g(1) * (-t12 * pkin(1) + t6) - g(2) * t14, -t2 * sin(pkin(6)), -t2 * cos(pkin(6)), t1, -g(1) * (t6 + (-pkin(1) - qJ(3)) * t12) - g(2) * (t13 * qJ(3) + t14), 0, 0, 0, 0, 0, -t2 * t3, -t2 * t4; 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t3 - t1 * t4, g(3) * t4 + t1 * t3;];
taug_reg = t5;
