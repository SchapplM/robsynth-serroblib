% Calculate minimal parameter regressor of gravitation load for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (61->16), mult. (70->24), div. (0->0), fcn. (65->6), ass. (0->15)
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t4 = g(1) * t13 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t13;
t10 = sin(qJ(2));
t12 = cos(qJ(2));
t14 = -g(3) * t12 + t4 * t10;
t9 = -qJ(3) - pkin(5);
t8 = qJ(2) + pkin(7) + qJ(4);
t7 = t12 * pkin(2) + pkin(1);
t6 = cos(t8);
t5 = sin(t8);
t2 = g(3) * t5 + t4 * t6;
t1 = -g(3) * t6 + t4 * t5;
t15 = [0, t3, t4, 0, 0, 0, 0, 0, t3 * t12, -t3 * t10, -t4, -g(1) * (-t11 * t7 - t13 * t9) - g(2) * (-t11 * t9 + t13 * t7), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, t14, g(3) * t10 + t4 * t12, 0, t14 * pkin(2), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t15;
