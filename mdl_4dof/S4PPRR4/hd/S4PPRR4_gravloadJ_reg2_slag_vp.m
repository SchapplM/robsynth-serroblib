% Calculate inertial parameters regressor of gravitation load for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (57->22), mult. (70->30), div. (0->0), fcn. (71->6), ass. (0->17)
t6 = sin(pkin(6));
t7 = cos(pkin(6));
t11 = g(1) * t7 + g(2) * t6;
t5 = pkin(7) + qJ(3);
t3 = sin(t5);
t4 = cos(t5);
t10 = -g(3) * t4 + t11 * t3;
t18 = g(3) * t3;
t8 = sin(qJ(4));
t16 = t6 * t8;
t9 = cos(qJ(4));
t15 = t6 * t9;
t14 = t7 * t8;
t13 = t7 * t9;
t2 = -g(1) * t6 + g(2) * t7;
t1 = t11 * t4 + t18;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t1, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t9, -t10 * t8, -t1, -g(3) * (t4 * pkin(3) + t3 * pkin(5)) + t11 * (pkin(3) * t3 - pkin(5) * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t14 + t15) - g(2) * (-t4 * t16 - t13) + t8 * t18, -g(1) * (-t4 * t13 - t16) - g(2) * (-t4 * t15 + t14) + t9 * t18, 0, 0;];
taug_reg = t12;
