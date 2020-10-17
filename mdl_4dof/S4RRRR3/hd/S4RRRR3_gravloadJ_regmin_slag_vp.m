% Calculate minimal parameter regressor of gravitation load for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:37
% DurationCPUTime: 0.07s
% Computational Cost: add. (92->13), mult. (88->22), div. (0->0), fcn. (88->8), ass. (0->17)
t10 = qJ(2) + qJ(3);
t12 = sin(qJ(1));
t14 = cos(qJ(1));
t16 = g(1) * t14 + g(2) * t12;
t15 = g(1) * t12 - g(2) * t14;
t13 = cos(qJ(2));
t11 = sin(qJ(2));
t9 = qJ(4) + t10;
t8 = cos(t10);
t7 = sin(t10);
t6 = cos(t9);
t5 = sin(t9);
t4 = g(3) * t7 + t16 * t8;
t3 = -g(3) * t8 + t16 * t7;
t2 = g(3) * t5 + t16 * t6;
t1 = -g(3) * t6 + t16 * t5;
t17 = [0, t15, t16, 0, 0, 0, 0, 0, t15 * t13, -t15 * t11, 0, 0, 0, 0, 0, t15 * t8, -t15 * t7, 0, 0, 0, 0, 0, t15 * t6, -t15 * t5; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t13 + t16 * t11, g(3) * t11 + t16 * t13, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t17;
