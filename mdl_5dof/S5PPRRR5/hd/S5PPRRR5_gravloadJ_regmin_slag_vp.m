% Calculate minimal parameter regressor of gravitation load for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (81->16), mult. (88->24), div. (0->0), fcn. (112->8), ass. (0->18)
t20 = sin(pkin(8));
t19 = qJ(3) + qJ(4);
t18 = cos(t19);
t17 = sin(t19);
t12 = cos(pkin(8));
t5 = -t12 * t18 - t20 * t17;
t6 = t12 * t17 - t20 * t18;
t4 = g(1) * t6 - g(2) * t5;
t3 = -g(1) * t5 - g(2) * t6;
t16 = cos(qJ(3));
t15 = cos(qJ(5));
t14 = sin(qJ(3));
t13 = sin(qJ(5));
t8 = t12 * t14 - t20 * t16;
t7 = -t12 * t16 - t20 * t14;
t2 = t4 * t15;
t1 = t4 * t13;
t9 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t20 + g(2) * t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, g(1) * t8 - g(2) * t7, -g(1) * t7 - g(2) * t8, 0, t4, t3, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t15 + t3 * t13, -g(3) * t13 + t3 * t15;];
taug_reg = t9;
