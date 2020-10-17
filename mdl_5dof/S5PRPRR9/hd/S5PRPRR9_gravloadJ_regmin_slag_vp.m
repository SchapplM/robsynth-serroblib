% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:45
% DurationCPUTime: 0.09s
% Computational Cost: add. (96->20), mult. (90->24), div. (0->0), fcn. (108->6), ass. (0->16)
t14 = pkin(8) + qJ(2);
t12 = sin(t14);
t13 = cos(t14);
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t1 = -t12 * t15 - t13 * t16;
t2 = -t12 * t16 + t13 * t15;
t11 = g(1) * t2 - g(2) * t1;
t8 = sin(qJ(5));
t18 = t11 * t8;
t9 = cos(qJ(5));
t17 = t11 * t9;
t10 = g(1) * t1 + g(2) * t2;
t4 = g(1) * t13 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t13;
t5 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, t3, -t4, -g(1) * (-t12 * pkin(2) + t13 * qJ(3)) - g(2) * (t13 * pkin(2) + t12 * qJ(3)), 0, -t11, t10, 0, 0, 0, 0, 0, -t17, t18; 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, 0, 0, 0, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t9 - t10 * t8, -g(3) * t8 - t10 * t9;];
taug_reg = t5;
