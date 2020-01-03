% Calculate minimal parameter regressor of gravitation load for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% 
% Output:
% taug_reg [5x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t13 = cos(pkin(7));
t12 = sin(pkin(7));
t2 = -t12 * t14 - t13 * t15;
t3 = -t12 * t15 + t13 * t14;
t11 = g(1) * t3 - g(2) * t2;
t10 = g(1) * t2 + g(2) * t3;
t9 = cos(qJ(5));
t8 = sin(qJ(5));
t4 = -g(1) * t12 + g(2) * t13;
t1 = [-g(3), -g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t11, -t10, -t11, t10, -g(1) * (-t3 * pkin(3) - t2 * qJ(4)) - g(2) * (t2 * pkin(3) - t3 * qJ(4)), 0, 0, 0, 0, 0, t10 * t8, t10 * t9; 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t8 - t11 * t9, -g(3) * t9 + t11 * t8;];
taug_reg = t1;
