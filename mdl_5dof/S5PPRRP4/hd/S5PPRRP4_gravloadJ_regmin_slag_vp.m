% Calculate minimal parameter regressor of gravitation load for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t17 = cos(qJ(3));
t16 = sin(qJ(3));
t15 = cos(pkin(7));
t14 = sin(pkin(7));
t1 = -t14 * t16 - t15 * t17;
t2 = -t14 * t17 + t15 * t16;
t13 = g(1) * t2 - g(2) * t1;
t12 = g(1) * t1 + g(2) * t2;
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t11 = g(3) * t10 - t12 * t9;
t8 = -qJ(5) - pkin(6);
t7 = t10 * pkin(4) + pkin(3);
t3 = -g(1) * t14 + g(2) * t15;
t4 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, t13, -t12, 0, 0, 0, 0, 0, t13 * t10, -t13 * t9, t12, -g(1) * (t1 * t8 - t2 * t7) - g(2) * (t1 * t7 + t2 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -g(3) * t9 - t12 * t10, 0, t11 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13;];
taug_reg = t4;
