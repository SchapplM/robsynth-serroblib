% Calculate minimal parameter regressor of gravitation load for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% taug_reg [4x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t4 = qJ(1) + pkin(6);
t2 = sin(t4);
t3 = cos(t4);
t13 = g(1) * t3 + g(2) * t2;
t12 = g(1) * t2 - g(2) * t3;
t7 = sin(qJ(1));
t9 = cos(qJ(1));
t11 = g(1) * t7 - g(2) * t9;
t6 = sin(qJ(3));
t8 = cos(qJ(3));
t10 = -g(3) * t8 + t13 * t6;
t5 = -qJ(4) - pkin(5);
t1 = t8 * pkin(3) + pkin(2);
t14 = [0, t11, g(1) * t9 + g(2) * t7, t11 * pkin(1), 0, 0, 0, 0, 0, t12 * t8, -t12 * t6, -t13, -g(1) * (-t7 * pkin(1) - t2 * t1 - t3 * t5) - g(2) * (t9 * pkin(1) + t3 * t1 - t2 * t5); 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, g(3) * t6 + t13 * t8, 0, t10 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12;];
taug_reg = t14;
