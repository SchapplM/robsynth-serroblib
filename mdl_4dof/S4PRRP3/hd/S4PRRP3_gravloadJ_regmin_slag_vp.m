% Calculate minimal parameter regressor of gravitation load for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x13]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t6 = pkin(6) + qJ(2);
t4 = sin(t6);
t5 = cos(t6);
t2 = g(1) * t5 + g(2) * t4;
t1 = g(1) * t4 - g(2) * t5;
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t10 = -g(3) * t9 + t2 * t8;
t7 = -qJ(4) - pkin(5);
t3 = t9 * pkin(3) + pkin(2);
t11 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t1, t2, 0, 0, 0, 0, 0, t1 * t9, -t1 * t8, -t2, -g(1) * (-t4 * t3 - t5 * t7) - g(2) * (t5 * t3 - t4 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, g(3) * t8 + t2 * t9, 0, t10 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
