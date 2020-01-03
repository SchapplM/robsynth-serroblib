% Calculate minimal parameter regressor of gravitation load for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t4 = sin(qJ(3));
t6 = cos(qJ(3));
t1 = sin(pkin(6));
t2 = cos(pkin(6));
t8 = g(1) * t1 - g(2) * t2;
t13 = -g(3) * t4 + t8 * t6;
t11 = g(3) * t6;
t3 = sin(qJ(4));
t10 = t3 * t4;
t5 = cos(qJ(4));
t9 = t4 * t5;
t7 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t13, t8 * t4 + t11, 0, 0, 0, 0, 0, -t13 * t5, t13 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t1 * t10 + t2 * t5) - g(2) * (t1 * t5 + t2 * t10) + t3 * t11, -g(1) * (-t1 * t9 - t2 * t3) - g(2) * (-t1 * t3 + t2 * t9) + t5 * t11;];
taug_reg = t7;
