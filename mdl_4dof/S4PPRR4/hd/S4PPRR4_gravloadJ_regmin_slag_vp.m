% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:40
% EndTime: 2019-12-31 16:18:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (37->16), mult. (48->24), div. (0->0), fcn. (52->6), ass. (0->15)
t3 = pkin(7) + qJ(3);
t1 = sin(t3);
t2 = cos(t3);
t4 = sin(pkin(6));
t5 = cos(pkin(6));
t9 = g(1) * t5 + g(2) * t4;
t8 = -g(3) * t2 + t9 * t1;
t15 = g(3) * t1;
t6 = sin(qJ(4));
t13 = t4 * t6;
t7 = cos(qJ(4));
t12 = t4 * t7;
t11 = t5 * t6;
t10 = t5 * t7;
t14 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t4 + g(2) * t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t8, t9 * t2 + t15, 0, 0, 0, 0, 0, t8 * t7, -t8 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t2 * t11 + t12) - g(2) * (-t2 * t13 - t10) + t6 * t15, -g(1) * (-t2 * t10 - t13) - g(2) * (-t2 * t12 + t11) + t7 * t15;];
taug_reg = t14;
