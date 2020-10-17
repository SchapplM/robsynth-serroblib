% Calculate minimal parameter regressor of gravitation load for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x14]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.08s
% Computational Cost: add. (37->16), mult. (40->23), div. (0->0), fcn. (36->6), ass. (0->11)
t4 = qJ(1) + pkin(6);
t2 = sin(t4);
t3 = cos(t4);
t11 = -g(1) * t3 - g(2) * t2;
t10 = g(1) * t2 - g(2) * t3;
t6 = sin(qJ(1));
t8 = cos(qJ(1));
t9 = g(1) * t6 - g(2) * t8;
t7 = cos(qJ(4));
t5 = sin(qJ(4));
t1 = [0, t9, g(1) * t8 + g(2) * t6, t9 * pkin(1), -t10, t11, -g(1) * (-pkin(1) * t6 - pkin(2) * t2 + qJ(3) * t3) - g(2) * (pkin(1) * t8 + t3 * pkin(2) + t2 * qJ(3)), 0, 0, 0, 0, 0, t11 * t5, t11 * t7; 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t5 - t10 * t7, g(3) * t7 + t10 * t5;];
taug_reg = t1;
