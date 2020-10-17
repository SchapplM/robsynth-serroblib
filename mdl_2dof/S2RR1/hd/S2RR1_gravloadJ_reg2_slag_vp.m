% Calculate inertial parameters regressor of gravitation load for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S2RR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_gravloadJ_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:10
% EndTime: 2020-01-03 11:19:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (10->6), mult. (28->11), div. (0->0), fcn. (26->4), ass. (0->7)
t3 = sin(qJ(1));
t5 = cos(qJ(1));
t6 = g(1) * t5 - g(3) * t3;
t1 = g(1) * t3 + g(3) * t5;
t4 = cos(qJ(2));
t2 = sin(qJ(2));
t7 = [0, 0, 0, 0, 0, 0, -t6, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t4, t6 * t2, t1, t1 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t4 + t1 * t2, -g(2) * t2 + t1 * t4, 0, 0;];
taug_reg = t7;
