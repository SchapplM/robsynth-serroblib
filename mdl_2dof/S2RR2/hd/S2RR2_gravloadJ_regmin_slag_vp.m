% Calculate minimal parameter regressor of gravitation load for
% S2RR2
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
% taug_reg [2x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S2RR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_gravloadJ_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_gravloadJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:19:15
% EndTime: 2019-05-04 18:19:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (8->5), mult. (22->10), div. (0->0), fcn. (22->4), ass. (0->7)
t2 = sin(qJ(1));
t4 = cos(qJ(1));
t6 = g(1) * t4 - g(3) * t2;
t5 = g(1) * t2 + g(3) * t4;
t3 = cos(qJ(2));
t1 = sin(qJ(2));
t7 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t3, -t5 * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 + t1 * t6, g(2) * t1 + t3 * t6;];
taug_reg  = t7;
