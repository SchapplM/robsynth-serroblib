% Calculate minimal parameter regressor of gravitation load for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% taug_reg [4x6]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:35:36
% EndTime: 2019-05-04 18:35:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (12->9), mult. (18->12), div. (0->0), fcn. (22->4), ass. (0->8)
t7 = cos(qJ(4));
t6 = sin(qJ(4));
t5 = cos(pkin(5));
t4 = sin(pkin(5));
t3 = -g(1) * t4 + g(2) * t5;
t2 = t4 * t7 + t5 * t6;
t1 = -t4 * t6 + t5 * t7;
t8 = [-g(3), -g(3), -g(3), 0, 0, 0; 0, t3, t3, 0, 0, 0; 0, 0, -g(1) * t5 - g(2) * t4, 0, 0, 0; 0, 0, 0, 0, -g(1) * t1 - g(2) * t2, g(1) * t2 - g(2) * t1;];
taug_reg  = t8;
