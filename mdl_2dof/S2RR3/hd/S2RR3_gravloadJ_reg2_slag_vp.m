% Calculate inertial parameters regressor of gravitation load for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% taug_reg [2x(2*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S2RR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_gravloadJ_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_gravloadJ_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.13s
% Computational Cost: add. (15->5), mult. (16->9), div. (0->0), fcn. (14->4), ass. (0->9)
t6 = sin(qJ(1));
t7 = cos(qJ(1));
t8 = g(1) * t6 - g(2) * t7;
t5 = qJ(1) + qJ(2);
t4 = cos(t5);
t3 = sin(t5);
t2 = g(1) * t4 + g(2) * t3;
t1 = g(1) * t3 - g(2) * t4;
t9 = [0, 0, 0, 0, 0, 0, t8, g(1) * t7 + g(2) * t6, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t8 * pkin(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
