% Calculate minimal parameter regressor of gravitation load for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x8]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:59:54
% EndTime: 2019-05-04 18:59:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (26->8), mult. (16->9), div. (0->0), fcn. (14->4), ass. (0->9)
t6 = sin(qJ(2));
t7 = cos(qJ(2));
t8 = g(1) * t6 - g(2) * t7;
t5 = qJ(2) + pkin(6) + qJ(4);
t4 = cos(t5);
t3 = sin(t5);
t2 = g(1) * t4 + g(2) * t3;
t1 = g(1) * t3 - g(2) * t4;
t9 = [-g(2), 0, 0, 0, -g(2), 0, 0, 0; 0, 0, t8, g(1) * t7 + g(2) * t6, t8 * pkin(2), 0, t1, t2; 0, 0, 0, 0, -g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t9;
