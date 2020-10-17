% Calculate minimal parameter regressor of gravitation load for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_gravloadJ_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:24
% EndTime: 2019-07-18 13:27:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (44->8), mult. (24->12), div. (0->0), fcn. (24->6), ass. (0->13)
t10 = qJ(2) + qJ(3);
t12 = cos(qJ(2));
t11 = sin(qJ(2));
t9 = qJ(4) + t10;
t8 = cos(t10);
t7 = sin(t10);
t6 = cos(t9);
t5 = sin(t9);
t4 = g(1) * t8 + g(3) * t7;
t3 = -g(1) * t7 + g(3) * t8;
t2 = g(1) * t6 + g(3) * t5;
t1 = -g(1) * t5 + g(3) * t6;
t13 = [g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t12 + g(3) * t11, -g(1) * t11 + g(3) * t12, 0, t4, t3, 0, t2, t1; 0, 0, 0, 0, 0, t4, t3, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg  = t13;
