% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (146->29), mult. (98->43), div. (0->0), fcn. (104->10), ass. (0->25)
t28 = g(3) * sin(pkin(9));
t18 = cos(pkin(9));
t19 = sin(qJ(5));
t27 = t18 * t19;
t21 = cos(qJ(5));
t26 = t18 * t21;
t16 = qJ(1) + pkin(8);
t15 = qJ(3) + t16;
t13 = sin(t15);
t14 = cos(t15);
t25 = t14 * pkin(3) + t13 * qJ(4);
t24 = -t13 * pkin(3) + t14 * qJ(4);
t8 = g(1) * t13 - g(2) * t14;
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t23 = g(1) * t20 - g(2) * t22;
t9 = g(1) * t14 + g(2) * t13;
t7 = t8 * t18;
t6 = t13 * t19 + t14 * t26;
t5 = t13 * t21 - t14 * t27;
t4 = -t13 * t26 + t14 * t19;
t3 = t13 * t27 + t14 * t21;
t2 = -g(1) * t4 - g(2) * t6;
t1 = -g(1) * t3 - g(2) * t5;
t10 = [0, t23, g(1) * t22 + g(2) * t20, t23 * pkin(1), 0, t8, t9, t7, -t9, -g(1) * (-pkin(2) * sin(t16) - t20 * pkin(1) + t24) - g(2) * (pkin(2) * cos(t16) + t22 * pkin(1) + t25), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t8, t9, t7, -t9, -g(1) * t24 - g(2) * t25, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t19 * t28, g(1) * t6 - g(2) * t4 + t21 * t28;];
taug_reg = t10;
