% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:02:55
% DurationCPUTime: 0.10s
% Computational Cost: add. (147->24), mult. (102->31), div. (0->0), fcn. (98->9), ass. (0->23)
t21 = qJ(1) + qJ(2);
t18 = sin(t21);
t19 = cos(t21);
t26 = t19 * pkin(2) + t18 * qJ(3);
t20 = pkin(9) + qJ(4);
t25 = -t18 * pkin(2) + t19 * qJ(3);
t9 = g(1) * t19 + g(2) * t18;
t8 = g(1) * t18 - g(2) * t19;
t24 = cos(qJ(1));
t23 = sin(qJ(1));
t17 = qJ(5) + t20;
t16 = cos(t20);
t15 = sin(t20);
t13 = cos(t17);
t12 = sin(t17);
t7 = t8 * cos(pkin(9));
t6 = t8 * t16;
t5 = t8 * t15;
t4 = t8 * t13;
t3 = t8 * t12;
t2 = g(3) * t12 + t9 * t13;
t1 = -g(3) * t13 + t9 * t12;
t10 = [0, g(1) * t23 - g(2) * t24, g(1) * t24 + g(2) * t23, 0, t8, t9, t7, -t9, -g(1) * (-t23 * pkin(1) + t25) - g(2) * (t24 * pkin(1) + t26), 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, t8, t9, t7, -t9, -g(1) * t25 - g(2) * t26, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t16 + t9 * t15, g(3) * t15 + t9 * t16, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t10;
