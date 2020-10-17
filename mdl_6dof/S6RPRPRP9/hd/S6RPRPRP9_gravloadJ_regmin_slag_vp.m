% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:06:53
% EndTime: 2019-05-05 18:06:53
% DurationCPUTime: 0.27s
% Computational Cost: add. (216->73), mult. (334->98), div. (0->0), fcn. (349->8), ass. (0->46)
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t55 = g(1) * t29;
t58 = g(2) * t31 - t55;
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t10 = -g(3) * t28 - t58 * t30;
t56 = -pkin(1) - pkin(7);
t52 = g(3) * t30;
t50 = t28 * t31;
t24 = pkin(9) + qJ(5);
t17 = sin(t24);
t49 = t29 * t17;
t18 = cos(t24);
t48 = t29 * t18;
t25 = sin(pkin(9));
t47 = t29 * t25;
t26 = cos(pkin(9));
t46 = t29 * t26;
t27 = -pkin(8) - qJ(4);
t45 = t30 * t27;
t44 = t31 * t17;
t43 = t31 * t18;
t42 = t31 * t25;
t41 = t31 * t26;
t40 = t31 * pkin(1) + t29 * qJ(2);
t39 = t30 * qJ(4);
t38 = t31 * pkin(7) + t40;
t5 = t28 * t49 - t43;
t7 = t28 * t44 + t48;
t37 = g(1) * t7 + g(2) * t5;
t13 = g(1) * t31 + g(2) * t29;
t35 = t28 * pkin(3) - t39;
t16 = t26 * pkin(4) + pkin(3);
t34 = pkin(5) * t18 + qJ(6) * t17 + t16;
t1 = g(1) * t5 - g(2) * t7 + t17 * t52;
t6 = t28 * t48 + t44;
t8 = t28 * t43 - t49;
t33 = g(1) * t6 - g(2) * t8 + t18 * t52;
t20 = t31 * qJ(2);
t11 = t13 * t30;
t9 = -g(2) * t50 + t28 * t55 + t52;
t4 = t10 * t18;
t3 = t10 * t17;
t2 = -g(1) * t8 - g(2) * t6;
t12 = [0, -t58, t13, t58, -t13, -g(1) * (-t29 * pkin(1) + t20) - g(2) * t40, 0, 0, 0, 0, 0, -t13 * t28, -t11, -g(1) * (t28 * t41 - t47) - g(2) * (t28 * t46 + t42) -g(1) * (-t28 * t42 - t46) - g(2) * (-t28 * t47 + t41) t11, -g(1) * (pkin(3) * t50 - t31 * t39 + t20) - g(2) * t38 + (-g(1) * t56 - g(2) * t35) * t29, 0, 0, 0, 0, 0, t2, t37, t2, t11, -t37, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t16 * t50 + t31 * t45 + t20) - g(2) * (pkin(4) * t42 + t6 * pkin(5) + t5 * qJ(6) + t38) + (-g(1) * (-pkin(4) * t25 + t56) - g(2) * (t28 * t16 + t45)) * t29; 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, -t10 * t26, t10 * t25, -t9, g(3) * t35 + t58 * (pkin(3) * t30 + qJ(4) * t28) 0, 0, 0, 0, 0, -t4, t3, -t4, -t9, -t3, -g(3) * (-t34 * t28 - t45) + t58 * (-t27 * t28 + t34 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t33, t1, 0, -t33, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (t7 * pkin(5) - t8 * qJ(6)) - (-pkin(5) * t17 + qJ(6) * t18) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
