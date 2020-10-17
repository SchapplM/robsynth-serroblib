% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:54:19
% EndTime: 2019-05-05 17:54:20
% DurationCPUTime: 0.25s
% Computational Cost: add. (189->61), mult. (240->76), div. (0->0), fcn. (233->8), ass. (0->42)
t26 = cos(pkin(9));
t24 = pkin(9) + qJ(3);
t21 = sin(t24);
t22 = cos(t24);
t41 = t22 * pkin(3) + t21 * qJ(4);
t59 = t26 * pkin(2) + pkin(1) + t41;
t32 = cos(qJ(1));
t58 = g(2) * t32;
t30 = sin(qJ(1));
t52 = g(1) * t32;
t14 = g(2) * t30 + t52;
t31 = cos(qJ(5));
t49 = g(3) * t22;
t42 = t32 * t31;
t29 = sin(qJ(5));
t45 = t30 * t29;
t5 = t21 * t42 - t45;
t43 = t32 * t29;
t44 = t30 * t31;
t7 = t21 * t44 + t43;
t57 = -g(1) * t5 - g(2) * t7 + t31 * t49;
t2 = g(3) * t21 + t14 * t22;
t54 = pkin(3) * t21;
t53 = pkin(5) * t29;
t27 = -qJ(6) - pkin(8);
t47 = t22 * t27;
t46 = t22 * t29;
t28 = -pkin(7) - qJ(2);
t40 = t31 * pkin(5) + pkin(4) - t28;
t39 = qJ(4) * t22;
t38 = pkin(5) * t46;
t36 = t59 * t58;
t13 = g(1) * t30 - t58;
t35 = t21 * t53 - t47;
t11 = t32 * t39;
t9 = t30 * t39;
t8 = -t21 * t45 + t42;
t6 = t21 * t43 + t44;
t4 = t13 * t22;
t3 = t13 * t21;
t1 = t14 * t21 - t49;
t10 = [0, t13, t14, t13 * t26, -t13 * sin(pkin(9)) -t14, -g(1) * (-t30 * pkin(1) + t32 * qJ(2)) - g(2) * (t32 * pkin(1) + t30 * qJ(2)) 0, 0, 0, 0, 0, t4, -t3, -t14, -t4, t3, t28 * t52 - t36 + (g(1) * t59 + g(2) * t28) * t30, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t4, -t36 + (-g(1) * t40 - g(2) * t35) * t32 + (-g(1) * (-t59 - t35) - g(2) * t40) * t30; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(1) * (-t32 * t54 + t11) - g(2) * (-t30 * t54 + t9) - g(3) * t41, 0, 0, 0, 0, 0, -t2 * t29, -t2 * t31, t1, -g(1) * (t32 * t38 + t11) - g(2) * (t30 * t38 + t9) - g(3) * (t41 - t47) + (-g(3) * t53 + t14 * (pkin(3) - t27)) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, g(1) * t6 - g(2) * t8 - g(3) * t46, 0, t57 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg  = t10;
