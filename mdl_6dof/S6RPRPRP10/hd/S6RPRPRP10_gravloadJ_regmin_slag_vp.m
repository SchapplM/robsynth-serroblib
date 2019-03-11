% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t38 = pkin(5) * t29 - qJ(6) * t32;
t34 = cos(qJ(1));
t52 = g(2) * t34;
t31 = sin(qJ(1));
t53 = g(1) * t31;
t13 = -t52 + t53;
t30 = sin(qJ(3));
t33 = cos(qJ(3));
t5 = g(3) * t33 + t13 * t30;
t56 = -pkin(1) - pkin(7);
t55 = -pkin(3) - pkin(8);
t51 = g(3) * t30;
t49 = t30 * pkin(3);
t48 = t30 * t34;
t47 = t31 * t33;
t46 = t34 * t29;
t45 = t34 * t32;
t42 = qJ(4) * t30;
t44 = pkin(3) * t47 + t31 * t42;
t43 = t34 * pkin(1) + t31 * qJ(2);
t23 = t33 * qJ(4);
t40 = t34 * pkin(7) + t31 * t49 + t43;
t7 = t31 * t29 - t33 * t45;
t9 = t32 * t47 + t46;
t39 = g(1) * t7 - g(2) * t9;
t14 = g(1) * t34 + g(2) * t31;
t24 = t34 * qJ(2);
t37 = pkin(3) * t48 - t34 * t23 + t24;
t36 = -g(1) * t9 - g(2) * t7 + t32 * t51;
t10 = -t29 * t47 + t45;
t8 = t31 * t32 + t33 * t46;
t35 = g(1) * t10 + g(2) * t8 + t29 * t51;
t12 = t14 * t33;
t11 = t14 * t30;
t6 = g(1) * t47 - t33 * t52 - t51;
t4 = t5 * t32;
t3 = t5 * t29;
t2 = g(1) * t8 - g(2) * t10;
t1 = [0, t13, t14, -t13, -t14, -g(1) * (-t31 * pkin(1) + t24) - g(2) * t43, 0, 0, 0, 0, 0, -t11, -t12, t13, t11, t12, -g(1) * (t56 * t31 + t37) - g(2) * (-t31 * t23 + t40) 0, 0, 0, 0, 0, t2, -t39, t2, -t11, t39, -g(1) * (-t8 * pkin(5) + pkin(8) * t48 - t7 * qJ(6) + t37) - g(2) * (t34 * pkin(4) + t10 * pkin(5) + t9 * qJ(6) + t40) + (-g(1) * (-pkin(4) + t56) - g(2) * (t30 * pkin(8) - t23)) * t31; 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t6, -t5, -g(1) * t44 - g(3) * (t23 - t49) - (-pkin(3) * t33 - t42) * t52, 0, 0, 0, 0, 0, -t3, -t4, -t3, -t6, t4, -g(1) * (pkin(8) * t47 + t44) - g(3) * (t38 * t33 + t23) + (-g(3) * t55 - t38 * t53) * t30 - (t55 * t33 + (-qJ(4) - t38) * t30) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35, -t36, 0, -t35, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (pkin(5) * t32 + qJ(6) * t29) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36;];
taug_reg  = t1;
