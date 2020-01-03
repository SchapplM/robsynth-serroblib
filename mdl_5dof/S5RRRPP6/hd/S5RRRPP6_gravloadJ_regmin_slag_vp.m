% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = cos(qJ(3));
t20 = t32 * pkin(3) + pkin(2);
t33 = cos(qJ(2));
t13 = t33 * t20;
t28 = -qJ(4) - pkin(7);
t30 = sin(qJ(2));
t57 = -t30 * t28 + t13;
t56 = -pkin(1) - t57;
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t42 = g(1) * t34 + g(2) * t31;
t5 = -g(3) * t33 + t42 * t30;
t53 = g(3) * t30;
t29 = sin(qJ(3));
t50 = t31 * t29;
t49 = t31 * t32;
t48 = t31 * t33;
t27 = qJ(3) + pkin(8);
t21 = sin(t27);
t47 = t34 * t21;
t22 = cos(t27);
t46 = t34 * t22;
t45 = t34 * t29;
t44 = t34 * t32;
t43 = t33 * t45;
t41 = g(1) * t31 - g(2) * t34;
t40 = pkin(4) * t22 + qJ(5) * t21;
t7 = t29 * t48 + t44;
t37 = pkin(3) * t50 + t31 * pkin(6) - t56 * t34;
t1 = t21 * t48 + t46;
t3 = -t31 * t22 + t33 * t47;
t36 = g(1) * t3 + g(2) * t1 + t21 * t53;
t35 = pkin(3) * t45 + t34 * pkin(6) + t56 * t31;
t6 = t42 * t33 + t53;
t18 = pkin(3) * t49;
t11 = t41 * t30;
t10 = t33 * t44 + t50;
t9 = -t43 + t49;
t8 = -t32 * t48 + t45;
t4 = t31 * t21 + t33 * t46;
t2 = t22 * t48 - t47;
t12 = [0, t41, t42, 0, 0, 0, 0, 0, t41 * t33, -t11, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11, -g(1) * t35 - g(2) * t37, g(1) * t2 - g(2) * t4, t11, g(1) * t1 - g(2) * t3, -g(1) * (-t2 * pkin(4) - t1 * qJ(5) + t35) - g(2) * (t4 * pkin(4) + t3 * qJ(5) + t37); 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t5 * t32, -t5 * t29, -t6, -g(3) * t57 + t42 * (t20 * t30 + t28 * t33), t5 * t22, -t6, t5 * t21, -g(3) * t13 + (-g(3) * t40 + t42 * t28) * t33 + (g(3) * t28 + t42 * (t20 + t40)) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t29 * t53, g(1) * t10 - g(2) * t8 + t32 * t53, 0, -g(1) * t18 + (g(2) * t44 + t6 * t29) * pkin(3), t36, 0, -g(1) * t4 - g(2) * t2 - t22 * t53, -g(1) * (-pkin(3) * t43 - t3 * pkin(4) + t4 * qJ(5) + t18) - g(2) * (-t7 * pkin(3) - t1 * pkin(4) + t2 * qJ(5)) - (-pkin(3) * t29 - pkin(4) * t21 + qJ(5) * t22) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36;];
taug_reg = t12;
