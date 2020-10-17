% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:53
% EndTime: 2019-12-31 21:05:54
% DurationCPUTime: 0.24s
% Computational Cost: add. (165->69), mult. (420->94), div. (0->0), fcn. (452->6), ass. (0->46)
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t45 = g(1) * t39 + g(2) * t36;
t65 = t45 * t35;
t64 = -pkin(3) - pkin(4);
t63 = g(1) * t36;
t29 = t35 * pkin(7);
t38 = cos(qJ(2));
t31 = t38 * pkin(2);
t34 = sin(qJ(3));
t60 = t34 * t35;
t37 = cos(qJ(3));
t59 = t35 * t37;
t58 = t35 * t39;
t57 = t36 * t38;
t56 = t37 * t38;
t55 = t38 * t39;
t54 = t39 * t34;
t53 = qJ(4) * t34;
t52 = qJ(5) * t39;
t51 = -pkin(1) - t31;
t50 = -pkin(2) - t53;
t14 = t34 * t57 + t37 * t39;
t15 = t36 * t56 - t54;
t49 = -t14 * pkin(3) + qJ(4) * t15;
t16 = -t36 * t37 + t38 * t54;
t17 = t36 * t34 + t37 * t55;
t48 = -t16 * pkin(3) + qJ(4) * t17;
t47 = pkin(3) * t56 + t38 * t53 + t29 + t31;
t46 = -t15 * pkin(3) + t39 * pkin(6) - t14 * qJ(4);
t4 = g(1) * t14 - g(2) * t16;
t44 = -g(2) * t39 + t63;
t41 = t39 * pkin(1) + pkin(2) * t55 + t17 * pkin(3) + t36 * pkin(6) + pkin(7) * t58 + t16 * qJ(4);
t2 = g(1) * t16 + g(2) * t14 + g(3) * t60;
t40 = g(1) * t17 + g(2) * t15 + g(3) * t59;
t8 = -g(3) * t38 + t65;
t24 = pkin(7) * t55;
t21 = pkin(7) * t57;
t19 = qJ(4) * t59;
t18 = -g(2) * t58 + t35 * t63;
t9 = g(3) * t35 + t45 * t38;
t7 = t8 * t37;
t6 = t8 * t34;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, t44, t45, 0, 0, 0, 0, 0, t44 * t38, -t18, 0, 0, 0, 0, 0, t5, -t4, t5, t18, t4, -g(1) * t46 - g(2) * t41 - (t51 - t29) * t63, t5, t4, -t18, -g(1) * (-t15 * pkin(4) + t46) - g(2) * (t17 * pkin(4) - t35 * t52 + t41) - ((-pkin(7) + qJ(5)) * t35 + t51) * t63; 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, t7, -t6, t7, -t9, t6, -g(1) * t24 - g(2) * t21 - g(3) * t47 + (pkin(3) * t37 - t50) * t65, t7, t6, t9, -g(1) * (-t38 * t52 + t24) - g(2) * (-qJ(5) * t57 + t21) - g(3) * (pkin(4) * t56 + t47) + (g(3) * qJ(5) + t45 * (-t64 * t37 - t50)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t40, t2, 0, -t40, -g(1) * t48 - g(2) * t49 - g(3) * (-pkin(3) * t60 + t19), t2, -t40, 0, -g(1) * (-pkin(4) * t16 + t48) - g(2) * (-pkin(4) * t14 + t49) - g(3) * (t64 * t60 + t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg = t1;
