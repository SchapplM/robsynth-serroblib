% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:11:21
% EndTime: 2019-05-05 18:11:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (162->87), mult. (383->97), div. (0->0), fcn. (387->6), ass. (0->48)
t33 = sin(qJ(5));
t36 = cos(qJ(5));
t45 = pkin(5) * t33 - qJ(6) * t36;
t38 = cos(qJ(1));
t63 = g(2) * t38;
t35 = sin(qJ(1));
t64 = g(1) * t35;
t13 = -t63 + t64;
t34 = sin(qJ(3));
t37 = cos(qJ(3));
t5 = g(3) * t37 + t13 * t34;
t67 = -pkin(1) - pkin(7);
t66 = -pkin(3) - pkin(8);
t62 = g(3) * t34;
t60 = t34 * t35;
t59 = t34 * t38;
t58 = t35 * t37;
t57 = t38 * t33;
t56 = t38 * t36;
t53 = qJ(4) * t34;
t55 = pkin(3) * t58 + t35 * t53;
t54 = t38 * pkin(1) + t35 * qJ(2);
t26 = t37 * qJ(4);
t51 = t38 * pkin(7) + t54;
t50 = t67 * t35;
t49 = t66 * t37;
t48 = pkin(3) * t60 + t51;
t47 = g(1) * (pkin(8) * t58 + t55);
t7 = t35 * t33 - t37 * t56;
t9 = t36 * t58 + t57;
t46 = g(1) * t7 - g(2) * t9;
t14 = g(1) * t38 + g(2) * t35;
t27 = t38 * qJ(2);
t44 = pkin(3) * t59 - t38 * t26 + t27;
t43 = t38 * pkin(4) + pkin(8) * t60 + t48;
t42 = pkin(8) * t59 + t44;
t41 = -g(1) * t9 - g(2) * t7 + t36 * t62;
t10 = -t33 * t58 + t56;
t8 = t35 * t36 + t37 * t57;
t40 = g(1) * t10 + g(2) * t8 + t33 * t62;
t39 = (-g(1) * (-pkin(4) + t67) + g(2) * t26) * t35;
t12 = t14 * t37;
t11 = t14 * t34;
t6 = g(1) * t58 - t37 * t63 - t62;
t4 = t5 * t36;
t3 = t5 * t33;
t2 = g(1) * t8 - g(2) * t10;
t1 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -g(1) * (-t35 * pkin(1) + t27) - g(2) * t54, 0, 0, 0, 0, 0, 0, -t11, -t12, t13, -g(1) * (t27 + t50) - g(2) * t51, 0, 0, 0, 0, 0, 0, t13, t11, t12, -g(1) * (t50 + t44) - g(2) * (-t35 * t26 + t48) 0, 0, 0, 0, 0, 0, t2, -t46, -t11, -g(1) * t42 - g(2) * t43 + t39, 0, 0, 0, 0, 0, 0, t2, -t11, t46, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t42) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t43) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -g(1) * t55 - g(3) * (-t34 * pkin(3) + t26) - (-pkin(3) * t37 - t53) * t63, 0, 0, 0, 0, 0, 0, -t3, -t4, -t6, -t47 - g(3) * (t66 * t34 + t26) - (t49 - t53) * t63, 0, 0, 0, 0, 0, 0, -t3, -t6, t4, -t47 - g(3) * (t45 * t37 + t26) + (-g(3) * t66 - t45 * t64) * t34 - (t49 + (-qJ(4) - t45) * t34) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t40, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, -t40, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (pkin(5) * t36 + qJ(6) * t33) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41;];
taug_reg  = t1;
