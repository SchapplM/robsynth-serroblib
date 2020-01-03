% Calculate inertial parameters regressor of potential energy for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:54
% EndTime: 2019-12-31 16:25:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (53->41), mult. (102->51), div. (0->0), fcn. (96->6), ass. (0->24)
t33 = sin(pkin(6));
t36 = sin(qJ(2));
t43 = qJ(3) * t36;
t38 = cos(qJ(2));
t48 = t33 * t38;
t51 = pkin(2) * t48 + t33 * t43;
t50 = g(3) * t38;
t49 = g(3) * qJ(1);
t34 = cos(pkin(6));
t47 = t34 * t38;
t35 = sin(qJ(4));
t46 = t35 * t36;
t37 = cos(qJ(4));
t45 = t36 * t37;
t44 = t34 * pkin(1) + t33 * pkin(4);
t29 = t33 * pkin(1);
t42 = -t34 * pkin(4) + t29;
t41 = pkin(2) * t47 + t34 * t43 + t44;
t40 = g(1) * t34 + g(2) * t33;
t39 = t36 * pkin(2) - t38 * qJ(3) + qJ(1);
t23 = g(1) * t33 - g(2) * t34;
t22 = g(3) * t36 + t40 * t38;
t21 = t40 * t36 - t50;
t1 = [0, 0, 0, 0, 0, 0, -t40, t23, -g(3), -t49, 0, 0, 0, 0, 0, 0, -t22, t21, -t23, -g(1) * t44 - g(2) * t42 - t49, 0, 0, 0, 0, 0, 0, -t23, t22, -t21, -g(1) * t41 - g(2) * (t42 + t51) - g(3) * t39, 0, 0, 0, 0, 0, 0, -g(1) * (t33 * t37 + t34 * t46) - g(2) * (t33 * t46 - t34 * t37) + t35 * t50, -g(1) * (-t33 * t35 + t34 * t45) - g(2) * (t33 * t45 + t34 * t35) + t37 * t50, -t22, -g(1) * (t33 * pkin(3) + pkin(5) * t47 + t41) - g(2) * (pkin(5) * t48 + t29 + (-pkin(3) - pkin(4)) * t34 + t51) - g(3) * (t36 * pkin(5) + t39);];
U_reg = t1;
