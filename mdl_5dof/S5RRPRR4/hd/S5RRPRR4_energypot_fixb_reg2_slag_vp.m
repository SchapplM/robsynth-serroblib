% Calculate inertial parameters regressor of potential energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:05
% EndTime: 2020-01-03 12:02:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->41), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t67 = pkin(6) + pkin(5);
t54 = qJ(3) + t67;
t66 = g(1) * t54;
t56 = qJ(1) + qJ(2);
t50 = sin(t56);
t58 = sin(qJ(1));
t65 = t58 * pkin(1) + pkin(2) * t50;
t52 = cos(t56);
t60 = cos(qJ(1));
t64 = -t60 * pkin(1) - pkin(2) * t52;
t48 = pkin(9) + t56;
t44 = sin(t48);
t45 = cos(t48);
t63 = g(2) * t44 - g(3) * t45;
t62 = -g(2) * t58 + g(3) * t60;
t61 = -pkin(8) - pkin(7);
t59 = cos(qJ(4));
t57 = sin(qJ(4));
t55 = qJ(4) + qJ(5);
t51 = cos(t55);
t49 = sin(t55);
t47 = t59 * pkin(4) + pkin(3);
t43 = g(2) * t45 + g(3) * t44;
t1 = [0, 0, 0, 0, 0, 0, t62, -g(2) * t60 - g(3) * t58, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t50 + g(3) * t52, -g(2) * t52 - g(3) * t50, -g(1), t62 * pkin(1) - g(1) * t67, 0, 0, 0, 0, 0, 0, -t63, -t43, -g(1), -g(2) * t65 - g(3) * t64 - t66, 0, 0, 0, 0, 0, 0, -g(1) * t57 - t63 * t59, -g(1) * t59 + t63 * t57, t43, -t66 - g(2) * (t44 * pkin(3) - t45 * pkin(7) + t65) - g(3) * (-t45 * pkin(3) - t44 * pkin(7) + t64), 0, 0, 0, 0, 0, 0, -g(1) * t49 - t63 * t51, -g(1) * t51 + t63 * t49, t43, -g(1) * (t57 * pkin(4) + t54) - g(2) * (t44 * t47 + t45 * t61 + t65) - g(3) * (t44 * t61 - t45 * t47 + t64);];
U_reg = t1;
