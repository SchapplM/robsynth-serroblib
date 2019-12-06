% Calculate inertial parameters regressor of potential energy for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:03
% EndTime: 2019-12-05 15:31:03
% DurationCPUTime: 0.11s
% Computational Cost: add. (128->49), mult. (130->62), div. (0->0), fcn. (125->8), ass. (0->29)
t62 = cos(qJ(4));
t48 = t62 * pkin(4) + pkin(3);
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t59 = -qJ(5) - pkin(6);
t76 = t48 * t57 - t55 * t59;
t75 = g(3) * t55;
t60 = pkin(5) + qJ(1);
t74 = g(3) * t60;
t54 = pkin(7) + qJ(2);
t49 = sin(t54);
t61 = sin(qJ(4));
t72 = t49 * t61;
t70 = t57 * t61;
t69 = t57 * t62;
t56 = sin(pkin(7));
t68 = t56 * pkin(1) + t49 * pkin(2);
t50 = cos(t54);
t58 = cos(pkin(7));
t67 = t58 * pkin(1) + t50 * pkin(2) + t49 * qJ(3);
t66 = -t50 * qJ(3) + t68;
t65 = pkin(3) * t57 + pkin(6) * t55;
t64 = g(1) * t50 + g(2) * t49;
t63 = -g(1) * t58 - g(2) * t56;
t44 = g(1) * t49 - g(2) * t50;
t43 = -g(3) * t57 + t64 * t55;
t42 = -g(1) * (t50 * t69 + t72) - g(2) * (t49 * t69 - t50 * t61) - t62 * t75;
t41 = -g(1) * (t49 * t62 - t50 * t70) - g(2) * (-t49 * t70 - t50 * t62) + t61 * t75;
t1 = [0, 0, 0, 0, 0, 0, t63, g(1) * t56 - g(2) * t58, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t64, t44, -g(3), t63 * pkin(1) - t74, 0, 0, 0, 0, 0, 0, -t64 * t57 - t75, t43, -t44, -g(1) * t67 - g(2) * t66 - t74, 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * (t65 * t50 + t67) - g(2) * (t65 * t49 + t66) - g(3) * (t55 * pkin(3) - t57 * pkin(6) + t60), 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * (pkin(4) * t72 + t67) - g(2) * (t76 * t49 + t68) - g(3) * (t55 * t48 + t57 * t59 + t60) + (-g(1) * t76 - g(2) * (-pkin(4) * t61 - qJ(3))) * t50;];
U_reg = t1;
