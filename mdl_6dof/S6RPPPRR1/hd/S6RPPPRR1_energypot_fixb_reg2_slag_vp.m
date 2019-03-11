% Calculate inertial parameters regressor of potential energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:13
% EndTime: 2019-03-09 01:30:13
% DurationCPUTime: 0.09s
% Computational Cost: add. (138->50), mult. (115->58), div. (0->0), fcn. (103->8), ass. (0->28)
t56 = qJ(2) + pkin(6);
t76 = g(3) * t56;
t61 = cos(qJ(5));
t75 = g(3) * t61;
t57 = sin(qJ(6));
t58 = sin(qJ(5));
t74 = t57 * t58;
t60 = cos(qJ(6));
t73 = t58 * t60;
t72 = pkin(3) + t56;
t55 = qJ(1) + pkin(9);
t51 = sin(t55);
t52 = cos(t55);
t62 = cos(qJ(1));
t71 = pkin(1) * t62 + pkin(2) * t52 + qJ(3) * t51;
t70 = pkin(4) + t72;
t69 = qJ(4) * t52 + t71;
t59 = sin(qJ(1));
t68 = pkin(1) * t59 + pkin(2) * t51 - t52 * qJ(3);
t67 = pkin(5) * t58 - pkin(8) * t61;
t43 = g(1) * t52 + g(2) * t51;
t66 = -g(1) * t62 - g(2) * t59;
t65 = qJ(4) * t51 + t68;
t64 = -t51 * pkin(7) + t69;
t63 = pkin(7) * t52 + t65;
t42 = g(1) * t51 - g(2) * t52;
t41 = -g(3) * t58 + t43 * t61;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t59 - g(2) * t62, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t43, t42, -g(3), pkin(1) * t66 - t76, 0, 0, 0, 0, 0, 0, -g(3), t43, -t42, -g(1) * t71 - g(2) * t68 - t76, 0, 0, 0, 0, 0, 0, -g(3), -t42, -t43, -g(1) * t69 - g(2) * t65 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t43 * t58 - t75, -t41, t42, -g(1) * t64 - g(2) * t63 - g(3) * t70, 0, 0, 0, 0, 0, 0, -g(1) * (-t51 * t57 + t52 * t73) - g(2) * (t51 * t73 + t52 * t57) - t60 * t75, -g(1) * (-t51 * t60 - t52 * t74) - g(2) * (-t51 * t74 + t52 * t60) + t57 * t75, t41, -g(1) * (t52 * t67 + t64) - g(2) * (t51 * t67 + t63) - g(3) * (pkin(5) * t61 + pkin(8) * t58 + t70);];
U_reg  = t1;
