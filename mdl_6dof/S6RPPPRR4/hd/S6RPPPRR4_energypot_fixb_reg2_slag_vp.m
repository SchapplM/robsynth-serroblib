% Calculate inertial parameters regressor of potential energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:51
% EndTime: 2019-03-09 01:35:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (126->53), mult. (211->68), div. (0->0), fcn. (243->8), ass. (0->32)
t83 = g(3) * pkin(6);
t59 = -qJ(3) + pkin(6);
t82 = g(3) * t59;
t63 = cos(qJ(5));
t81 = g(3) * t63;
t80 = cos(qJ(1));
t79 = sin(qJ(1));
t60 = sin(qJ(6));
t61 = sin(qJ(5));
t78 = t60 * t61;
t62 = cos(qJ(6));
t77 = t61 * t62;
t76 = t80 * pkin(1) + t79 * qJ(2);
t75 = cos(pkin(9));
t74 = sin(pkin(9));
t73 = -pkin(4) + t59;
t72 = t80 * pkin(2) + t76;
t47 = -t79 * t74 - t80 * t75;
t71 = -t47 * pkin(3) + t72;
t48 = t80 * t74 - t79 * t75;
t70 = g(1) * t48 - g(2) * t47;
t44 = g(1) * t47 + g(2) * t48;
t69 = t79 * pkin(1) - t80 * qJ(2);
t68 = pkin(5) * t61 - pkin(8) * t63 + qJ(4);
t67 = t79 * pkin(2) + t69;
t66 = t48 * qJ(4) + t71;
t65 = -t48 * pkin(3) + t67;
t64 = -t47 * qJ(4) + t65;
t50 = -g(1) * t80 - g(2) * t79;
t49 = g(1) * t79 - g(2) * t80;
t43 = g(3) * t61 + t70 * t63;
t1 = [0, 0, 0, 0, 0, 0, t50, t49, -g(3), -t83, 0, 0, 0, 0, 0, 0, t50, -g(3), -t49, -g(1) * t76 - g(2) * t69 - t83, 0, 0, 0, 0, 0, 0, t44, t70, g(3), -g(1) * t72 - g(2) * t67 - t82, 0, 0, 0, 0, 0, 0, g(3), -t44, -t70, -g(1) * t66 - g(2) * t64 - t82, 0, 0, 0, 0, 0, 0, -t70 * t61 + t81, -t43, t44, -g(1) * (-t47 * pkin(7) + t66) - g(2) * (-t48 * pkin(7) + t64) - g(3) * t73, 0, 0, 0, 0, 0, 0, -g(1) * (-t47 * t60 + t48 * t77) - g(2) * (-t47 * t77 - t48 * t60) + t62 * t81, -g(1) * (-t47 * t62 - t48 * t78) - g(2) * (t47 * t78 - t48 * t62) - t60 * t81, t43, -g(1) * t71 - g(2) * t65 - g(3) * (-t63 * pkin(5) - t61 * pkin(8) + t73) + (g(2) * pkin(7) - g(1) * t68) * t48 + (g(1) * pkin(7) + g(2) * t68) * t47;];
U_reg  = t1;
