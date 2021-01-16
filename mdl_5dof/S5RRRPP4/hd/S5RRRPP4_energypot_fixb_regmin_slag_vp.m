% Calculate minimal parameter regressor of potential energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:47
% EndTime: 2021-01-15 22:24:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (102->30), mult. (84->36), div. (0->0), fcn. (78->8), ass. (0->21)
t75 = qJ(2) + qJ(3);
t71 = cos(t75);
t78 = cos(qJ(2));
t63 = t78 * pkin(2) + pkin(3) * t71 + pkin(1);
t74 = -qJ(4) - pkin(7) - pkin(6);
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t84 = t77 * t63 + t79 * t74;
t70 = sin(t75);
t76 = sin(qJ(2));
t83 = t76 * pkin(2) + pkin(3) * t70 + pkin(5);
t82 = t79 * t63 - t77 * t74;
t81 = g(1) * t79 + g(2) * t77;
t69 = pkin(8) + t75;
t66 = sin(t69);
t67 = cos(t69);
t80 = pkin(4) * t67 + qJ(5) * t66;
t64 = g(1) * t77 - g(2) * t79;
t60 = -g(3) * t66 - t81 * t67;
t59 = -g(3) * t67 + t81 * t66;
t1 = [0, -t81, t64, 0, 0, 0, 0, 0, -g(3) * t76 - t81 * t78, -g(3) * t78 + t81 * t76, 0, 0, 0, 0, 0, -g(3) * t70 - t81 * t71, -g(3) * t71 + t81 * t70, t60, t59, -t64, -g(1) * t82 - g(2) * t84 - g(3) * t83, t60, -t64, -t59, -g(1) * (t80 * t79 + t82) - g(2) * (t80 * t77 + t84) - g(3) * (t66 * pkin(4) - t67 * qJ(5) + t83);];
U_reg = t1;
