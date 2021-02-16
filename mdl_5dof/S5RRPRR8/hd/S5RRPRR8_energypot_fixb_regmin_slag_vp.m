% Calculate minimal parameter regressor of potential energy for
% S5RRPRR8
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
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:10
% EndTime: 2021-01-15 21:35:10
% DurationCPUTime: 0.05s
% Computational Cost: add. (62->27), mult. (66->39), div. (0->0), fcn. (67->10), ass. (0->22)
t73 = qJ(2) + pkin(9);
t72 = qJ(4) + t73;
t67 = sin(t72);
t86 = g(3) * t67;
t75 = sin(qJ(5));
t77 = sin(qJ(1));
t85 = t77 * t75;
t78 = cos(qJ(5));
t84 = t77 * t78;
t80 = cos(qJ(1));
t83 = t80 * t75;
t82 = t80 * t78;
t81 = g(1) * t80 + g(2) * t77;
t79 = cos(qJ(2));
t76 = sin(qJ(2));
t74 = pkin(6) + qJ(3);
t71 = cos(t73);
t70 = sin(t73);
t69 = t79 * pkin(2) + pkin(1);
t68 = cos(t72);
t66 = g(1) * t77 - g(2) * t80;
t1 = [0, -t81, t66, 0, 0, 0, 0, 0, -g(3) * t76 - t81 * t79, -g(3) * t79 + t81 * t76, -g(3) * t70 - t81 * t71, -g(3) * t71 + t81 * t70, -t66, -g(1) * (t80 * t69 + t74 * t77) - g(2) * (t77 * t69 - t80 * t74) - g(3) * (t76 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -t81 * t68 - t86, -g(3) * t68 + t81 * t67, 0, 0, 0, 0, 0, -g(1) * (t68 * t82 + t85) - g(2) * (t68 * t84 - t83) - t78 * t86, -g(1) * (-t68 * t83 + t84) - g(2) * (-t68 * t85 - t82) + t75 * t86;];
U_reg = t1;
