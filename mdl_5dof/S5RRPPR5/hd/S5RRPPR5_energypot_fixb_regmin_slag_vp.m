% Calculate minimal parameter regressor of potential energy for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:29
% EndTime: 2021-01-15 19:36:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (80->33), mult. (92->42), div. (0->0), fcn. (92->10), ass. (0->24)
t78 = sin(qJ(2));
t85 = t78 * pkin(2) + pkin(5);
t84 = cos(qJ(5));
t83 = sin(qJ(5));
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t82 = g(1) * t81 + g(2) * t79;
t80 = cos(qJ(2));
t77 = pkin(6) + qJ(3);
t76 = cos(pkin(8));
t75 = sin(pkin(8));
t74 = qJ(2) + pkin(8);
t71 = cos(t74);
t70 = sin(t74);
t69 = t80 * pkin(2) + pkin(1);
t68 = t77 * t79;
t67 = t81 * t77;
t64 = g(1) * t79 - g(2) * t81;
t63 = t70 * t84 - t71 * t83;
t62 = -t70 * t83 - t71 * t84;
t61 = -g(3) * t70 - t82 * t71;
t60 = -g(3) * t71 + t82 * t70;
t59 = (pkin(3) * t76 + qJ(4) * t75 + pkin(2)) * t80 + (-t75 * pkin(3) + qJ(4) * t76) * t78 + pkin(1);
t1 = [0, -t82, t64, 0, 0, 0, 0, 0, -g(3) * t78 - t82 * t80, -g(3) * t80 + t82 * t78, t61, t60, -t64, -g(1) * (t81 * t69 + t68) - g(2) * (t79 * t69 - t67) - g(3) * t85, t61, -t64, -t60, -g(1) * (t59 * t81 + t68) - g(2) * (t59 * t79 - t67) - g(3) * (t70 * pkin(3) - t71 * qJ(4) + t85), 0, 0, 0, 0, 0, -g(3) * t63 + t82 * t62, -g(3) * t62 - t82 * t63;];
U_reg = t1;
