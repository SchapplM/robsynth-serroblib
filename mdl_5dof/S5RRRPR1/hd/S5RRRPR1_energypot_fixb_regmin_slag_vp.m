% Calculate minimal parameter regressor of potential energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:48:56
% EndTime: 2021-01-15 22:48:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (72->24), mult. (59->31), div. (0->0), fcn. (56->10), ass. (0->18)
t78 = qJ(2) + qJ(3);
t74 = pkin(9) + t78;
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t83 = g(1) * t82 + g(2) * t80;
t81 = cos(qJ(2));
t79 = sin(qJ(2));
t77 = -qJ(4) - pkin(7) - pkin(6);
t76 = cos(t78);
t75 = sin(t78);
t73 = qJ(5) + t74;
t72 = cos(t74);
t71 = sin(t74);
t70 = cos(t73);
t69 = sin(t73);
t68 = g(1) * t80 - g(2) * t82;
t67 = t81 * pkin(2) + pkin(3) * t76 + pkin(1);
t1 = [0, -t83, t68, 0, 0, 0, 0, 0, -g(3) * t79 - t83 * t81, -g(3) * t81 + t83 * t79, 0, 0, 0, 0, 0, -g(3) * t75 - t83 * t76, -g(3) * t76 + t83 * t75, -g(3) * t71 - t83 * t72, -g(3) * t72 + t83 * t71, -t68, -g(1) * (t82 * t67 - t80 * t77) - g(2) * (t80 * t67 + t82 * t77) - g(3) * (t79 * pkin(2) + pkin(3) * t75 + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t69 - t83 * t70, -g(3) * t70 + t83 * t69;];
U_reg = t1;
