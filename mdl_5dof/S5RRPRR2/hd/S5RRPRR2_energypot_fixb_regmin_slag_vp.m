% Calculate minimal parameter regressor of potential energy for
% S5RRPRR2
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
% Datum: 2021-01-15 21:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:23:34
% EndTime: 2021-01-15 21:23:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->22), mult. (56->29), div. (0->0), fcn. (53->10), ass. (0->18)
t77 = qJ(2) + pkin(9);
t76 = qJ(4) + t77;
t80 = sin(qJ(1));
t82 = cos(qJ(1));
t83 = g(1) * t82 + g(2) * t80;
t81 = cos(qJ(2));
t79 = sin(qJ(2));
t78 = pkin(6) + qJ(3);
t75 = cos(t77);
t74 = sin(t77);
t73 = t81 * pkin(2) + pkin(1);
t72 = qJ(5) + t76;
t71 = cos(t76);
t70 = sin(t76);
t69 = cos(t72);
t68 = sin(t72);
t67 = g(1) * t80 - g(2) * t82;
t1 = [0, -t83, t67, 0, 0, 0, 0, 0, -g(3) * t79 - t83 * t81, -g(3) * t81 + t83 * t79, -g(3) * t74 - t83 * t75, -g(3) * t75 + t83 * t74, -t67, -g(1) * (t82 * t73 + t78 * t80) - g(2) * (t80 * t73 - t82 * t78) - g(3) * (t79 * pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t70 - t83 * t71, -g(3) * t71 + t83 * t70, 0, 0, 0, 0, 0, -g(3) * t68 - t83 * t69, -g(3) * t69 + t83 * t68;];
U_reg = t1;
