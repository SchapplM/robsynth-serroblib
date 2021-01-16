% Calculate minimal parameter regressor of potential energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:19
% EndTime: 2021-01-16 00:21:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->36), mult. (108->52), div. (0->0), fcn. (117->8), ass. (0->24)
t77 = sin(qJ(2));
t90 = g(3) * t77;
t75 = qJ(3) + qJ(4);
t72 = sin(t75);
t76 = sin(qJ(3));
t89 = t76 * pkin(3) + pkin(4) * t72 + pkin(6);
t78 = sin(qJ(1));
t80 = cos(qJ(2));
t88 = t78 * t80;
t81 = cos(qJ(1));
t87 = t81 * t72;
t73 = cos(t75);
t86 = t81 * t73;
t85 = t81 * t76;
t79 = cos(qJ(3));
t84 = t81 * t79;
t83 = g(1) * t81 + g(2) * t78;
t70 = t79 * pkin(3) + pkin(4) * t73 + pkin(2);
t74 = -qJ(5) - pkin(8) - pkin(7);
t82 = t70 * t80 - t74 * t77 + pkin(1);
t69 = -g(3) * t80 + t83 * t77;
t68 = -g(1) * (t78 * t72 + t80 * t86) - g(2) * (t73 * t88 - t87) - t73 * t90;
t67 = -g(1) * (t78 * t73 - t80 * t87) - g(2) * (-t72 * t88 - t86) + t72 * t90;
t1 = [0, -t83, g(1) * t78 - g(2) * t81, 0, 0, 0, 0, 0, -t83 * t80 - t90, t69, 0, 0, 0, 0, 0, -g(1) * (t78 * t76 + t80 * t84) - g(2) * (t79 * t88 - t85) - t79 * t90, -g(1) * (t78 * t79 - t80 * t85) - g(2) * (-t76 * t88 - t84) + t76 * t90, 0, 0, 0, 0, 0, t68, t67, t68, t67, -t69, -g(3) * (t77 * t70 + t80 * t74 + pkin(5)) + (-g(1) * t82 + g(2) * t89) * t81 + (-g(1) * t89 - g(2) * t82) * t78;];
U_reg = t1;
