% Calculate minimal parameter regressor of potential energy for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:29:56
% EndTime: 2021-01-15 20:29:56
% DurationCPUTime: 0.06s
% Computational Cost: add. (88->36), mult. (108->47), div. (0->0), fcn. (110->8), ass. (0->28)
t77 = qJ(2) + pkin(8);
t74 = sin(t77);
t95 = g(3) * t74;
t81 = sin(qJ(2));
t94 = t81 * pkin(2) + pkin(5);
t80 = sin(qJ(4));
t82 = sin(qJ(1));
t93 = t82 * t80;
t83 = cos(qJ(4));
t92 = t82 * t83;
t85 = cos(qJ(1));
t91 = t85 * t80;
t90 = t85 * t83;
t84 = cos(qJ(2));
t73 = t84 * pkin(2) + pkin(1);
t79 = pkin(6) + qJ(3);
t89 = t82 * t73 - t85 * t79;
t88 = t85 * t73 + t79 * t82;
t87 = g(1) * t85 + g(2) * t82;
t72 = t83 * pkin(4) + pkin(3);
t75 = cos(t77);
t78 = -qJ(5) - pkin(7);
t86 = t72 * t75 - t74 * t78;
t67 = g(1) * t82 - g(2) * t85;
t66 = -g(3) * t75 + t87 * t74;
t65 = -g(1) * (t75 * t90 + t93) - g(2) * (t75 * t92 - t91) - t83 * t95;
t64 = -g(1) * (-t75 * t91 + t92) - g(2) * (-t75 * t93 - t90) + t80 * t95;
t1 = [0, -t87, t67, 0, 0, 0, 0, 0, -g(3) * t81 - t87 * t84, -g(3) * t84 + t87 * t81, -t87 * t75 - t95, t66, -t67, -g(1) * t88 - g(2) * t89 - g(3) * t94, 0, 0, 0, 0, 0, t65, t64, t65, t64, -t66, -g(1) * (pkin(4) * t93 + t86 * t85 + t88) - g(2) * (-pkin(4) * t91 + t86 * t82 + t89) - g(3) * (t74 * t72 + t75 * t78 + t94);];
U_reg = t1;
