% Calculate minimal parameter regressor of potential energy for
% S5RRPRP7
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
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:41:22
% EndTime: 2021-01-15 20:41:22
% DurationCPUTime: 0.07s
% Computational Cost: add. (96->43), mult. (123->54), div. (0->0), fcn. (129->10), ass. (0->32)
t78 = qJ(2) + pkin(8);
t75 = sin(t78);
t95 = g(3) * t75;
t83 = sin(qJ(2));
t94 = t83 * pkin(2) + pkin(5);
t82 = sin(qJ(4));
t84 = sin(qJ(1));
t93 = t84 * t82;
t85 = cos(qJ(4));
t92 = t84 * t85;
t87 = cos(qJ(1));
t91 = t87 * t82;
t90 = t87 * t85;
t89 = g(1) * t87 + g(2) * t84;
t76 = cos(t78);
t66 = t76 * t93 + t90;
t68 = t76 * t91 - t92;
t88 = g(1) * t68 + g(2) * t66 + t82 * t95;
t86 = cos(qJ(2));
t81 = pkin(6) + qJ(3);
t80 = cos(pkin(8));
t79 = sin(pkin(8));
t74 = t86 * pkin(2) + pkin(1);
t72 = t81 * t84;
t71 = t87 * t81;
t70 = g(1) * t84 - g(2) * t87;
t69 = t76 * t90 + t93;
t67 = t76 * t92 - t91;
t65 = -g(3) * t76 + t89 * t75;
t64 = (t80 * pkin(3) + t79 * pkin(7) + pkin(2)) * t86 + (-t79 * pkin(3) + t80 * pkin(7)) * t83 + pkin(1);
t63 = -g(1) * t69 - g(2) * t67 - t85 * t95;
t1 = [0, -t89, t70, 0, 0, 0, 0, 0, -g(3) * t83 - t89 * t86, -g(3) * t86 + t89 * t83, -t89 * t76 - t95, t65, -t70, -g(1) * (t87 * t74 + t72) - g(2) * (t84 * t74 - t71) - g(3) * t94, 0, 0, 0, 0, 0, t63, t88, t63, -t65, -t88, -g(1) * (t69 * pkin(4) + t68 * qJ(5) + t64 * t87 + t72) - g(2) * (t67 * pkin(4) + t66 * qJ(5) + t64 * t84 - t71) - g(3) * (-t76 * pkin(7) + t94) - (pkin(4) * t85 + qJ(5) * t82 + pkin(3)) * t95;];
U_reg = t1;
