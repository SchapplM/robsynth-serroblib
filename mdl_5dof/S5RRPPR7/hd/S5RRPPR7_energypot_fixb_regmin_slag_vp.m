% Calculate minimal parameter regressor of potential energy for
% S5RRPPR7
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
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:00
% EndTime: 2021-01-15 19:59:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (72->38), mult. (90->48), div. (0->0), fcn. (88->10), ass. (0->27)
t68 = qJ(2) + pkin(8);
t65 = cos(t68);
t84 = g(3) * t65;
t73 = sin(qJ(2));
t83 = t73 * pkin(2) + pkin(5);
t72 = sin(qJ(5));
t74 = sin(qJ(1));
t82 = t74 * t72;
t75 = cos(qJ(5));
t81 = t74 * t75;
t77 = cos(qJ(1));
t80 = t77 * t72;
t79 = t77 * t75;
t78 = g(1) * t77 + g(2) * t74;
t76 = cos(qJ(2));
t71 = pkin(6) + qJ(3);
t70 = cos(pkin(8));
t69 = sin(pkin(8));
t64 = sin(t68);
t63 = t76 * pkin(2) + pkin(1);
t62 = t71 * t74;
t61 = t77 * t71;
t60 = g(1) * t74 - g(2) * t77;
t59 = g(3) * t64 + t78 * t65;
t58 = t78 * t64 - t84;
t57 = (pkin(3) * t70 + qJ(4) * t69 + pkin(2)) * t76 + (-t69 * pkin(3) + qJ(4) * t70) * t73 + pkin(1);
t1 = [0, -t78, t60, 0, 0, 0, 0, 0, -g(3) * t73 - t78 * t76, -g(3) * t76 + t78 * t73, -t59, t58, -t60, -g(1) * (t77 * t63 + t62) - g(2) * (t74 * t63 - t61) - g(3) * t83, -t60, t59, -t58, -g(1) * (t57 * t77 + t62) - g(2) * (t57 * t74 - t61) - g(3) * (t64 * pkin(3) - t65 * qJ(4) + t83), 0, 0, 0, 0, 0, -g(1) * (t64 * t80 + t81) - g(2) * (t64 * t82 - t79) + t72 * t84, -g(1) * (t64 * t79 - t82) - g(2) * (t64 * t81 + t80) + t75 * t84;];
U_reg = t1;
