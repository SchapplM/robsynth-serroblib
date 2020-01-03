% Calculate inertial parameters regressor of potential energy for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:25
% EndTime: 2019-12-31 21:20:25
% DurationCPUTime: 0.09s
% Computational Cost: add. (113->52), mult. (130->62), div. (0->0), fcn. (121->8), ass. (0->32)
t78 = cos(qJ(1));
t72 = qJ(2) + qJ(3);
t68 = sin(t72);
t84 = qJ(4) * t68;
t69 = cos(t72);
t90 = t69 * t78;
t95 = pkin(3) * t90 + t78 * t84;
t94 = g(3) * pkin(5);
t93 = g(3) * t69;
t74 = sin(qJ(2));
t92 = t74 * pkin(2) + pkin(5);
t75 = sin(qJ(1));
t91 = t69 * t75;
t73 = sin(qJ(5));
t89 = t75 * t73;
t76 = cos(qJ(5));
t88 = t75 * t76;
t87 = t78 * t73;
t86 = t78 * t76;
t77 = cos(qJ(2));
t66 = t77 * pkin(2) + pkin(1);
t79 = -pkin(7) - pkin(6);
t85 = t75 * t66 + t78 * t79;
t62 = t78 * t66;
t83 = -t75 * t79 + t62;
t82 = pkin(3) * t91 + t75 * t84 + t85;
t81 = g(1) * t78 + g(2) * t75;
t80 = t68 * pkin(3) - t69 * qJ(4) + t92;
t58 = g(1) * t75 - g(2) * t78;
t57 = g(3) * t68 + t81 * t69;
t56 = t81 * t68 - t93;
t1 = [0, 0, 0, 0, 0, 0, -t81, t58, -g(3), -t94, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t81 * t77, -g(3) * t77 + t81 * t74, -t58, -g(1) * (t78 * pkin(1) + t75 * pkin(6)) - g(2) * (t75 * pkin(1) - t78 * pkin(6)) - t94, 0, 0, 0, 0, 0, 0, -t57, t56, -t58, -g(1) * t83 - g(2) * t85 - g(3) * t92, 0, 0, 0, 0, 0, 0, -t58, t57, -t56, -g(1) * (t83 + t95) - g(2) * t82 - g(3) * t80, 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t87 + t88) - g(2) * (t68 * t89 - t86) + t73 * t93, -g(1) * (t68 * t86 - t89) - g(2) * (t68 * t88 + t87) + t76 * t93, -t57, -g(1) * (pkin(8) * t90 + t62 + (pkin(4) - t79) * t75 + t95) - g(2) * (-t78 * pkin(4) + pkin(8) * t91 + t82) - g(3) * (t68 * pkin(8) + t80);];
U_reg = t1;
