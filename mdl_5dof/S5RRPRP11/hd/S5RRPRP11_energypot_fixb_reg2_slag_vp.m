% Calculate inertial parameters regressor of potential energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:03
% EndTime: 2019-12-31 20:14:03
% DurationCPUTime: 0.09s
% Computational Cost: add. (88->53), mult. (171->59), div. (0->0), fcn. (170->6), ass. (0->36)
t75 = sin(qJ(1));
t74 = sin(qJ(2));
t86 = qJ(3) * t74;
t77 = cos(qJ(2));
t91 = t75 * t77;
t97 = pkin(2) * t91 + t75 * t86;
t96 = g(3) * pkin(5);
t95 = g(3) * t77;
t94 = t74 * pkin(2) + pkin(5);
t73 = sin(qJ(4));
t93 = t75 * t73;
t76 = cos(qJ(4));
t92 = t75 * t76;
t78 = cos(qJ(1));
t90 = t77 * t78;
t89 = t78 * t73;
t88 = t78 * t76;
t87 = t78 * pkin(1) + t75 * pkin(6);
t70 = t75 * pkin(1);
t85 = -t78 * pkin(6) + t70;
t84 = pkin(2) * t90 + t78 * t86 + t87;
t83 = -t77 * qJ(3) + t94;
t82 = g(1) * t78 + g(2) * t75;
t81 = t75 * pkin(3) + pkin(7) * t90 + t84;
t80 = pkin(7) * t91 + t70 + (-pkin(3) - pkin(6)) * t78 + t97;
t54 = -t74 * t88 + t93;
t56 = t74 * t92 + t89;
t79 = g(1) * t54 - g(2) * t56 + t76 * t95;
t66 = t74 * pkin(7);
t58 = g(1) * t75 - g(2) * t78;
t57 = t74 * t93 - t88;
t55 = t74 * t89 + t92;
t53 = g(3) * t74 + t82 * t77;
t52 = t82 * t74 - t95;
t51 = -g(1) * t55 - g(2) * t57 + t73 * t95;
t1 = [0, 0, 0, 0, 0, 0, -t82, t58, -g(3), -t96, 0, 0, 0, 0, 0, 0, -t53, t52, -t58, -g(1) * t87 - g(2) * t85 - t96, 0, 0, 0, 0, 0, 0, -t58, t53, -t52, -g(1) * t84 - g(2) * (t85 + t97) - g(3) * t83, 0, 0, 0, 0, 0, 0, t51, t79, -t53, -g(1) * t81 - g(2) * t80 - g(3) * (t66 + t83), 0, 0, 0, 0, 0, 0, t51, -t53, -t79, -g(1) * (t55 * pkin(4) + t54 * qJ(5) + t81) - g(2) * (t57 * pkin(4) - t56 * qJ(5) + t80) - g(3) * (t66 + t94) - (-pkin(4) * t73 + qJ(5) * t76 - qJ(3)) * t95;];
U_reg = t1;
