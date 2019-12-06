% Calculate inertial parameters regressor of potential energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:31
% EndTime: 2019-12-05 16:13:31
% DurationCPUTime: 0.13s
% Computational Cost: add. (125->58), mult. (268->81), div. (0->0), fcn. (301->8), ass. (0->41)
t78 = cos(qJ(2));
t97 = pkin(2) * t78;
t96 = g(3) * qJ(1);
t72 = sin(pkin(7));
t76 = sin(qJ(2));
t95 = t72 * t76;
t74 = cos(pkin(7));
t94 = t74 * t76;
t75 = sin(qJ(3));
t93 = t75 * t76;
t92 = t75 * t78;
t77 = cos(qJ(3));
t91 = t76 * t77;
t90 = t77 * t78;
t89 = t74 * pkin(1) + t72 * pkin(5);
t88 = t72 * pkin(1) - t74 * pkin(5);
t87 = pkin(6) * t94 + t74 * t97 + t89;
t86 = t76 * pkin(2) - t78 * pkin(6) + qJ(1);
t85 = g(1) * t74 + g(2) * t72;
t84 = pkin(6) * t95 + t72 * t97 + t88;
t83 = pkin(3) * t91 + qJ(4) * t93 + t86;
t52 = t72 * t90 - t74 * t75;
t71 = sin(pkin(8));
t73 = cos(pkin(8));
t44 = t52 * t71 - t73 * t95;
t56 = t72 * t75 + t74 * t90;
t46 = t56 * t71 - t73 * t94;
t53 = t71 * t91 + t78 * t73;
t82 = g(1) * t46 + g(2) * t44 + g(3) * t53;
t55 = -t72 * t77 + t74 * t92;
t81 = t56 * pkin(3) + t55 * qJ(4) + t87;
t51 = t72 * t92 + t74 * t77;
t80 = g(1) * t55 + g(2) * t51 + g(3) * t93;
t79 = t52 * pkin(3) + t51 * qJ(4) + t84;
t57 = g(1) * t72 - g(2) * t74;
t54 = -t78 * t71 + t73 * t91;
t48 = -g(3) * t78 + t85 * t76;
t47 = t56 * t73 + t71 * t94;
t45 = t52 * t73 + t71 * t95;
t42 = -g(1) * t47 - g(2) * t45 - g(3) * t54;
t1 = [0, 0, 0, 0, 0, 0, -t85, t57, -g(3), -t96, 0, 0, 0, 0, 0, 0, -g(3) * t76 - t85 * t78, t48, -t57, -g(1) * t89 - g(2) * t88 - t96, 0, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t52 - g(3) * t91, t80, -t48, -g(1) * t87 - g(2) * t84 - g(3) * t86, 0, 0, 0, 0, 0, 0, t42, t82, -t80, -g(1) * t81 - g(2) * t79 - g(3) * t83, 0, 0, 0, 0, 0, 0, t42, -t80, -t82, -g(1) * (t47 * pkin(4) + t46 * qJ(5) + t81) - g(2) * (t45 * pkin(4) + t44 * qJ(5) + t79) - g(3) * (t54 * pkin(4) + t53 * qJ(5) + t83);];
U_reg = t1;
