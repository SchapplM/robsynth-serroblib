% Calculate inertial parameters regressor of potential energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:52
% EndTime: 2019-12-31 18:32:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (113->52), mult. (130->62), div. (0->0), fcn. (121->8), ass. (0->32)
t75 = cos(qJ(1));
t68 = pkin(8) + qJ(3);
t64 = sin(t68);
t80 = qJ(4) * t64;
t65 = cos(t68);
t86 = t65 * t75;
t91 = pkin(3) * t86 + t75 * t80;
t90 = g(3) * pkin(5);
t89 = g(3) * t65;
t69 = sin(pkin(8));
t88 = t69 * pkin(2) + pkin(5);
t73 = sin(qJ(1));
t87 = t65 * t73;
t72 = sin(qJ(5));
t85 = t73 * t72;
t74 = cos(qJ(5));
t84 = t73 * t74;
t83 = t75 * t72;
t82 = t75 * t74;
t70 = cos(pkin(8));
t62 = t70 * pkin(2) + pkin(1);
t71 = -pkin(6) - qJ(2);
t81 = t73 * t62 + t75 * t71;
t57 = t75 * t62;
t79 = -t73 * t71 + t57;
t78 = pkin(3) * t87 + t73 * t80 + t81;
t77 = g(1) * t75 + g(2) * t73;
t76 = t64 * pkin(3) - t65 * qJ(4) + t88;
t58 = g(1) * t73 - g(2) * t75;
t53 = g(3) * t64 + t77 * t65;
t52 = t77 * t64 - t89;
t1 = [0, 0, 0, 0, 0, 0, -t77, t58, -g(3), -t90, 0, 0, 0, 0, 0, 0, -g(3) * t69 - t77 * t70, -g(3) * t70 + t77 * t69, -t58, -g(1) * (t75 * pkin(1) + t73 * qJ(2)) - g(2) * (t73 * pkin(1) - t75 * qJ(2)) - t90, 0, 0, 0, 0, 0, 0, -t53, t52, -t58, -g(1) * t79 - g(2) * t81 - g(3) * t88, 0, 0, 0, 0, 0, 0, -t58, t53, -t52, -g(1) * (t79 + t91) - g(2) * t78 - g(3) * t76, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t83 + t84) - g(2) * (t64 * t85 - t82) + t72 * t89, -g(1) * (t64 * t82 - t85) - g(2) * (t64 * t84 + t83) + t74 * t89, -t53, -g(1) * (pkin(7) * t86 + t57 + (pkin(4) - t71) * t73 + t91) - g(2) * (-t75 * pkin(4) + pkin(7) * t87 + t78) - g(3) * (t64 * pkin(7) + t76);];
U_reg = t1;
