% Calculate minimal parameter regressor of potential energy for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:17
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:17:02
% EndTime: 2021-01-15 12:17:02
% DurationCPUTime: 0.05s
% Computational Cost: add. (49->30), mult. (67->42), div. (0->0), fcn. (65->8), ass. (0->19)
t61 = qJ(3) + pkin(8);
t59 = cos(t61);
t72 = g(3) * t59;
t62 = sin(qJ(5));
t64 = sin(qJ(1));
t71 = t64 * t62;
t65 = cos(qJ(5));
t70 = t64 * t65;
t67 = cos(qJ(1));
t69 = t67 * t62;
t68 = t67 * t65;
t55 = g(1) * t64 - g(2) * t67;
t66 = cos(qJ(3));
t63 = sin(qJ(3));
t60 = pkin(1) + pkin(6) + qJ(4);
t58 = sin(t61);
t57 = t63 * pkin(3) + qJ(2);
t56 = g(1) * t67 + g(2) * t64;
t1 = [0, -t56, t55, t56, -t55, -g(1) * (t67 * pkin(1) + t64 * qJ(2)) - g(2) * (t64 * pkin(1) - t67 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t66 - t55 * t63, g(3) * t63 - t55 * t66, -t55 * t58 - t72, g(3) * t58 - t55 * t59, -t56, -g(1) * (t57 * t64 + t60 * t67) - g(2) * (-t57 * t67 + t60 * t64) - g(3) * (t66 * pkin(3) + pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(1) * (t58 * t70 + t69) - g(2) * (-t58 * t68 + t71) - t65 * t72, -g(1) * (-t58 * t71 + t68) - g(2) * (t58 * t69 + t70) + t62 * t72;];
U_reg = t1;
