% Calculate minimal parameter regressor of potential energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:13
% EndTime: 2021-01-15 11:44:13
% DurationCPUTime: 0.05s
% Computational Cost: add. (66->24), mult. (53->33), div. (0->0), fcn. (47->10), ass. (0->19)
t75 = qJ(2) + pkin(5);
t66 = qJ(3) + pkin(9);
t67 = qJ(1) + pkin(8);
t62 = sin(t67);
t64 = cos(t67);
t74 = g(2) * t62 - g(3) * t64;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t73 = -g(2) * t70 + g(3) * t72;
t71 = cos(qJ(3));
t69 = sin(qJ(3));
t68 = -qJ(4) - pkin(6);
t65 = qJ(5) + t66;
t63 = cos(t66);
t61 = sin(t66);
t60 = t71 * pkin(3) + pkin(2);
t59 = cos(t65);
t58 = sin(t65);
t1 = [0, t73, -g(2) * t72 - g(3) * t70, t73 * pkin(1) - g(1) * t75, 0, 0, 0, 0, 0, -g(1) * t69 - t74 * t71, -g(1) * t71 + t74 * t69, -g(1) * t61 - t74 * t63, -g(1) * t63 + t74 * t61, g(2) * t64 + g(3) * t62, -g(1) * (t69 * pkin(3) + t75) - g(2) * (t70 * pkin(1) + t62 * t60 + t64 * t68) - g(3) * (-t72 * pkin(1) - t64 * t60 + t62 * t68), 0, 0, 0, 0, 0, -g(1) * t58 - t74 * t59, -g(1) * t59 + t74 * t58;];
U_reg = t1;
