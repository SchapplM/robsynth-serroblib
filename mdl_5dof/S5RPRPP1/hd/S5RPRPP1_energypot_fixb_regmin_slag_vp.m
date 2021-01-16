% Calculate minimal parameter regressor of potential energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:03
% EndTime: 2021-01-15 11:14:03
% DurationCPUTime: 0.05s
% Computational Cost: add. (95->29), mult. (77->38), div. (0->0), fcn. (68->8), ass. (0->23)
t77 = qJ(2) + pkin(5);
t69 = cos(qJ(3));
t56 = t69 * pkin(3) + pkin(2);
t65 = qJ(1) + pkin(7);
t58 = sin(t65);
t60 = cos(t65);
t66 = -qJ(4) - pkin(6);
t68 = sin(qJ(1));
t76 = t68 * pkin(1) + t58 * t56 + t60 * t66;
t67 = sin(qJ(3));
t75 = t67 * pkin(3) + t77;
t70 = cos(qJ(1));
t74 = t70 * pkin(1) + t60 * t56 - t58 * t66;
t73 = g(1) * t60 + g(2) * t58;
t72 = -g(1) * t70 - g(2) * t68;
t64 = qJ(3) + pkin(8);
t57 = sin(t64);
t59 = cos(t64);
t71 = pkin(4) * t59 + qJ(5) * t57;
t52 = -g(1) * t58 + g(2) * t60;
t51 = -g(3) * t57 - t73 * t59;
t50 = -g(3) * t59 + t73 * t57;
t1 = [0, t72, g(1) * t68 - g(2) * t70, t72 * pkin(1) - g(3) * t77, 0, 0, 0, 0, 0, -g(3) * t67 - t73 * t69, -g(3) * t69 + t73 * t67, t51, t50, t52, -g(1) * t74 - g(2) * t76 - g(3) * t75, t51, t52, -t50, -g(1) * (t71 * t60 + t74) - g(2) * (t71 * t58 + t76) - g(3) * (t57 * pkin(4) - t59 * qJ(5) + t75);];
U_reg = t1;
