% Calculate inertial parameters regressor of potential energy for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:40
% EndTime: 2019-12-31 20:53:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (109->42), mult. (104->41), div. (0->0), fcn. (91->6), ass. (0->23)
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t75 = pkin(3) * t60 + qJ(4) * t58;
t57 = qJ(1) + qJ(2);
t51 = sin(t57);
t74 = t75 * t51;
t62 = pkin(6) + pkin(5);
t72 = g(3) * t62;
t59 = sin(qJ(1));
t71 = t59 * pkin(1) + t51 * pkin(2);
t69 = qJ(5) * t60;
t52 = cos(t57);
t61 = cos(qJ(1));
t68 = t61 * pkin(1) + t52 * pkin(2) + t51 * pkin(7);
t67 = -t52 * pkin(7) + t71;
t66 = t75 * t52 + t68;
t65 = g(1) * t52 + g(2) * t51;
t64 = -g(1) * t61 - g(2) * t59;
t63 = t58 * pkin(3) - t60 * qJ(4) + t62;
t42 = g(1) * t51 - g(2) * t52;
t41 = g(3) * t58 + t65 * t60;
t40 = -g(3) * t60 + t65 * t58;
t1 = [0, 0, 0, 0, 0, 0, t64, g(1) * t59 - g(2) * t61, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t65, t42, -g(3), t64 * pkin(1) - t72, 0, 0, 0, 0, 0, 0, -t41, t40, -t42, -g(1) * t68 - g(2) * t67 - t72, 0, 0, 0, 0, 0, 0, -t42, t41, -t40, -g(1) * t66 - g(2) * (t67 + t74) - g(3) * t63, 0, 0, 0, 0, 0, 0, -t42, -t40, -t41, -g(1) * (t51 * pkin(4) + t52 * t69 + t66) - g(2) * (t51 * t69 + (-pkin(4) - pkin(7)) * t52 + t71 + t74) - g(3) * (t58 * qJ(5) + t63);];
U_reg = t1;
