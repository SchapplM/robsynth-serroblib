% Calculate inertial parameters regressor of potential energy for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:45
% EndTime: 2019-12-31 20:55:45
% DurationCPUTime: 0.07s
% Computational Cost: add. (119->43), mult. (106->50), div. (0->0), fcn. (93->8), ass. (0->25)
t82 = g(3) * pkin(5);
t75 = -pkin(7) - pkin(6);
t71 = sin(qJ(2));
t81 = t71 * pkin(2) + pkin(5);
t73 = cos(qJ(2));
t62 = t73 * pkin(2) + pkin(1);
t70 = qJ(2) + qJ(3);
t65 = cos(t70);
t56 = pkin(3) * t65 + t62;
t69 = -qJ(4) + t75;
t72 = sin(qJ(1));
t74 = cos(qJ(1));
t80 = t72 * t56 + t74 * t69;
t64 = sin(t70);
t79 = pkin(3) * t64 + t81;
t78 = t74 * t56 - t72 * t69;
t77 = g(1) * t74 + g(2) * t72;
t63 = pkin(8) + t70;
t59 = sin(t63);
t60 = cos(t63);
t76 = pkin(4) * t60 + qJ(5) * t59;
t57 = g(1) * t72 - g(2) * t74;
t53 = -g(3) * t59 - t77 * t60;
t52 = -g(3) * t60 + t77 * t59;
t1 = [0, 0, 0, 0, 0, 0, -t77, t57, -g(3), -t82, 0, 0, 0, 0, 0, 0, -g(3) * t71 - t77 * t73, -g(3) * t73 + t77 * t71, -t57, -g(1) * (t74 * pkin(1) + t72 * pkin(6)) - g(2) * (t72 * pkin(1) - t74 * pkin(6)) - t82, 0, 0, 0, 0, 0, 0, -g(3) * t64 - t77 * t65, -g(3) * t65 + t77 * t64, -t57, -g(1) * (t74 * t62 - t72 * t75) - g(2) * (t72 * t62 + t74 * t75) - g(3) * t81, 0, 0, 0, 0, 0, 0, t53, t52, -t57, -g(1) * t78 - g(2) * t80 - g(3) * t79, 0, 0, 0, 0, 0, 0, t53, -t57, -t52, -g(1) * (t76 * t74 + t78) - g(2) * (t76 * t72 + t80) - g(3) * (t59 * pkin(4) - t60 * qJ(5) + t79);];
U_reg = t1;
