% Calculate inertial parameters regressor of potential energy for
% S5RRRPP2
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:51:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (109->40), mult. (104->43), div. (0->0), fcn. (91->6), ass. (0->24)
t59 = qJ(1) + qJ(2);
t53 = sin(t59);
t60 = sin(qJ(3));
t71 = qJ(4) * t60;
t62 = cos(qJ(3));
t74 = t53 * t62;
t76 = pkin(3) * t74 + t53 * t71;
t64 = pkin(6) + pkin(5);
t75 = g(3) * t64;
t54 = cos(t59);
t73 = t54 * t62;
t61 = sin(qJ(1));
t72 = t61 * pkin(1) + t53 * pkin(2);
t63 = cos(qJ(1));
t70 = t63 * pkin(1) + t54 * pkin(2) + t53 * pkin(7);
t69 = -t54 * pkin(7) + t72;
t68 = pkin(3) * t73 + t54 * t71 + t70;
t67 = g(1) * t54 + g(2) * t53;
t66 = -g(1) * t63 - g(2) * t61;
t65 = t60 * pkin(3) - t62 * qJ(4) + t64;
t44 = g(1) * t53 - g(2) * t54;
t43 = -g(3) * t60 - t67 * t62;
t42 = -g(3) * t62 + t67 * t60;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t61 - g(2) * t63, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t67, t44, -g(3), t66 * pkin(1) - t75, 0, 0, 0, 0, 0, 0, t43, t42, -t44, -g(1) * t70 - g(2) * t69 - t75, 0, 0, 0, 0, 0, 0, t43, -t44, -t42, -g(1) * t68 - g(2) * (t69 + t76) - g(3) * t65, 0, 0, 0, 0, 0, 0, t43, -t42, t44, -g(1) * (pkin(4) * t73 - t53 * qJ(5) + t68) - g(2) * (pkin(4) * t74 + (-pkin(7) + qJ(5)) * t54 + t72 + t76) - g(3) * (t60 * pkin(4) + t65);];
U_reg = t1;
