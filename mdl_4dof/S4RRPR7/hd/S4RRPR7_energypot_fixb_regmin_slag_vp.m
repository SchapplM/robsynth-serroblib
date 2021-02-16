% Calculate minimal parameter regressor of potential energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:37
% EndTime: 2021-01-15 10:56:37
% DurationCPUTime: 0.15s
% Computational Cost: add. (40->24), mult. (56->35), div. (0->0), fcn. (57->8), ass. (0->19)
t58 = qJ(2) + pkin(7);
t56 = sin(t58);
t71 = g(3) * t56;
t60 = sin(qJ(4));
t62 = sin(qJ(1));
t70 = t62 * t60;
t63 = cos(qJ(4));
t69 = t62 * t63;
t65 = cos(qJ(1));
t68 = t65 * t60;
t67 = t65 * t63;
t66 = g(1) * t65 + g(2) * t62;
t64 = cos(qJ(2));
t61 = sin(qJ(2));
t59 = pkin(5) + qJ(3);
t57 = cos(t58);
t55 = t64 * pkin(2) + pkin(1);
t54 = g(1) * t62 - g(2) * t65;
t1 = [0, -t66, t54, 0, 0, 0, 0, 0, -g(3) * t61 - t66 * t64, -g(3) * t64 + t66 * t61, -t66 * t57 - t71, -g(3) * t57 + t66 * t56, -t54, -g(1) * (t65 * t55 + t59 * t62) - g(2) * (t62 * t55 - t65 * t59) - g(3) * (t61 * pkin(2) + pkin(4)), 0, 0, 0, 0, 0, -g(1) * (t57 * t67 + t70) - g(2) * (t57 * t69 - t68) - t63 * t71, -g(1) * (-t57 * t68 + t69) - g(2) * (-t57 * t70 - t67) + t60 * t71;];
U_reg = t1;
