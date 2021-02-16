% Calculate minimal parameter regressor of potential energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:43
% EndTime: 2021-01-15 23:52:43
% DurationCPUTime: 0.04s
% Computational Cost: add. (77->23), mult. (62->29), div. (0->0), fcn. (59->8), ass. (0->17)
t69 = qJ(2) + qJ(3);
t71 = sin(qJ(1));
t73 = cos(qJ(1));
t74 = g(1) * t73 + g(2) * t71;
t72 = cos(qJ(2));
t70 = sin(qJ(2));
t68 = qJ(4) + t69;
t67 = -qJ(5) - pkin(8) - pkin(7) - pkin(6);
t66 = cos(t69);
t65 = sin(t69);
t64 = cos(t68);
t63 = sin(t68);
t62 = g(1) * t71 - g(2) * t73;
t61 = t72 * pkin(2) + pkin(3) * t66 + pkin(4) * t64 + pkin(1);
t60 = -g(3) * t63 - t74 * t64;
t59 = -g(3) * t64 + t74 * t63;
t1 = [0, -t74, t62, 0, 0, 0, 0, 0, -g(3) * t70 - t74 * t72, -g(3) * t72 + t74 * t70, 0, 0, 0, 0, 0, -g(3) * t65 - t74 * t66, -g(3) * t66 + t74 * t65, 0, 0, 0, 0, 0, t60, t59, t60, t59, -t62, -g(1) * (t73 * t61 - t71 * t67) - g(2) * (t71 * t61 + t73 * t67) - g(3) * (t70 * pkin(2) + pkin(3) * t65 + pkin(4) * t63 + pkin(5));];
U_reg = t1;
