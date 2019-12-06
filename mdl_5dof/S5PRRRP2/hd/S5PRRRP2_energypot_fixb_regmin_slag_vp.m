% Calculate minimal parameter regressor of potential energy for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:56
% EndTime: 2019-12-05 16:41:57
% DurationCPUTime: 0.05s
% Computational Cost: add. (80->26), mult. (52->30), div. (0->0), fcn. (48->8), ass. (0->14)
t95 = pkin(8) + qJ(2);
t94 = qJ(3) + t95;
t90 = sin(t94);
t91 = cos(t94);
t99 = g(1) * t91 + g(2) * t90;
t96 = sin(qJ(4));
t97 = cos(qJ(4));
t98 = pkin(4) * t97 + qJ(5) * t96 + pkin(3);
t93 = cos(t95);
t92 = sin(t95);
t89 = g(1) * t90 - g(2) * t91;
t88 = -g(3) * t96 - t99 * t97;
t87 = -g(3) * t97 + t99 * t96;
t1 = [-g(3) * qJ(1), 0, -g(1) * t93 - g(2) * t92, g(1) * t92 - g(2) * t93, 0, -t99, t89, 0, 0, 0, 0, 0, t88, t87, t88, -t89, -t87, -g(1) * (cos(pkin(8)) * pkin(1) + pkin(2) * t93) - g(2) * (sin(pkin(8)) * pkin(1) + pkin(2) * t92) - g(3) * (t96 * pkin(4) - t97 * qJ(5) + pkin(5) + pkin(6) + qJ(1)) + (g(2) * pkin(7) - g(1) * t98) * t91 + (-g(1) * pkin(7) - g(2) * t98) * t90;];
U_reg = t1;
