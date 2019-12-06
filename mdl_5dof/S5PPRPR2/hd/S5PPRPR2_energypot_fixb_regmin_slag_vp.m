% Calculate minimal parameter regressor of potential energy for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 0.05s
% Computational Cost: add. (60->30), mult. (68->40), div. (0->0), fcn. (65->8), ass. (0->19)
t92 = pkin(8) + qJ(3);
t91 = cos(t92);
t105 = g(3) * t91;
t104 = g(3) * qJ(1);
t93 = sin(pkin(7));
t96 = sin(qJ(5));
t103 = t93 * t96;
t97 = cos(qJ(5));
t102 = t93 * t97;
t94 = cos(pkin(7));
t101 = t94 * t96;
t100 = t94 * t97;
t99 = g(1) * t94 + g(2) * t93;
t90 = sin(t92);
t98 = pkin(3) * t91 + qJ(4) * t90 + cos(pkin(8)) * pkin(2) + pkin(1);
t95 = -pkin(5) - qJ(2);
t88 = g(3) * t90 + t99 * t91;
t87 = t99 * t90 - t105;
t1 = [-t104, -g(1) * (t94 * pkin(1) + t93 * qJ(2)) - g(2) * (t93 * pkin(1) - t94 * qJ(2)) - t104, 0, -t88, t87, t88, -t87, -g(3) * (t90 * pkin(3) - t91 * qJ(4) + sin(pkin(8)) * pkin(2) + qJ(1)) + (-g(1) * t98 - g(2) * t95) * t94 + (g(1) * t95 - g(2) * t98) * t93, 0, 0, 0, 0, 0, -g(1) * (t90 * t101 + t102) - g(2) * (t90 * t103 - t100) + t96 * t105, -g(1) * (t90 * t100 - t103) - g(2) * (t90 * t102 + t101) + t97 * t105;];
U_reg = t1;
