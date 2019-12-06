% Calculate minimal parameter regressor of potential energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:49
% EndTime: 2019-12-05 16:15:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (72->24), mult. (42->28), div. (0->0), fcn. (38->10), ass. (0->14)
t103 = pkin(8) + qJ(2);
t101 = qJ(3) + t103;
t95 = sin(t101);
t96 = cos(t101);
t106 = g(1) * t96 + g(2) * t95;
t105 = cos(pkin(9));
t104 = sin(pkin(9));
t102 = pkin(9) + qJ(5);
t100 = cos(t103);
t99 = cos(t102);
t98 = sin(t103);
t97 = sin(t102);
t94 = g(1) * t95 - g(2) * t96;
t1 = [-g(3) * qJ(1), 0, -g(1) * t100 - g(2) * t98, g(1) * t98 - g(2) * t100, 0, -t106, t94, -g(3) * t104 - t106 * t105, -g(3) * t105 + t106 * t104, -t94, -g(1) * (t96 * pkin(3) + t95 * qJ(4) + pkin(2) * t100 + cos(pkin(8)) * pkin(1)) - g(2) * (t95 * pkin(3) - t96 * qJ(4) + pkin(2) * t98 + sin(pkin(8)) * pkin(1)) - g(3) * (pkin(6) + pkin(5) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t97 - t106 * t99, -g(3) * t99 + t106 * t97;];
U_reg = t1;
