% Calculate minimal parameter regressor of potential energy for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (48->21), mult. (56->25), div. (0->0), fcn. (53->6), ass. (0->14)
t103 = sin(qJ(1));
t105 = cos(qJ(1));
t108 = g(1) * t105 + g(2) * t103;
t101 = qJ(2) + qJ(3);
t100 = cos(t101);
t104 = cos(qJ(2));
t99 = sin(t101);
t107 = t104 * pkin(2) + pkin(3) * t100 + qJ(4) * t99 + pkin(1);
t106 = -pkin(6) - pkin(5);
t102 = sin(qJ(2));
t97 = g(1) * t103 - g(2) * t105;
t96 = -g(3) * t99 - t108 * t100;
t95 = -g(3) * t100 + t108 * t99;
t1 = [0, -t108, t97, 0, 0, 0, 0, 0, -g(3) * t102 - t108 * t104, -g(3) * t104 + t108 * t102, 0, 0, 0, 0, 0, t96, t95, t96, -t97, -t95, -g(3) * (t102 * pkin(2) + t99 * pkin(3) - t100 * qJ(4) + pkin(4)) + (-g(1) * t107 - g(2) * t106) * t105 + (g(1) * t106 - g(2) * t107) * t103;];
U_reg = t1;
