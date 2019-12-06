% Calculate minimal parameter regressor of potential energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:50
% EndTime: 2019-12-05 15:12:50
% DurationCPUTime: 0.05s
% Computational Cost: add. (51->20), mult. (48->31), div. (0->0), fcn. (48->8), ass. (0->18)
t95 = pkin(9) + qJ(3);
t94 = qJ(4) + t95;
t90 = sin(t94);
t106 = g(3) * t90;
t105 = g(3) * qJ(1);
t96 = sin(pkin(8));
t98 = sin(qJ(5));
t104 = t96 * t98;
t99 = cos(qJ(5));
t103 = t96 * t99;
t97 = cos(pkin(8));
t102 = t97 * t98;
t101 = t97 * t99;
t100 = g(1) * t97 + g(2) * t96;
t93 = cos(t95);
t92 = sin(t95);
t91 = cos(t94);
t1 = [-t105, -g(1) * (t97 * pkin(1) + t96 * qJ(2)) - g(2) * (t96 * pkin(1) - t97 * qJ(2)) - t105, 0, -g(3) * t92 - t100 * t93, -g(3) * t93 + t100 * t92, 0, -t100 * t91 - t106, -g(3) * t91 + t100 * t90, 0, 0, 0, 0, 0, -g(1) * (t91 * t101 + t104) - g(2) * (t91 * t103 - t102) - t99 * t106, -g(1) * (-t91 * t102 + t103) - g(2) * (-t91 * t104 - t101) + t98 * t106;];
U_reg = t1;
