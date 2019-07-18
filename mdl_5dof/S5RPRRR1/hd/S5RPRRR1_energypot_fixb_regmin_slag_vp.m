% Calculate minimal parameter regressor of potential energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:01
% EndTime: 2019-07-18 13:26:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (31->26), mult. (78->46), div. (0->0), fcn. (90->8), ass. (0->19)
t96 = sin(qJ(3));
t108 = g(3) * t96;
t97 = sin(qJ(1));
t107 = t96 * t97;
t99 = cos(qJ(4));
t106 = t96 * t99;
t100 = cos(qJ(3));
t105 = t100 * t97;
t101 = cos(qJ(1));
t104 = t101 * t96;
t103 = t100 * t101;
t92 = g(1) * t97 - g(2) * t101;
t102 = g(1) * t101 + g(2) * t97;
t98 = cos(qJ(5));
t95 = sin(qJ(4));
t94 = sin(qJ(5));
t91 = t99 * t103 + t97 * t95;
t90 = -t101 * t95 + t99 * t105;
t1 = [0, -t102, t92, -t102, -t92, -t92 * qJ(2), 0, 0, 0, 0, 0, -t102 * t100 - t108, -g(3) * t100 + t102 * t96, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t90 - g(3) * t106, -g(1) * (-t95 * t103 + t97 * t99) - g(2) * (-t101 * t99 - t95 * t105) + t95 * t108, 0, 0, 0, 0, 0, -g(1) * (t94 * t104 + t91 * t98) - g(2) * (t94 * t107 + t90 * t98) - g(3) * (-t100 * t94 + t98 * t106), -g(1) * (t98 * t104 - t91 * t94) - g(2) * (t98 * t107 - t90 * t94) - g(3) * (-t100 * t98 - t94 * t106);];
U_reg  = t1;
