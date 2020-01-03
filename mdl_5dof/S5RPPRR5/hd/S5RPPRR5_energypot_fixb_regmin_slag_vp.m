% Calculate minimal parameter regressor of potential energy for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (52->19), mult. (52->30), div. (0->0), fcn. (54->8), ass. (0->15)
t107 = g(3) * (qJ(2) + pkin(5));
t106 = cos(qJ(4));
t105 = sin(qJ(4));
t97 = qJ(1) + pkin(8);
t95 = sin(t97);
t96 = cos(t97);
t91 = -t95 * t105 - t96 * t106;
t92 = -t96 * t105 + t95 * t106;
t104 = g(1) * t91 - g(2) * t92;
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t103 = -g(1) * t102 - g(2) * t100;
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t1 = [0, t103, g(1) * t100 - g(2) * t102, t103 * pkin(1) - t107, -g(1) * t96 - g(2) * t95, -g(1) * t95 + g(2) * t96, -g(1) * (t102 * pkin(1) + t96 * pkin(2) + t95 * qJ(3)) - g(2) * (t100 * pkin(1) + t95 * pkin(2) - t96 * qJ(3)) - t107, 0, t104, -g(1) * t92 - g(2) * t91, 0, 0, 0, 0, 0, g(3) * t99 + t104 * t101, g(3) * t101 - t104 * t99;];
U_reg = t1;
