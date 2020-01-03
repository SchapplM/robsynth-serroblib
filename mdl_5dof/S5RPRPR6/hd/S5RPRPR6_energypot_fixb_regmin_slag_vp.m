% Calculate minimal parameter regressor of potential energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->21), mult. (38->25), div. (0->0), fcn. (32->8), ass. (0->13)
t97 = qJ(2) + pkin(5);
t91 = qJ(1) + pkin(8);
t90 = qJ(3) + t91;
t88 = sin(t90);
t89 = cos(t90);
t86 = g(1) * t88 - g(2) * t89;
t93 = sin(qJ(1));
t95 = cos(qJ(1));
t96 = -g(1) * t95 - g(2) * t93;
t94 = cos(qJ(5));
t92 = sin(qJ(5));
t87 = g(1) * t89 + g(2) * t88;
t1 = [0, t96, g(1) * t93 - g(2) * t95, t96 * pkin(1) - g(3) * t97, 0, -t87, t86, t87, -t86, -g(1) * (t89 * pkin(3) + t88 * qJ(4) + pkin(2) * cos(t91) + t95 * pkin(1)) - g(2) * (t88 * pkin(3) - t89 * qJ(4) + pkin(2) * sin(t91) + t93 * pkin(1)) - g(3) * (pkin(6) + t97), 0, 0, 0, 0, 0, -g(3) * t94 - t86 * t92, g(3) * t92 - t86 * t94;];
U_reg = t1;
