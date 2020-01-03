% Calculate minimal parameter regressor of potential energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:41
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->16), mult. (43->23), div. (0->0), fcn. (40->8), ass. (0->13)
t112 = pkin(7) + qJ(3);
t115 = sin(qJ(1));
t116 = cos(qJ(1));
t117 = g(1) * t116 + g(2) * t115;
t114 = cos(pkin(7));
t113 = sin(pkin(7));
t111 = qJ(4) + t112;
t110 = cos(t112);
t109 = sin(t112);
t108 = cos(t111);
t107 = sin(t111);
t106 = g(1) * t115 - g(2) * t116;
t1 = [0, -t117, t106, -g(3) * t113 - t117 * t114, -g(3) * t114 + t117 * t113, -t106, -g(1) * (t116 * pkin(1) + t115 * qJ(2)) - g(2) * (t115 * pkin(1) - t116 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, -g(3) * t109 - t117 * t110, -g(3) * t110 + t117 * t109, 0, 0, 0, 0, 0, -g(3) * t107 - t117 * t108, -g(3) * t108 + t117 * t107;];
U_reg = t1;
