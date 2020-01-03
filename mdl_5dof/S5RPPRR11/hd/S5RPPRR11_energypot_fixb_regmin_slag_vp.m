% Calculate minimal parameter regressor of potential energy for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:53
% EndTime: 2019-12-31 18:05:53
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->25), mult. (58->34), div. (0->0), fcn. (56->6), ass. (0->16)
t106 = cos(qJ(4));
t114 = g(3) * t106;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t113 = t107 * pkin(1) + t104 * qJ(2);
t102 = sin(qJ(5));
t112 = t104 * t102;
t105 = cos(qJ(5));
t111 = t104 * t105;
t110 = t107 * t102;
t109 = t107 * t105;
t108 = t104 * pkin(1) - t107 * qJ(2);
t97 = g(1) * t107 + g(2) * t104;
t103 = sin(qJ(4));
t96 = g(1) * t104 - g(2) * t107;
t1 = [0, -t97, t96, t97, -t96, -g(3) * pkin(5) - g(1) * t113 - g(2) * t108, -t96, -t97, -g(1) * (t107 * qJ(3) + t113) - g(2) * (t104 * qJ(3) + t108) - g(3) * (pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -t97 * t103 - t114, g(3) * t103 - t97 * t106, 0, 0, 0, 0, 0, -g(1) * (t103 * t109 - t112) - g(2) * (t103 * t111 + t110) - t105 * t114, -g(1) * (-t103 * t110 - t111) - g(2) * (-t103 * t112 + t109) + t102 * t114;];
U_reg = t1;
