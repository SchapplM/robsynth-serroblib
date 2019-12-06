% Calculate minimal parameter regressor of potential energy for
% S5PRRRP3
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
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:12
% EndTime: 2019-12-05 16:44:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (56->23), mult. (42->26), div. (0->0), fcn. (38->8), ass. (0->13)
t109 = pkin(8) + qJ(2);
t104 = sin(t109);
t105 = cos(t109);
t113 = g(1) * t105 + g(2) * t104;
t112 = cos(qJ(3));
t111 = sin(qJ(3));
t110 = qJ(3) + qJ(4);
t108 = -qJ(5) - pkin(7) - pkin(6);
t107 = cos(t110);
t106 = sin(t110);
t103 = t112 * pkin(3) + pkin(4) * t107 + pkin(2);
t102 = g(1) * t104 - g(2) * t105;
t1 = [-g(3) * qJ(1), 0, -t113, t102, 0, 0, 0, 0, 0, -g(3) * t111 - t113 * t112, -g(3) * t112 + t113 * t111, 0, 0, 0, 0, 0, -g(3) * t106 - t113 * t107, -g(3) * t107 + t113 * t106, -t102, -g(1) * (t105 * t103 - t104 * t108 + cos(pkin(8)) * pkin(1)) - g(2) * (t104 * t103 + t105 * t108 + sin(pkin(8)) * pkin(1)) - g(3) * (t111 * pkin(3) + pkin(4) * t106 + pkin(5) + qJ(1));];
U_reg = t1;
