% Calculate minimal parameter regressor of potential energy for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:28
% EndTime: 2020-01-03 11:45:28
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->24), mult. (39->27), div. (0->0), fcn. (33->8), ass. (0->15)
t115 = qJ(2) + pkin(5);
t107 = qJ(1) + pkin(8);
t106 = qJ(3) + t107;
t103 = sin(t106);
t104 = cos(t106);
t114 = g(2) * t103 - g(3) * t104;
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t113 = -g(2) * t110 + g(3) * t112;
t111 = cos(qJ(4));
t109 = sin(qJ(4));
t108 = -qJ(5) - pkin(7);
t105 = t111 * pkin(4) + pkin(3);
t102 = g(2) * t104 + g(3) * t103;
t1 = [0, t113, -g(2) * t112 - g(3) * t110, t113 * pkin(1) - g(1) * t115, 0, -t114, -t102, 0, 0, 0, 0, 0, -g(1) * t109 - t114 * t111, -g(1) * t111 + t114 * t109, t102, -g(1) * (t109 * pkin(4) + pkin(6) + t115) - g(2) * (t103 * t105 + t104 * t108 + pkin(2) * sin(t107) + t110 * pkin(1)) - g(3) * (-t104 * t105 + t103 * t108 - pkin(2) * cos(t107) - t112 * pkin(1));];
U_reg = t1;
