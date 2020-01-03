% Calculate minimal parameter regressor of potential energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->24), mult. (62->34), div. (0->0), fcn. (64->8), ass. (0->18)
t114 = pkin(6) + pkin(5);
t113 = cos(pkin(8));
t112 = sin(pkin(8));
t104 = qJ(1) + qJ(2);
t100 = sin(t104);
t101 = cos(t104);
t108 = cos(qJ(1));
t111 = t108 * pkin(1) + t101 * pkin(2) + t100 * qJ(3);
t91 = -t100 * t112 - t101 * t113;
t92 = t100 * t113 - t101 * t112;
t110 = g(1) * t91 - g(2) * t92;
t106 = sin(qJ(1));
t109 = t106 * pkin(1) + t100 * pkin(2) - t101 * qJ(3);
t107 = cos(qJ(5));
t105 = sin(qJ(5));
t94 = -g(1) * t101 - g(2) * t100;
t93 = g(1) * t100 - g(2) * t101;
t1 = [0, -g(1) * t108 - g(2) * t106, g(1) * t106 - g(2) * t108, 0, t94, t93, t94, -t93, -g(1) * t111 - g(2) * t109 - g(3) * t114, t110, -g(1) * t92 - g(2) * t91, -g(1) * (t101 * pkin(3) + t111) - g(2) * (t100 * pkin(3) + t109) - g(3) * (-qJ(4) + t114), 0, 0, 0, 0, 0, g(3) * t105 + t110 * t107, g(3) * t107 - t110 * t105;];
U_reg = t1;
