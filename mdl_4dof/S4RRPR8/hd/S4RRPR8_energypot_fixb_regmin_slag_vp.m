% Calculate minimal parameter regressor of potential energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:22
% EndTime: 2019-12-31 17:08:22
% DurationCPUTime: 0.05s
% Computational Cost: add. (31->19), mult. (69->27), div. (0->0), fcn. (72->6), ass. (0->14)
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t116 = g(1) * t112 + g(2) * t109;
t107 = sin(qJ(4));
t108 = sin(qJ(2));
t110 = cos(qJ(4));
t111 = cos(qJ(2));
t115 = t111 * t107 - t108 * t110;
t114 = t108 * t107 + t111 * t110;
t113 = pkin(2) * t111 + qJ(3) * t108 + pkin(1);
t106 = g(1) * t109 - g(2) * t112;
t105 = -g(3) * t108 - t116 * t111;
t104 = -g(3) * t111 + t116 * t108;
t1 = [0, -t116, t106, 0, 0, 0, 0, 0, t105, t104, t105, -t106, -t104, -g(3) * (t108 * pkin(2) - t111 * qJ(3) + pkin(4)) + (g(2) * pkin(5) - g(1) * t113) * t112 + (-g(1) * pkin(5) - g(2) * t113) * t109, 0, 0, 0, 0, 0, g(3) * t115 - t116 * t114, g(3) * t114 + t116 * t115;];
U_reg = t1;
