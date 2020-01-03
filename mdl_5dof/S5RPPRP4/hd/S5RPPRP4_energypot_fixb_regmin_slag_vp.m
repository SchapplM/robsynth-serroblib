% Calculate minimal parameter regressor of potential energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:19
% EndTime: 2019-12-31 17:52:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->28), mult. (84->37), div. (0->0), fcn. (89->6), ass. (0->20)
t121 = sin(qJ(1));
t120 = -qJ(3) + pkin(5);
t119 = cos(pkin(7));
t118 = sin(pkin(7));
t111 = cos(qJ(1));
t117 = t111 * pkin(1) + t121 * qJ(2);
t116 = t111 * pkin(2) + t117;
t115 = t121 * pkin(1) - t111 * qJ(2);
t95 = -t111 * t119 - t121 * t118;
t96 = t111 * t118 - t121 * t119;
t114 = g(1) * t96 - g(2) * t95;
t113 = g(1) * t95 + g(2) * t96;
t112 = t121 * pkin(2) + t115;
t110 = cos(qJ(4));
t109 = sin(qJ(4));
t108 = -qJ(5) - pkin(6);
t102 = t110 * pkin(4) + pkin(3);
t98 = -g(1) * t111 - g(2) * t121;
t97 = g(1) * t121 - g(2) * t111;
t1 = [0, t98, t97, t98, -t97, -g(3) * pkin(5) - g(1) * t117 - g(2) * t115, t113, t114, -g(1) * t116 - g(2) * t112 - g(3) * t120, 0, 0, 0, 0, 0, g(3) * t109 + t113 * t110, g(3) * t110 - t113 * t109, -t114, -g(1) * (-t95 * t102 - t96 * t108 + t116) - g(2) * (-t96 * t102 + t95 * t108 + t112) - g(3) * (-t109 * pkin(4) + t120);];
U_reg = t1;
