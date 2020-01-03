% Calculate minimal parameter regressor of potential energy for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:53
% EndTime: 2019-12-31 17:04:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (32->16), mult. (36->21), div. (0->0), fcn. (33->6), ass. (0->12)
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t122 = g(1) * t121 + g(2) * t119;
t120 = cos(qJ(2));
t118 = sin(qJ(2));
t117 = -pkin(5) - qJ(3);
t116 = qJ(2) + pkin(7) + qJ(4);
t115 = t120 * pkin(2) + pkin(1);
t114 = cos(t116);
t113 = sin(t116);
t112 = g(1) * t119 - g(2) * t121;
t1 = [0, -t122, t112, 0, 0, 0, 0, 0, -g(3) * t118 - t122 * t120, -g(3) * t120 + t122 * t118, -t112, -g(1) * (t121 * t115 - t119 * t117) - g(2) * (t119 * t115 + t121 * t117) - g(3) * (t118 * pkin(2) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t113 - t122 * t114, -g(3) * t114 + t122 * t113;];
U_reg = t1;
