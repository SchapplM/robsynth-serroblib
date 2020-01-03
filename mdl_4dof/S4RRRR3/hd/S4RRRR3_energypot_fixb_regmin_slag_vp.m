% Calculate minimal parameter regressor of potential energy for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:36
% EndTime: 2019-12-31 17:24:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (32->11), mult. (34->16), div. (0->0), fcn. (34->8), ass. (0->12)
t116 = qJ(2) + qJ(3);
t118 = sin(qJ(1));
t120 = cos(qJ(1));
t121 = g(1) * t120 + g(2) * t118;
t119 = cos(qJ(2));
t117 = sin(qJ(2));
t115 = qJ(4) + t116;
t114 = cos(t116);
t113 = sin(t116);
t112 = cos(t115);
t111 = sin(t115);
t1 = [0, -t121, g(1) * t118 - g(2) * t120, 0, 0, 0, 0, 0, -g(3) * t117 - t121 * t119, -g(3) * t119 + t121 * t117, 0, 0, 0, 0, 0, -g(3) * t113 - t121 * t114, -g(3) * t114 + t121 * t113, 0, 0, 0, 0, 0, -g(3) * t111 - t121 * t112, -g(3) * t112 + t121 * t111;];
U_reg = t1;
