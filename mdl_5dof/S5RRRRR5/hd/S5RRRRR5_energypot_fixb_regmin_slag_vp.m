% Calculate minimal parameter regressor of potential energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:05
% EndTime: 2022-01-20 12:02:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (48->14), mult. (32->20), div. (0->0), fcn. (32->10), ass. (0->15)
t124 = qJ(1) + qJ(2);
t122 = qJ(3) + t124;
t116 = sin(t122);
t117 = cos(t122);
t129 = g(1) * t117 + g(2) * t116;
t128 = cos(qJ(1));
t127 = cos(qJ(4));
t126 = sin(qJ(1));
t125 = sin(qJ(4));
t123 = qJ(4) + qJ(5);
t121 = cos(t124);
t120 = cos(t123);
t119 = sin(t124);
t118 = sin(t123);
t1 = [0, -g(1) * t128 - g(2) * t126, g(1) * t126 - g(2) * t128, 0, -g(1) * t121 - g(2) * t119, g(1) * t119 - g(2) * t121, 0, -t129, g(1) * t116 - g(2) * t117, 0, 0, 0, 0, 0, -g(3) * t125 - t129 * t127, -g(3) * t127 + t129 * t125, 0, 0, 0, 0, 0, -g(3) * t118 - t129 * t120, -g(3) * t120 + t129 * t118;];
U_reg = t1;
