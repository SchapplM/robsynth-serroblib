% Calculate minimal parameter regressor of potential energy for
% S4RRRR1
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
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:13
% EndTime: 2019-12-31 17:22:13
% DurationCPUTime: 0.02s
% Computational Cost: add. (30->11), mult. (22->16), div. (0->0), fcn. (22->8), ass. (0->12)
t82 = qJ(1) + qJ(2);
t81 = qJ(3) + t82;
t77 = sin(t81);
t78 = cos(t81);
t87 = g(1) * t78 + g(2) * t77;
t86 = cos(qJ(1));
t85 = cos(qJ(4));
t84 = sin(qJ(1));
t83 = sin(qJ(4));
t80 = cos(t82);
t79 = sin(t82);
t1 = [0, -g(1) * t86 - g(2) * t84, g(1) * t84 - g(2) * t86, 0, -g(1) * t80 - g(2) * t79, g(1) * t79 - g(2) * t80, 0, -t87, g(1) * t77 - g(2) * t78, 0, 0, 0, 0, 0, -g(3) * t83 - t87 * t85, -g(3) * t85 + t87 * t83;];
U_reg = t1;
