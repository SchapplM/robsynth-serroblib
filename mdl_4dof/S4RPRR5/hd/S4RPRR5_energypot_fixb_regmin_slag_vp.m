% Calculate minimal parameter regressor of potential energy for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:39
% EndTime: 2019-12-31 16:51:39
% DurationCPUTime: 0.07s
% Computational Cost: add. (22->13), mult. (45->23), div. (0->0), fcn. (50->6), ass. (0->12)
t87 = cos(qJ(3));
t86 = sin(qJ(3));
t82 = sin(qJ(1));
t84 = cos(qJ(1));
t75 = -t82 * t86 - t84 * t87;
t76 = t82 * t87 - t84 * t86;
t85 = g(1) * t75 - g(2) * t76;
t83 = cos(qJ(4));
t81 = sin(qJ(4));
t78 = -g(1) * t84 - g(2) * t82;
t77 = g(1) * t82 - g(2) * t84;
t1 = [0, t78, t77, t78, -t77, -g(1) * (pkin(1) * t84 + qJ(2) * t82) - g(2) * (pkin(1) * t82 - qJ(2) * t84) - g(3) * pkin(4), 0, t85, -g(1) * t76 - g(2) * t75, 0, 0, 0, 0, 0, g(3) * t81 + t83 * t85, g(3) * t83 - t81 * t85;];
U_reg = t1;
