% Calculate minimal parameter regressor of potential energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:18
% EndTime: 2019-07-18 13:30:18
% DurationCPUTime: 0.02s
% Computational Cost: add. (31->12), mult. (23->17), div. (0->0), fcn. (22->8), ass. (0->12)
t82 = qJ(2) + qJ(3);
t81 = qJ(4) + t82;
t77 = sin(t81);
t78 = cos(t81);
t87 = g(1) * t78 + g(2) * t77;
t86 = cos(qJ(2));
t85 = cos(qJ(5));
t84 = sin(qJ(2));
t83 = sin(qJ(5));
t80 = cos(t82);
t79 = sin(t82);
t1 = [-g(3) * qJ(1), 0, -g(1) * t86 - g(2) * t84, g(1) * t84 - g(2) * t86, 0, -g(1) * t80 - g(2) * t79, g(1) * t79 - g(2) * t80, 0, -t87, g(1) * t77 - g(2) * t78, 0, 0, 0, 0, 0, -g(3) * t83 - t87 * t85, -g(3) * t85 + t87 * t83;];
U_reg  = t1;
