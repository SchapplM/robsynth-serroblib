% Calculate minimal parameter regressor of potential energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x10]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:42
% EndTime: 2019-03-08 18:29:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (45->18), mult. (28->21), div. (0->0), fcn. (22->6), ass. (0->11)
t69 = qJ(2) + pkin(4);
t65 = qJ(1) + pkin(6);
t66 = sin(qJ(1));
t67 = cos(qJ(1));
t68 = -g(1) * t67 - g(2) * t66;
t64 = qJ(3) + t65;
t63 = cos(t64);
t62 = sin(t64);
t61 = -g(1) * t63 - g(2) * t62;
t60 = g(1) * t62 - g(2) * t63;
t1 = [0, t68, g(1) * t66 - g(2) * t67, t68 * pkin(1) - g(3) * t69, 0, t61, t60, t61, -t60, -g(1) * (t63 * pkin(3) + t62 * qJ(4) + pkin(2) * cos(t65) + t67 * pkin(1)) - g(2) * (t62 * pkin(3) - t63 * qJ(4) + pkin(2) * sin(t65) + t66 * pkin(1)) - g(3) * (pkin(5) + t69);];
U_reg  = t1;
