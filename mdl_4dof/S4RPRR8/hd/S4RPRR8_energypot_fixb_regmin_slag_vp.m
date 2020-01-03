% Calculate minimal parameter regressor of potential energy for
% S4RPRR8
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
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.03s
% Computational Cost: add. (22->13), mult. (35->19), div. (0->0), fcn. (32->6), ass. (0->10)
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t69 = g(1) * t75 - g(2) * t77;
t76 = cos(qJ(3));
t74 = sin(qJ(3));
t73 = qJ(3) + qJ(4);
t72 = cos(t73);
t71 = sin(t73);
t70 = g(1) * t77 + g(2) * t75;
t1 = [0, -t70, t69, t70, -t69, -g(1) * (t77 * pkin(1) + t75 * qJ(2)) - g(2) * (t75 * pkin(1) - t77 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, -g(3) * t76 - t69 * t74, g(3) * t74 - t69 * t76, 0, 0, 0, 0, 0, -g(3) * t72 - t69 * t71, g(3) * t71 - t69 * t72;];
U_reg = t1;
