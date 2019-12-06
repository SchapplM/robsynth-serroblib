% Calculate minimal parameter regressor of potential energy for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (47->13), mult. (23->17), div. (0->0), fcn. (22->8), ass. (0->13)
t90 = pkin(9) + qJ(2);
t89 = qJ(3) + t90;
t86 = qJ(4) + t89;
t82 = sin(t86);
t83 = cos(t86);
t93 = g(1) * t83 + g(2) * t82;
t92 = cos(qJ(5));
t91 = sin(qJ(5));
t88 = cos(t90);
t87 = sin(t90);
t85 = cos(t89);
t84 = sin(t89);
t1 = [-g(3) * qJ(1), 0, -g(1) * t88 - g(2) * t87, g(1) * t87 - g(2) * t88, 0, -g(1) * t85 - g(2) * t84, g(1) * t84 - g(2) * t85, 0, -t93, g(1) * t82 - g(2) * t83, 0, 0, 0, 0, 0, -g(3) * t91 - t93 * t92, -g(3) * t92 + t93 * t91;];
U_reg = t1;
