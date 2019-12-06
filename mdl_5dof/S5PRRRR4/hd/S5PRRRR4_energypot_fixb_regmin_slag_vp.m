% Calculate minimal parameter regressor of potential energy for
% S5PRRRR4
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:58
% EndTime: 2019-12-05 17:07:58
% DurationCPUTime: 0.03s
% Computational Cost: add. (47->13), mult. (29->17), div. (0->0), fcn. (28->8), ass. (0->13)
t106 = pkin(9) + qJ(2);
t103 = qJ(3) + t106;
t100 = cos(t103);
t99 = sin(t103);
t110 = g(1) * t100 + g(2) * t99;
t109 = cos(qJ(4));
t108 = sin(qJ(4));
t107 = qJ(4) + qJ(5);
t105 = cos(t107);
t104 = sin(t107);
t102 = cos(t106);
t101 = sin(t106);
t1 = [-g(3) * qJ(1), 0, -g(1) * t102 - g(2) * t101, g(1) * t101 - g(2) * t102, 0, -t110, g(1) * t99 - g(2) * t100, 0, 0, 0, 0, 0, -g(3) * t108 - t110 * t109, -g(3) * t109 + t110 * t108, 0, 0, 0, 0, 0, -g(3) * t104 - t110 * t105, -g(3) * t105 + t110 * t104;];
U_reg = t1;
