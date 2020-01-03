% Calculate minimal parameter regressor of potential energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:19
% EndTime: 2020-01-03 11:52:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->13), mult. (27->18), div. (0->0), fcn. (24->8), ass. (0->13)
t100 = qJ(1) + pkin(9) + qJ(3);
t99 = qJ(4) + t100;
t95 = sin(t99);
t96 = cos(t99);
t106 = g(2) * t95 - g(3) * t96;
t102 = sin(qJ(1));
t104 = cos(qJ(1));
t105 = -g(2) * t102 + g(3) * t104;
t103 = cos(qJ(5));
t101 = sin(qJ(5));
t98 = cos(t100);
t97 = sin(t100);
t1 = [0, t105, -g(2) * t104 - g(3) * t102, -g(1) * (qJ(2) + pkin(5)) + t105 * pkin(1), 0, -g(2) * t97 + g(3) * t98, -g(2) * t98 - g(3) * t97, 0, -t106, -g(2) * t96 - g(3) * t95, 0, 0, 0, 0, 0, -g(1) * t101 - t106 * t103, -g(1) * t103 + t106 * t101;];
U_reg = t1;
