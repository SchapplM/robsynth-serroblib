% Calculate minimal parameter regressor of potential energy for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:36
% EndTime: 2020-01-03 12:07:36
% DurationCPUTime: 0.03s
% Computational Cost: add. (49->21), mult. (31->27), div. (0->0), fcn. (28->10), ass. (0->13)
t102 = qJ(1) + qJ(2);
t101 = qJ(3) + t102;
t96 = pkin(9) + t101;
t107 = g(2) * sin(t96) - g(3) * cos(t96);
t106 = cos(qJ(1));
t105 = cos(qJ(5));
t104 = sin(qJ(1));
t103 = sin(qJ(5));
t100 = cos(t102);
t99 = sin(t102);
t98 = cos(t101);
t97 = sin(t101);
t1 = [0, -g(2) * t104 + g(3) * t106, -g(2) * t106 - g(3) * t104, 0, -g(2) * t99 + g(3) * t100, -g(2) * t100 - g(3) * t99, 0, -g(2) * t97 + g(3) * t98, -g(2) * t98 - g(3) * t97, -g(1) * (qJ(4) + pkin(7) + pkin(6) + pkin(5)) - g(2) * (t104 * pkin(1) + pkin(2) * t99 + pkin(3) * t97) - g(3) * (-t106 * pkin(1) - pkin(2) * t100 - pkin(3) * t98), 0, 0, 0, 0, 0, -g(1) * t103 - t107 * t105, -g(1) * t105 + t107 * t103;];
U_reg = t1;
