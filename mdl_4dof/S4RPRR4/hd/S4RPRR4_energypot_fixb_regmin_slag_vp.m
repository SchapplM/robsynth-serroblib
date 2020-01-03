% Calculate minimal parameter regressor of potential energy for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->16), mult. (39->28), div. (0->0), fcn. (40->8), ass. (0->15)
t100 = sin(qJ(3));
t109 = g(3) * t100;
t103 = cos(qJ(3));
t99 = sin(qJ(4));
t108 = t103 * t99;
t102 = cos(qJ(4));
t107 = t102 * t103;
t98 = qJ(1) + pkin(7);
t96 = sin(t98);
t97 = cos(t98);
t106 = g(1) * t97 + g(2) * t96;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t105 = -g(1) * t104 - g(2) * t101;
t1 = [0, t105, g(1) * t101 - g(2) * t104, -g(3) * (qJ(2) + pkin(4)) + t105 * pkin(1), 0, 0, 0, 0, 0, -t106 * t103 - t109, -g(3) * t103 + t106 * t100, 0, 0, 0, 0, 0, -g(1) * (t97 * t107 + t96 * t99) - g(2) * (t96 * t107 - t97 * t99) - t102 * t109, -g(1) * (t96 * t102 - t97 * t108) - g(2) * (-t97 * t102 - t96 * t108) + t99 * t109;];
U_reg = t1;
