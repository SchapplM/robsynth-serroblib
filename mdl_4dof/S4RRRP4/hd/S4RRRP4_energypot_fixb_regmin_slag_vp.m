% Calculate minimal parameter regressor of potential energy for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:41
% DurationCPUTime: 0.03s
% Computational Cost: add. (34->18), mult. (39->23), div. (0->0), fcn. (36->6), ass. (0->12)
t101 = sin(qJ(1));
t103 = cos(qJ(1));
t104 = g(1) * t103 + g(2) * t101;
t102 = cos(qJ(2));
t100 = sin(qJ(2));
t99 = qJ(2) + qJ(3);
t98 = -qJ(4) - pkin(6) - pkin(5);
t97 = cos(t99);
t96 = sin(t99);
t95 = g(1) * t101 - g(2) * t103;
t94 = t102 * pkin(2) + pkin(3) * t97 + pkin(1);
t1 = [0, -t104, t95, 0, 0, 0, 0, 0, -g(3) * t100 - t104 * t102, -g(3) * t102 + t104 * t100, 0, 0, 0, 0, 0, -g(3) * t96 - t104 * t97, -g(3) * t97 + t104 * t96, -t95, -g(1) * (-t101 * t98 + t103 * t94) - g(2) * (t101 * t94 + t103 * t98) - g(3) * (t100 * pkin(2) + pkin(3) * t96 + pkin(4));];
U_reg = t1;
