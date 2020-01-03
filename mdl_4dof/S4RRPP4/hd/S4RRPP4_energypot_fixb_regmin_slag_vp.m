% Calculate minimal parameter regressor of potential energy for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:17
% EndTime: 2019-12-31 16:59:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (41->25), mult. (79->30), div. (0->0), fcn. (73->4), ass. (0->15)
t93 = sin(qJ(2));
t104 = qJ(3) * t93 + pkin(1);
t94 = sin(qJ(1));
t95 = cos(qJ(2));
t103 = t94 * t95;
t96 = cos(qJ(1));
t102 = t95 * t96;
t100 = pkin(2) * t103 + t104 * t94;
t99 = pkin(2) * t102 + t94 * pkin(5) + t104 * t96;
t98 = t93 * pkin(2) - t95 * qJ(3) + pkin(4);
t97 = g(1) * t96 + g(2) * t94;
t83 = g(1) * t94 - g(2) * t96;
t82 = -g(3) * t93 - t97 * t95;
t81 = -g(3) * t95 + t97 * t93;
t1 = [0, -t97, t83, 0, 0, 0, 0, 0, t82, t81, t82, -t83, -t81, -g(1) * t99 - g(2) * (-t96 * pkin(5) + t100) - g(3) * t98, t82, -t81, t83, -g(1) * (pkin(3) * t102 - t94 * qJ(4) + t99) - g(2) * (pkin(3) * t103 + (-pkin(5) + qJ(4)) * t96 + t100) - g(3) * (t93 * pkin(3) + t98);];
U_reg = t1;
