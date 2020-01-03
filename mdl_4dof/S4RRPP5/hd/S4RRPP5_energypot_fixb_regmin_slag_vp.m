% Calculate minimal parameter regressor of potential energy for
% S4RRPP5
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:36
% EndTime: 2019-12-31 17:00:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (41->27), mult. (79->28), div. (0->0), fcn. (73->4), ass. (0->14)
t91 = sin(qJ(2));
t93 = cos(qJ(2));
t102 = pkin(2) * t93 + qJ(3) * t91 + pkin(1);
t99 = qJ(4) * t93;
t92 = sin(qJ(1));
t98 = t102 * t92;
t94 = cos(qJ(1));
t97 = t92 * pkin(5) + t102 * t94;
t96 = t91 * pkin(2) - t93 * qJ(3) + pkin(4);
t95 = g(1) * t94 + g(2) * t92;
t80 = g(1) * t92 - g(2) * t94;
t79 = g(3) * t91 + t93 * t95;
t78 = -g(3) * t93 + t95 * t91;
t1 = [0, -t95, t80, 0, 0, 0, 0, 0, -t79, t78, -t80, t79, -t78, -g(1) * t97 - g(2) * (-t94 * pkin(5) + t98) - g(3) * t96, -t80, -t78, -t79, -g(1) * (t92 * pkin(3) + t94 * t99 + t97) - g(2) * (t92 * t99 + (-pkin(3) - pkin(5)) * t94 + t98) - g(3) * (t91 * qJ(4) + t96);];
U_reg = t1;
