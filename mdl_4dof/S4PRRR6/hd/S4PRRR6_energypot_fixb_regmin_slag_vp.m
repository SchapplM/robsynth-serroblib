% Calculate minimal parameter regressor of potential energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (31->21), mult. (51->39), div. (0->0), fcn. (58->8), ass. (0->16)
t97 = sin(qJ(2));
t105 = g(3) * t97;
t94 = sin(pkin(7));
t99 = cos(qJ(2));
t104 = t94 * t99;
t95 = cos(pkin(7));
t103 = t95 * t99;
t96 = sin(qJ(3));
t102 = t96 * t99;
t98 = cos(qJ(3));
t101 = t98 * t99;
t100 = g(1) * t95 + g(2) * t94;
t93 = qJ(3) + qJ(4);
t92 = cos(t93);
t91 = sin(t93);
t1 = [-g(3) * qJ(1), 0, -t100 * t99 - t105, -g(3) * t99 + t100 * t97, 0, 0, 0, 0, 0, -g(1) * (t95 * t101 + t94 * t96) - g(2) * (t94 * t101 - t95 * t96) - t98 * t105, -g(1) * (-t95 * t102 + t94 * t98) - g(2) * (-t94 * t102 - t95 * t98) + t96 * t105, 0, 0, 0, 0, 0, -g(1) * (t92 * t103 + t94 * t91) - g(2) * (t92 * t104 - t95 * t91) - t92 * t105, -g(1) * (-t91 * t103 + t94 * t92) - g(2) * (-t91 * t104 - t95 * t92) + t91 * t105;];
U_reg = t1;
