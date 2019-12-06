% Calculate minimal parameter regressor of potential energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:14
% EndTime: 2019-12-05 17:38:14
% DurationCPUTime: 0.03s
% Computational Cost: add. (31->20), mult. (48->24), div. (0->0), fcn. (42->6), ass. (0->12)
t95 = sin(qJ(1));
t97 = cos(qJ(1));
t99 = t97 * pkin(1) + t95 * qJ(2);
t98 = t95 * pkin(1) - t97 * qJ(2);
t86 = g(1) * t97 + g(2) * t95;
t96 = cos(qJ(4));
t94 = sin(qJ(4));
t93 = qJ(4) + qJ(5);
t88 = cos(t93);
t87 = sin(t93);
t85 = g(1) * t95 - g(2) * t97;
t1 = [0, -t86, t85, t86, -t85, -g(3) * pkin(5) - g(1) * t99 - g(2) * t98, -t85, -t86, -g(1) * (t97 * qJ(3) + t99) - g(2) * (t95 * qJ(3) + t98) - g(3) * (pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t96 - t86 * t94, g(3) * t94 - t86 * t96, 0, 0, 0, 0, 0, -g(3) * t88 - t86 * t87, g(3) * t87 - t86 * t88;];
U_reg = t1;
