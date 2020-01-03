% Calculate minimal parameter regressor of potential energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:41
% DurationCPUTime: 0.07s
% Computational Cost: add. (42->24), mult. (68->30), div. (0->0), fcn. (75->6), ass. (0->16)
t102 = cos(qJ(3));
t101 = sin(qJ(3));
t100 = g(3) * qJ(1);
t91 = cos(pkin(7));
t98 = sin(pkin(7));
t99 = t91 * pkin(1) + t98 * qJ(2);
t97 = t98 * pkin(1) - t91 * qJ(2);
t82 = -t98 * t101 - t91 * t102;
t83 = t91 * t101 - t98 * t102;
t96 = g(1) * t83 - g(2) * t82;
t95 = g(1) * t82 + g(2) * t83;
t94 = cos(qJ(4));
t93 = sin(qJ(4));
t92 = -qJ(5) - pkin(6);
t87 = t94 * pkin(4) + pkin(3);
t1 = [-t100, -g(1) * t99 - g(2) * t97 - t100, 0, t95, t96, 0, 0, 0, 0, 0, g(3) * t93 + t95 * t94, g(3) * t94 - t95 * t93, -t96, -g(1) * (t91 * pkin(2) - t82 * t87 - t83 * t92 + t99) - g(2) * (t98 * pkin(2) + t82 * t92 - t83 * t87 + t97) - g(3) * (-t93 * pkin(4) - pkin(5) + qJ(1));];
U_reg = t1;
