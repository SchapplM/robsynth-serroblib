% Calculate minimal parameter regressor of potential energy for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:11
% EndTime: 2019-12-05 16:40:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->24), mult. (35->26), div. (0->0), fcn. (31->8), ass. (0->13)
t92 = pkin(8) + qJ(2);
t91 = qJ(3) + t92;
t86 = sin(t91);
t87 = cos(t91);
t96 = g(1) * t87 + g(2) * t86;
t95 = cos(qJ(4));
t94 = sin(qJ(4));
t93 = -qJ(5) - pkin(7);
t90 = cos(t92);
t89 = sin(t92);
t88 = t95 * pkin(4) + pkin(3);
t85 = g(1) * t86 - g(2) * t87;
t1 = [-g(3) * qJ(1), 0, -g(1) * t90 - g(2) * t89, g(1) * t89 - g(2) * t90, 0, -t96, t85, 0, 0, 0, 0, 0, -g(3) * t94 - t96 * t95, -g(3) * t95 + t96 * t94, -t85, -g(1) * (t87 * t88 - t86 * t93 + pkin(2) * t90 + cos(pkin(8)) * pkin(1)) - g(2) * (t86 * t88 + t87 * t93 + pkin(2) * t89 + sin(pkin(8)) * pkin(1)) - g(3) * (t94 * pkin(4) + pkin(5) + pkin(6) + qJ(1));];
U_reg = t1;
