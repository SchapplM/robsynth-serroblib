% Calculate minimal parameter regressor of potential energy for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:42
% EndTime: 2019-12-31 16:30:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (43->28), mult. (93->42), div. (0->0), fcn. (101->6), ass. (0->19)
t92 = sin(qJ(2));
t100 = g(3) * t92;
t91 = sin(qJ(3));
t94 = cos(qJ(2));
t99 = t91 * t94;
t93 = cos(qJ(3));
t98 = t93 * t94;
t89 = sin(pkin(6));
t90 = cos(pkin(6));
t97 = g(1) * t90 + g(2) * t89;
t96 = pkin(2) * t94 + pkin(5) * t92 + pkin(1);
t84 = t89 * t99 + t90 * t93;
t86 = -t89 * t93 + t90 * t99;
t95 = g(1) * t86 + g(2) * t84 + t91 * t100;
t87 = t89 * t91 + t90 * t98;
t85 = t89 * t98 - t90 * t91;
t83 = -g(3) * t94 + t97 * t92;
t82 = -g(1) * t87 - g(2) * t85 - t93 * t100;
t1 = [-g(3) * qJ(1), 0, -t97 * t94 - t100, t83, 0, 0, 0, 0, 0, t82, t95, t82, -t83, -t95, -g(1) * (t87 * pkin(3) + t86 * qJ(4)) - g(2) * (t85 * pkin(3) + t84 * qJ(4)) - g(3) * (-t94 * pkin(5) + qJ(1)) - (pkin(3) * t93 + qJ(4) * t91 + pkin(2)) * t100 + (g(2) * pkin(4) - g(1) * t96) * t90 + (-g(1) * pkin(4) - g(2) * t96) * t89;];
U_reg = t1;
