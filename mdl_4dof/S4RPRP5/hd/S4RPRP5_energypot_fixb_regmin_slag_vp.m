% Calculate minimal parameter regressor of potential energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (53->26), mult. (65->32), div. (0->0), fcn. (59->6), ass. (0->14)
t100 = cos(qJ(1));
t99 = sin(qJ(1));
t102 = g(1) * t100 + g(2) * t99;
t95 = pkin(6) + qJ(3);
t92 = sin(t95);
t93 = cos(t95);
t97 = cos(pkin(6));
t101 = t97 * pkin(2) + pkin(3) * t93 + qJ(4) * t92 + pkin(1);
t98 = -pkin(5) - qJ(2);
t96 = sin(pkin(6));
t90 = g(1) * t99 - g(2) * t100;
t89 = -g(3) * t92 - t102 * t93;
t88 = -g(3) * t93 + t102 * t92;
t1 = [0, -t102, t90, -g(3) * t96 - t102 * t97, -g(3) * t97 + t102 * t96, -t90, -g(1) * (t100 * pkin(1) + t99 * qJ(2)) - g(2) * (t99 * pkin(1) - t100 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, t89, t88, t89, -t90, -t88, -g(3) * (t96 * pkin(2) + t92 * pkin(3) - t93 * qJ(4) + pkin(4)) + (g(1) * t98 - g(2) * t101) * t99 + (-g(1) * t101 - g(2) * t98) * t100;];
U_reg = t1;
