% Calculate minimal parameter regressor of potential energy for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:55
% EndTime: 2019-12-31 19:52:55
% DurationCPUTime: 0.05s
% Computational Cost: add. (70->29), mult. (64->33), div. (0->0), fcn. (58->6), ass. (0->16)
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t110 = pkin(4) * t101 - qJ(5) * t103;
t109 = pkin(6) + pkin(5);
t102 = sin(qJ(1));
t100 = qJ(1) + qJ(2);
t96 = sin(t100);
t107 = t102 * pkin(1) + t96 * pkin(2);
t104 = cos(qJ(1));
t97 = cos(t100);
t105 = t104 * pkin(1) + t97 * pkin(2) + t96 * qJ(3);
t91 = g(1) * t96 - g(2) * t97;
t92 = g(1) * t97 + g(2) * t96;
t90 = -g(3) * t101 + t91 * t103;
t89 = -g(3) * t103 - t91 * t101;
t1 = [0, -g(1) * t104 - g(2) * t102, g(1) * t102 - g(2) * t104, 0, -t92, t91, t92, -t91, -g(1) * t105 - g(2) * (-t97 * qJ(3) + t107) - g(3) * t109, 0, 0, 0, 0, 0, t89, -t90, t89, -t92, t90, -g(1) * (t110 * t96 + t105) - g(2) * (t96 * pkin(7) + t107) - g(3) * (t103 * pkin(4) + t101 * qJ(5) + pkin(3) + t109) + (-g(1) * pkin(7) - g(2) * (-qJ(3) - t110)) * t97;];
U_reg = t1;
