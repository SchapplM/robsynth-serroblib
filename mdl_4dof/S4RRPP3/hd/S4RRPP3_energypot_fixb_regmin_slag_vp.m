% Calculate minimal parameter regressor of potential energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:50
% EndTime: 2019-12-31 16:57:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (48->24), mult. (58->30), div. (0->0), fcn. (52->6), ass. (0->16)
t111 = sin(qJ(2));
t119 = t111 * pkin(2) + pkin(4);
t113 = cos(qJ(2));
t104 = t113 * pkin(2) + pkin(1);
t110 = -pkin(5) - qJ(3);
t112 = sin(qJ(1));
t114 = cos(qJ(1));
t118 = t112 * t104 + t114 * t110;
t117 = t114 * t104 - t112 * t110;
t116 = g(1) * t114 + g(2) * t112;
t109 = qJ(2) + pkin(6);
t105 = sin(t109);
t106 = cos(t109);
t115 = pkin(3) * t106 + qJ(4) * t105;
t100 = g(1) * t112 - g(2) * t114;
t1 = [0, -t116, t100, 0, 0, 0, 0, 0, -g(3) * t111 - t116 * t113, -g(3) * t113 + t116 * t111, -t100, -g(1) * t117 - g(2) * t118 - g(3) * t119, -g(3) * t105 - t116 * t106, -t100, g(3) * t106 - t116 * t105, -g(1) * (t115 * t114 + t117) - g(2) * (t115 * t112 + t118) - g(3) * (t105 * pkin(3) - t106 * qJ(4) + t119);];
U_reg = t1;
