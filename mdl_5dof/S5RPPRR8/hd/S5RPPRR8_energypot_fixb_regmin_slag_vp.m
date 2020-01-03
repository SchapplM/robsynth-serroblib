% Calculate minimal parameter regressor of potential energy for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:13
% EndTime: 2019-12-31 18:01:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (51->23), mult. (66->36), div. (0->0), fcn. (72->8), ass. (0->19)
t107 = sin(qJ(1));
t109 = cos(qJ(1));
t115 = t109 * pkin(1) + t107 * qJ(2);
t114 = pkin(8) + qJ(4);
t113 = cos(t114);
t112 = sin(t114);
t111 = t107 * pkin(1) - t109 * qJ(2);
t93 = -t107 * t112 - t109 * t113;
t94 = t107 * t113 - t109 * t112;
t110 = g(1) * t93 - g(2) * t94;
t108 = cos(qJ(5));
t106 = sin(qJ(5));
t105 = cos(pkin(8));
t104 = sin(pkin(8));
t98 = -g(1) * t109 - g(2) * t107;
t97 = g(1) * t107 - g(2) * t109;
t96 = -t109 * t104 + t107 * t105;
t95 = -t107 * t104 - t109 * t105;
t1 = [0, t98, t97, t98, -t97, -g(3) * pkin(5) - g(1) * t115 - g(2) * t111, g(1) * t95 - g(2) * t96, -g(1) * t96 - g(2) * t95, -g(1) * (t109 * pkin(2) + t115) - g(2) * (t107 * pkin(2) + t111) - g(3) * (-qJ(3) + pkin(5)), 0, t110, -g(1) * t94 - g(2) * t93, 0, 0, 0, 0, 0, g(3) * t106 + t110 * t108, g(3) * t108 - t110 * t106;];
U_reg = t1;
