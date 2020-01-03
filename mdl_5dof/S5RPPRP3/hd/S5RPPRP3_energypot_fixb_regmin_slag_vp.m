% Calculate minimal parameter regressor of potential energy for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:01
% EndTime: 2019-12-31 17:51:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (56->26), mult. (50->31), div. (0->0), fcn. (41->6), ass. (0->17)
t101 = sin(qJ(4));
t110 = pkin(4) * t101;
t100 = qJ(2) + pkin(5);
t109 = g(3) * t100;
t102 = sin(qJ(1));
t98 = qJ(1) + pkin(7);
t94 = sin(t98);
t108 = t102 * pkin(1) + t94 * pkin(2);
t104 = cos(qJ(1));
t95 = cos(t98);
t107 = t104 * pkin(1) + t95 * pkin(2) + t94 * qJ(3);
t106 = -g(1) * t94 + g(2) * t95;
t105 = -g(1) * t104 - g(2) * t102;
t103 = cos(qJ(4));
t99 = -qJ(5) - pkin(6);
t90 = g(1) * t95 + g(2) * t94;
t1 = [0, t105, g(1) * t102 - g(2) * t104, t105 * pkin(1) - t109, t90, t106, -g(1) * t107 - g(2) * (-t95 * qJ(3) + t108) - t109, 0, 0, 0, 0, 0, -g(3) * t103 + t106 * t101, g(3) * t101 + t106 * t103, -t90, -g(1) * (t94 * t110 - t95 * t99 + t107) - g(2) * (-t94 * t99 + (-qJ(3) - t110) * t95 + t108) - g(3) * (t103 * pkin(4) + pkin(3) + t100);];
U_reg = t1;
