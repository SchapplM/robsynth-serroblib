% Calculate minimal parameter regressor of potential energy for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:13
% EndTime: 2019-12-05 18:12:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->19), mult. (53->27), div. (0->0), fcn. (50->10), ass. (0->16)
t162 = pkin(9) + qJ(3);
t161 = qJ(4) + t162;
t165 = sin(qJ(1));
t166 = cos(qJ(1));
t167 = g(1) * t166 + g(2) * t165;
t164 = cos(pkin(9));
t163 = sin(pkin(9));
t160 = cos(t162);
t159 = sin(t162);
t158 = qJ(5) + t161;
t157 = cos(t161);
t156 = sin(t161);
t155 = cos(t158);
t154 = sin(t158);
t153 = g(1) * t165 - g(2) * t166;
t1 = [0, -t167, t153, -g(3) * t163 - t167 * t164, -g(3) * t164 + t167 * t163, -t153, -g(1) * (t166 * pkin(1) + t165 * qJ(2)) - g(2) * (t165 * pkin(1) - t166 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t159 - t167 * t160, -g(3) * t160 + t167 * t159, 0, 0, 0, 0, 0, -g(3) * t156 - t167 * t157, -g(3) * t157 + t167 * t156, 0, 0, 0, 0, 0, -g(3) * t154 - t167 * t155, -g(3) * t155 + t167 * t154;];
U_reg = t1;
