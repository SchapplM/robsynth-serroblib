% Calculate minimal parameter regressor of potential energy for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:14
% EndTime: 2019-12-31 18:30:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (85->42), mult. (98->60), div. (0->0), fcn. (100->10), ass. (0->27)
t169 = pkin(8) + qJ(3);
t165 = sin(t169);
t187 = g(3) * t165;
t168 = pkin(9) + qJ(5);
t164 = sin(t168);
t175 = sin(qJ(1));
t186 = t175 * t164;
t166 = cos(t168);
t185 = t175 * t166;
t170 = sin(pkin(9));
t184 = t175 * t170;
t172 = cos(pkin(9));
t183 = t175 * t172;
t176 = cos(qJ(1));
t182 = t176 * t164;
t181 = t176 * t166;
t180 = t176 * t170;
t179 = t176 * t172;
t178 = g(1) * t176 + g(2) * t175;
t167 = cos(t169);
t173 = cos(pkin(8));
t177 = t173 * pkin(2) + pkin(3) * t167 + qJ(4) * t165 + pkin(1);
t174 = -pkin(6) - qJ(2);
t171 = sin(pkin(8));
t162 = g(1) * t175 - g(2) * t176;
t161 = -g(3) * t167 + t178 * t165;
t1 = [0, -t178, t162, -g(3) * t171 - t178 * t173, -g(3) * t173 + t178 * t171, -t162, -g(1) * (t176 * pkin(1) + t175 * qJ(2)) - g(2) * (t175 * pkin(1) - t176 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t178 * t167 - t187, t161, -g(1) * (t167 * t179 + t184) - g(2) * (t167 * t183 - t180) - t172 * t187, -g(1) * (-t167 * t180 + t183) - g(2) * (-t167 * t184 - t179) + t170 * t187, -t161, -g(3) * (t171 * pkin(2) + t165 * pkin(3) - t167 * qJ(4) + pkin(5)) + (-g(1) * t177 - g(2) * t174) * t176 + (g(1) * t174 - g(2) * t177) * t175, 0, 0, 0, 0, 0, -g(1) * (t167 * t181 + t186) - g(2) * (t167 * t185 - t182) - t166 * t187, -g(1) * (-t167 * t182 + t185) - g(2) * (-t167 * t186 - t181) + t164 * t187;];
U_reg = t1;
