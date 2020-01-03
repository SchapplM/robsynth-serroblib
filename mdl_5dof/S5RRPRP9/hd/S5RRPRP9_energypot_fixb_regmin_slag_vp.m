% Calculate minimal parameter regressor of potential energy for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:27
% EndTime: 2019-12-31 20:07:27
% DurationCPUTime: 0.11s
% Computational Cost: add. (102->51), mult. (145->69), div. (0->0), fcn. (155->8), ass. (0->30)
t174 = sin(qJ(2));
t189 = g(3) * t174;
t171 = sin(pkin(8));
t175 = sin(qJ(1));
t188 = t175 * t171;
t176 = cos(qJ(2));
t187 = t175 * t176;
t170 = pkin(8) + qJ(4);
t164 = sin(t170);
t177 = cos(qJ(1));
t186 = t177 * t164;
t165 = cos(t170);
t185 = t177 * t165;
t184 = t177 * t171;
t172 = cos(pkin(8));
t183 = t177 * t172;
t182 = t177 * pkin(1) + t175 * pkin(6);
t181 = t175 * pkin(1) - t177 * pkin(6);
t180 = g(1) * t177 + g(2) * t175;
t179 = pkin(2) * t176 + qJ(3) * t174;
t157 = t164 * t187 + t185;
t159 = -t175 * t165 + t176 * t186;
t178 = g(1) * t159 + g(2) * t157 + t164 * t189;
t173 = -pkin(7) - qJ(3);
t163 = t172 * pkin(3) + pkin(2);
t161 = -g(3) * t176 + t180 * t174;
t160 = t175 * t164 + t176 * t185;
t158 = t165 * t187 - t186;
t156 = -g(1) * t160 - g(2) * t158 - t165 * t189;
t1 = [0, -t180, g(1) * t175 - g(2) * t177, 0, 0, 0, 0, 0, -t180 * t176 - t189, t161, -g(1) * (t176 * t183 + t188) - g(2) * (t172 * t187 - t184) - t172 * t189, -g(1) * (t175 * t172 - t176 * t184) - g(2) * (-t171 * t187 - t183) + t171 * t189, -t161, -g(1) * (t179 * t177 + t182) - g(2) * (t179 * t175 + t181) - g(3) * (t174 * pkin(2) - t176 * qJ(3) + pkin(5)), 0, 0, 0, 0, 0, t156, t178, t156, -t161, -t178, -g(1) * (t177 * t176 * t163 + pkin(3) * t188 + t160 * pkin(4) + t159 * qJ(5) + t182) - g(2) * (-pkin(3) * t184 + t158 * pkin(4) + t157 * qJ(5) + t163 * t187 + t181) - g(3) * (t176 * t173 + pkin(5)) + (-g(3) * (pkin(4) * t165 + qJ(5) * t164 + t163) + t180 * t173) * t174;];
U_reg = t1;
