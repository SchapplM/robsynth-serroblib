% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:19:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (56->42), mult. (96->57), div. (0->0), fcn. (98->8), ass. (0->23)
t177 = sin(qJ(3));
t189 = pkin(3) * t177;
t188 = g(3) * t177;
t178 = sin(qJ(1));
t180 = cos(qJ(3));
t187 = t178 * t180;
t175 = qJ(5) + qJ(6);
t170 = sin(t175);
t181 = cos(qJ(1));
t186 = t181 * t170;
t171 = cos(t175);
t185 = t181 * t171;
t176 = sin(qJ(5));
t184 = t181 * t176;
t179 = cos(qJ(5));
t183 = t181 * t179;
t182 = t181 * pkin(1) + t178 * qJ(2);
t168 = g(1) * t178 - g(2) * t181;
t173 = t178 * pkin(1);
t169 = g(1) * t181 + g(2) * t178;
t167 = t168 * t180 - t188;
t166 = g(3) * t180 + t168 * t177;
t1 = [0, -t169, t168, t169, -t168, -g(1) * t182 - g(2) * (-t181 * qJ(2) + t173) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t166, -t167, -t169, t166, t167, -g(1) * (-qJ(4) * t187 + t178 * t189 + t182) - g(2) * (t178 * pkin(7) + t173) - g(3) * (t180 * pkin(3) + t177 * qJ(4) + pkin(2) + pkin(6)) + (-g(1) * pkin(7) - g(2) * (qJ(4) * t180 - qJ(2) - t189)) * t181, 0, 0, 0, 0, 0, -g(1) * (-t176 * t187 + t183) - g(2) * (t178 * t179 + t180 * t184) - t176 * t188, -g(1) * (-t179 * t187 - t184) - g(2) * (-t178 * t176 + t180 * t183) - t179 * t188, 0, 0, 0, 0, 0, -g(1) * (-t170 * t187 + t185) - g(2) * (t178 * t171 + t180 * t186) - t170 * t188, -g(1) * (-t171 * t187 - t186) - g(2) * (-t178 * t170 + t180 * t185) - t171 * t188;];
U_reg  = t1;
