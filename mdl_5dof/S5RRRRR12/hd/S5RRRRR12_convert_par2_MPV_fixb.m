% Return the minimum parameter vector for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRR12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t167 = sin(pkin(6));
t170 = (m(5) + m(6));
t164 = (m(4) + t170);
t180 = pkin(9) * t164 + mrSges(4,3);
t189 = t180 * t167;
t169 = cos(pkin(6));
t188 = t180 * t169;
t187 = (pkin(9) * mrSges(4,3));
t159 = m(3) + t164;
t186 = t159 * pkin(8);
t176 = (pkin(3) ^ 2);
t158 = (t176 * t170 + Ifges(4,2));
t175 = (pkin(4) ^ 2);
t185 = (t175 * m(6) + Ifges(5,2));
t184 = 2 * pkin(11) * mrSges(6,3) + Ifges(6,2);
t183 = 2 * pkin(10) * mrSges(5,3) + t185;
t182 = pkin(11) * m(6) + mrSges(6,3);
t181 = t158 + 2 * t187;
t179 = pkin(10) * t170 + mrSges(5,3);
t174 = (pkin(9) ^ 2);
t178 = t174 * t164 + t181;
t177 = pkin(2) ^ 2;
t173 = pkin(10) ^ 2;
t172 = pkin(11) ^ 2;
t168 = sin(pkin(5));
t163 = t169 ^ 2;
t157 = t174 * t163 + t177;
t156 = mrSges(3,3) + t188;
t1 = [pkin(1) ^ 2 * t159 + Ifges(2,3) + (t157 * t164 + Ifges(3,2) + t181 * t163 + (0.2e1 * t156 + t186) * pkin(8)) * t168 ^ 2; pkin(1) * t159 + mrSges(2,1); mrSges(2,2) + (-t156 - t186) * t168; -t163 * t158 + Ifges(3,1) - Ifges(3,2) + (-t157 + t174) * t164 + (-0.2e1 * t163 + 0.2e1) * t187 + t158; pkin(2) * t189 + Ifges(3,4); -pkin(2) * t188 + Ifges(3,5); t178 * t169 * t167 + Ifges(3,6); t178 * t167 ^ 2 + t177 * t164 + Ifges(3,3); pkin(2) * t164 + mrSges(3,1); mrSges(3,2) - t189; t173 * t170 + Ifges(4,1) - t158 + t183; t179 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t173 + t176) * t170 + t183; pkin(3) * t170 + mrSges(4,1); mrSges(4,2) - t179; m(6) * t172 + Ifges(5,1) + t184 - t185; t182 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t172 + t175) * m(6) + t184; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t182; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
