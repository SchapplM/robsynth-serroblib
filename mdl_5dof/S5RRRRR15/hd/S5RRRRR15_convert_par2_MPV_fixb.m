% Return the minimum parameter vector for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRR15_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR15_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t171 = cos(pkin(6));
t182 = m(6) * pkin(11) + mrSges(6,3);
t188 = t182 * t171;
t169 = sin(pkin(6));
t187 = t182 * t169;
t157 = mrSges(5,3) + t188;
t172 = m(5) + m(6);
t186 = -pkin(10) * t172 - t157;
t184 = (mrSges(6,3) * pkin(11));
t168 = m(4) + t172;
t161 = m(3) + t168;
t183 = t161 * pkin(8);
t181 = Ifges(6,2) + 2 * t184;
t156 = mrSges(4,3) - t186;
t180 = -pkin(9) * t168 - t156;
t173 = (pkin(11) ^ 2);
t179 = m(6) * t173 + t181;
t178 = pkin(2) ^ 2;
t177 = pkin(3) ^ 2;
t176 = pkin(4) ^ 2;
t175 = pkin(9) ^ 2;
t174 = pkin(10) ^ 2;
t170 = sin(pkin(5));
t167 = t171 ^ 2;
t164 = t175 + t178;
t163 = t174 + t177;
t158 = t173 * t167 + t176;
t155 = mrSges(3,3) - t180;
t1 = [pkin(1) ^ 2 * t161 + Ifges(2,3) + (t158 * m(6) + 0.2e1 * pkin(9) * t156 + 0.2e1 * pkin(10) * t157 + t163 * t172 + t164 * t168 + Ifges(3,2) + Ifges(4,2) + Ifges(5,2) + t181 * t167 + (0.2e1 * t155 + t183) * pkin(8)) * t170 ^ 2; pkin(1) * t161 + mrSges(2,1); mrSges(2,2) + (-t155 - t183) * t170; Ifges(3,1) - Ifges(3,2) + (-t164 + t175) * t168; Ifges(3,4); t180 * pkin(2) + Ifges(3,5); Ifges(3,6); t178 * t168 + Ifges(3,3); pkin(2) * t168 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (-t163 + t174) * t172; Ifges(4,4); t186 * pkin(3) + Ifges(4,5); Ifges(4,6); t177 * t172 + Ifges(4,3); pkin(3) * t172 + mrSges(4,1); mrSges(4,2); -t167 * Ifges(6,2) + Ifges(5,1) - Ifges(5,2) + Ifges(6,2) + (-0.2e1 * t167 + 0.2e1) * t184 + (-t158 + t173) * m(6); pkin(4) * t187 + Ifges(5,4); -pkin(4) * t188 + Ifges(5,5); t179 * t171 * t169 + Ifges(5,6); t179 * t169 ^ 2 + t176 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t187; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
