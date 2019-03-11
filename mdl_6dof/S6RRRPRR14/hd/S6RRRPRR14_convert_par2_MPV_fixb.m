% Return the minimum parameter vector for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR14_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR14_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t168 = -pkin(11) * m(7) - mrSges(7,3);
t150 = mrSges(6,3) - t168;
t173 = (pkin(10) * t150);
t174 = (pkin(11) * mrSges(7,3));
t175 = 2 * t173 + 2 * t174;
t157 = (m(6) + m(7));
t158 = m(3) + m(4);
t172 = (t158 * pkin(8));
t159 = pkin(11) ^ 2;
t162 = pkin(5) ^ 2;
t171 = Ifges(6,2) + (t159 + t162) * m(7);
t170 = (Ifges(7,2) + t171);
t169 = pkin(9) * m(4) + mrSges(4,3);
t167 = -pkin(10) * t157 - t150;
t160 = (pkin(10) ^ 2);
t163 = (pkin(4) ^ 2);
t166 = (Ifges(4,2) + Ifges(5,3) + (t160 + t163) * t157 + t170);
t165 = 2 * pkin(9) * mrSges(4,3) + t166 + t175;
t164 = (pkin(2) ^ 2);
t161 = pkin(9) ^ 2;
t156 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t158 + Ifges(2,3) + (t164 * m(4) + Ifges(3,2) + (2 * mrSges(3,3) + t172) * pkin(8)) * t156 ^ 2; pkin(1) * t158 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t172) * t156; Ifges(3,1) - Ifges(3,2) + (t161 - t164) * m(4) + t165; pkin(2) * t169 + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t161 + t164) * m(4) + t165; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t169; t163 * t157 + Ifges(4,1) + Ifges(5,2) - t166 - 2 * t173 - 2 * t174; Ifges(4,4) + Ifges(5,6); pkin(4) * t167 - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5); t160 * t157 + Ifges(5,1) + Ifges(4,3) + t170 + t175; mrSges(4,1); mrSges(4,2); pkin(4) * t157 + mrSges(5,1); mrSges(5,2) + t167; mrSges(5,3); m(5) + t157; m(7) * t159 + Ifges(6,1) - t171; Ifges(6,4); pkin(5) * t168 + Ifges(6,5); Ifges(6,6); t162 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
