% Return the minimum parameter vector for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRP10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t164 = (m(5) + m(6));
t160 = (m(4) + t164);
t156 = m(3) + t160;
t179 = (t156 * pkin(8));
t178 = (-Ifges(6,2) - Ifges(7,3));
t169 = (pkin(3) ^ 2);
t177 = (t169 * t164 + Ifges(4,2));
t165 = (pkin(11) ^ 2);
t168 = (pkin(4) ^ 2);
t176 = (Ifges(5,2) + (t165 + t168) * m(6));
t175 = 2 * pkin(9) * mrSges(4,3) + t177;
t174 = -pkin(11) * m(6) - mrSges(6,3);
t173 = pkin(9) * t160 + mrSges(4,3);
t155 = (mrSges(5,3) - t174);
t172 = pkin(10) * t164 + t155;
t171 = 2 * pkin(11) * mrSges(6,3) + 2 * pkin(10) * t155 + t176 - t178;
t170 = (pkin(2) ^ 2);
t167 = pkin(9) ^ 2;
t166 = pkin(10) ^ 2;
t163 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * t156 + Ifges(2,3) + (t170 * t160 + Ifges(3,2) + (2 * mrSges(3,3) + t179) * pkin(8)) * t163 ^ 2; pkin(1) * t156 + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t179) * t163; Ifges(3,1) - Ifges(3,2) + (t167 - t170) * t160 + t175; t173 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t167 + t170) * t160 + t175; pkin(2) * t160 + mrSges(3,1); mrSges(3,2) - t173; t166 * t164 + Ifges(4,1) + t171 - t177; pkin(3) * t172 + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t166 + t169) * t164 + t171; pkin(3) * t164 + mrSges(4,1); mrSges(4,2) - t172; t165 * m(6) + Ifges(5,1) - t176; Ifges(5,4); t174 * pkin(4) + Ifges(5,5); Ifges(5,6); t168 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t178; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
