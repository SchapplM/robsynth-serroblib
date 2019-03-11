% Return the minimum parameter vector for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRRRP5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t169 = (m(5) + m(6));
t165 = m(4) + t169;
t181 = (t165 * pkin(9));
t180 = (-Ifges(6,2) - Ifges(7,2));
t173 = (pkin(4) ^ 2);
t179 = (t173 * m(6) + Ifges(5,2));
t178 = 2 * pkin(11) * mrSges(6,3) - t180;
t177 = 2 * pkin(10) * mrSges(5,3) + t179;
t176 = pkin(11) * m(6) + mrSges(6,3);
t175 = pkin(10) * t169 + mrSges(5,3);
t174 = (pkin(3) ^ 2);
t172 = pkin(10) ^ 2;
t171 = pkin(11) ^ 2;
t168 = sin(pkin(7));
t1 = [m(2) + m(3) + t165; pkin(2) ^ 2 * t165 + Ifges(3,3) + (t174 * t169 + Ifges(4,2) + (2 * mrSges(4,3) + t181) * pkin(9)) * t168 ^ 2; pkin(2) * t165 + mrSges(3,1); mrSges(3,2) + (-mrSges(4,3) - t181) * t168; Ifges(4,1) - Ifges(4,2) + (t172 - t174) * t169 + t177; t175 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t172 + t174) * t169 + t177; pkin(3) * t169 + mrSges(4,1); mrSges(4,2) - t175; t171 * m(6) + Ifges(5,1) + t178 - t179; t176 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t171 + t173) * m(6) + t178; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t176; Ifges(6,1) + Ifges(7,1) + t180; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
