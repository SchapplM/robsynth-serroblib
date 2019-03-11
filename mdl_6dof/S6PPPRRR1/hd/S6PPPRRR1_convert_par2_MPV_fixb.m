% Return the minimum parameter vector for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% MPV [20x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PPPRRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t171 = (m(6) + m(7));
t175 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t174 = pkin(11) * m(7) + mrSges(7,3);
t173 = (pkin(5) ^ 2);
t172 = pkin(11) ^ 2;
t1 = [m(2); m(3); m(4) + m(5) + t171; Ifges(5,3) + t173 * m(7) + Ifges(6,2) + 2 * pkin(10) * mrSges(6,3) + (pkin(4) ^ 2 + pkin(10) ^ 2) * t171; pkin(4) * t171 + mrSges(5,1); -pkin(10) * t171 + mrSges(5,2) - mrSges(6,3); Ifges(6,1) - Ifges(6,2) + (t172 - t173) * m(7) + t175; t174 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t172 + t173) * m(7) + t175; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t174; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
