% Return the minimum parameter vector for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% MPV [29x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t130 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t122 = sin(pkin(9));
t123 = cos(pkin(9));
t129 = t122 * t123;
t128 = Ifges(6,4) * t129;
t127 = -pkin(8) * m(7) - mrSges(7,3);
t124 = pkin(8) ^ 2;
t125 = (pkin(5) ^ 2);
t115 = Ifges(6,2) + (t124 + t125) * m(7) + t130;
t117 = m(7) * t124 + Ifges(6,1) + t130;
t119 = t122 ^ 2;
t120 = t123 ^ 2;
t126 = t119 * t115 + t120 * t117 + Ifges(5,1) + Ifges(3,2) + Ifges(4,3);
t118 = t127 * pkin(5) + Ifges(6,5);
t116 = 0.2e1 * t128;
t1 = [Ifges(2,3) - 0.2e1 * t128 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)) + t126; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); (t125 * m(7)) + Ifges(3,1) + Ifges(4,1) + Ifges(5,2) + Ifges(6,3) + t116 - t126; t122 * Ifges(6,6) - t123 * t118 + Ifges(3,4) + Ifges(5,4) - Ifges(4,5); -t123 * Ifges(6,6) - t122 * t118 + Ifges(4,4) + Ifges(3,5) + Ifges(5,6); Ifges(3,6) - Ifges(4,6) + Ifges(5,5) + (t120 - t119) * Ifges(6,4) + (-t115 + t117) * t129; t120 * t115 + t119 * t117 + Ifges(4,2) + Ifges(3,3) + Ifges(5,3) + t116; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t127; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
