% Return the minimum parameter vector for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPPR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t128 = (pkin(9) * mrSges(7,3));
t127 = -pkin(8) * m(4) - mrSges(4,3);
t126 = pkin(9) * m(7) + mrSges(7,3);
t125 = (-Ifges(6,1) - Ifges(4,2) - Ifges(7,2) - Ifges(5,3));
t124 = (mrSges(3,3) - t127);
t123 = (pkin(2) ^ 2);
t122 = (pkin(5) ^ 2);
t121 = (pkin(8) ^ 2);
t120 = (pkin(9) ^ 2);
t119 = (m(3) + m(4));
t118 = 2 * t128;
t117 = (t121 + t123);
t1 = [Ifges(2,3) + Ifges(3,2) + m(7) * t120 + t118 + 2 * pkin(8) * mrSges(4,3) + t117 * m(4) + 2 * pkin(7) * t124 + (pkin(1) ^ 2 + pkin(7) ^ 2) * t119 - t125; pkin(1) * t119 + mrSges(2,1); -pkin(7) * t119 + mrSges(2,2) - t124; Ifges(3,1) - Ifges(3,2) + (-t117 + t121) * m(4); Ifges(3,4); t127 * pkin(2) + Ifges(3,5); Ifges(3,6); t123 * m(4) + Ifges(3,3); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); -2 * t128 + Ifges(4,1) + Ifges(5,1) + Ifges(6,2) + (-t120 + t122) * m(7) + t125; t126 * pkin(5) + Ifges(4,4) + Ifges(6,4) - Ifges(5,5); Ifges(4,5) + Ifges(5,4) + Ifges(6,6); Ifges(4,6) - Ifges(5,6) + Ifges(6,5); Ifges(4,3) + Ifges(5,2) + Ifges(6,3) + Ifges(7,2) + t118 + (t120 + t122) * m(7); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t126; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
