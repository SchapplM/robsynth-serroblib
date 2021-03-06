% Return the minimum parameter vector for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRRR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t114 = (m(6) + m(7));
t128 = -pkin(8) * t114 - mrSges(6,3);
t117 = (pkin(8) ^ 2);
t120 = (pkin(4) ^ 2);
t127 = (Ifges(5,2) + (t117 + t120) * t114);
t119 = (pkin(5) ^ 2);
t126 = (t119 * m(7) + Ifges(6,2));
t125 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t107 = (m(5) + t114);
t124 = (mrSges(5,3) - t128);
t123 = pkin(9) * m(7) + mrSges(7,3);
t111 = sin(pkin(10));
t113 = cos(pkin(10));
t122 = t113 * mrSges(3,1) - t111 * mrSges(3,2);
t121 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(7) * t124 + t126 + t127;
t118 = pkin(7) ^ 2;
t116 = pkin(9) ^ 2;
t112 = cos(pkin(11));
t110 = sin(pkin(11));
t1 = [Ifges(2,3) + Ifges(3,3) + t112 ^ 2 * (Ifges(4,2) + (pkin(3) ^ 2 + t118) * t107 + t121) + (0.2e1 * t112 * Ifges(4,4) + (t118 * t107 + Ifges(4,1) + t121) * t110) * t110 + 0.2e1 * t122 * pkin(1); mrSges(2,1) + t122; t111 * mrSges(3,1) + t113 * mrSges(3,2) + mrSges(2,2); m(3); pkin(3) * t107 + mrSges(4,1); mrSges(4,2); pkin(7) * t107 + mrSges(4,3) + t124; m(4) + t107; t117 * t114 + Ifges(5,1) - t127; Ifges(5,4); t128 * pkin(4) + Ifges(5,5); Ifges(5,6); t120 * t114 + Ifges(5,3); pkin(4) * t114 + mrSges(5,1); mrSges(5,2); m(7) * t116 + Ifges(6,1) + t125 - t126; t123 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t116 + t119) * m(7) + t125; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t123; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
