% Return the minimum parameter vector for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t112 = (m(6) + m(7));
t117 = (pkin(4) ^ 2);
t123 = (t117 * t112 + Ifges(5,2));
t113 = (pkin(9) ^ 2);
t116 = (pkin(5) ^ 2);
t122 = (Ifges(6,2) + (t113 + t116) * m(7));
t107 = (m(5) + t112);
t121 = 2 * pkin(7) * mrSges(5,3) + t123;
t120 = -pkin(9) * m(7) - mrSges(7,3);
t104 = (mrSges(6,3) - t120);
t119 = pkin(8) * t112 + t104;
t118 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t104 + Ifges(7,2) + t122;
t115 = pkin(7) ^ 2;
t114 = pkin(8) ^ 2;
t111 = cos(pkin(10));
t110 = sin(pkin(10));
t1 = [Ifges(2,3) + Ifges(3,1) + t111 ^ 2 * (t115 * t107 + Ifges(4,1) + t121) + (-0.2e1 * t111 * Ifges(4,4) + (Ifges(4,2) + (pkin(3) ^ 2 + t115) * t107 + t121) * t110) * t110; mrSges(2,1); mrSges(2,2); mrSges(3,2); mrSges(3,3); m(3); pkin(3) * t107 + mrSges(4,1); mrSges(4,2); pkin(7) * t107 + mrSges(4,3) + mrSges(5,3); m(4) + t107; t114 * t112 + Ifges(5,1) + t118 - t123; t119 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t114 + t117) * t112 + t118; pkin(4) * t112 + mrSges(5,1); mrSges(5,2) - t119; m(7) * t113 + Ifges(6,1) - t122; Ifges(6,4); t120 * pkin(5) + Ifges(6,5); Ifges(6,6); t116 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
