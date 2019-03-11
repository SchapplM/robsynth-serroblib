% Return the minimum parameter vector for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t122 = (pkin(8) * mrSges(6,3));
t121 = -Ifges(6,2) - Ifges(7,3);
t120 = -pkin(8) * m(6) - mrSges(6,3);
t119 = -Ifges(4,2) - Ifges(5,3) + t121;
t110 = -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3);
t111 = m(4) * pkin(2) + mrSges(3,1);
t114 = sin(pkin(9));
t115 = cos(pkin(9));
t118 = -t114 * t110 + t115 * t111;
t117 = pkin(4) ^ 2;
t116 = pkin(8) ^ 2;
t113 = 2 * t122;
t112 = t116 + t117;
t1 = [Ifges(2,3) + Ifges(3,3) + t113 + t112 * m(6) + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)) + 0.2e1 * t118 * pkin(1) - t119; mrSges(2,1) + t118; t115 * t110 + t114 * t111 + mrSges(2,2); m(3) + m(4); -(2 * t122) + Ifges(4,1) + Ifges(5,2) + (-t112 + t117) * m(6) + t119; Ifges(4,4) + Ifges(5,6); t120 * pkin(4) - Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,5); t116 * m(6) + Ifges(5,1) + Ifges(4,3) + t113 - t121; mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) + t120; mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t121; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
