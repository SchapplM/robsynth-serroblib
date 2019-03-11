% Return the minimum parameter vector for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t113 = (pkin(8) ^ 2);
t114 = (pkin(5) ^ 2);
t120 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t103 = Ifges(5,2) + Ifges(6,3) + (t113 + t114) * m(7) + t120;
t105 = t114 * m(7) + Ifges(5,1) + Ifges(6,2);
t122 = -t103 + t105;
t121 = (pkin(7) * m(4));
t109 = sin(pkin(9));
t110 = cos(pkin(9));
t119 = t109 * t110;
t106 = t109 ^ 2;
t107 = t110 ^ 2;
t118 = t107 - t106;
t117 = -pkin(8) * m(7) - mrSges(7,3);
t112 = Ifges(5,4) + Ifges(6,6);
t116 = t112 * t119;
t115 = t110 * mrSges(5,1) - t109 * mrSges(5,2);
t111 = Ifges(5,6) - Ifges(6,5);
t104 = t117 * pkin(5) - Ifges(6,4) + Ifges(5,5);
t1 = [0.2e1 * t116 + t107 * t103 + t106 * t105 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t121) * pkin(7)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t121; mrSges(3,3); m(3) + m(4); t122 * t118 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t116; t118 * t112 + t122 * t119 + Ifges(4,4); t110 * t104 - t109 * t111 + Ifges(4,5); t109 * t104 + t110 * t111 + Ifges(4,6); (m(7) * t113) + 0.2e1 * pkin(3) * t115 + Ifges(6,1) + Ifges(4,3) + Ifges(5,3) + t120; mrSges(4,1) + t115; t109 * mrSges(5,1) + t110 * mrSges(5,2) + mrSges(4,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) + t117; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
