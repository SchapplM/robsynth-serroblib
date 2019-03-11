% Return the minimum parameter vector for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPPRPR4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t112 = (pkin(8) ^ 2);
t119 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t102 = m(7) * t112 + Ifges(6,1) + t119;
t113 = (pkin(5) ^ 2);
t105 = t113 * m(7) + Ifges(6,2);
t120 = t102 - t105;
t110 = sin(pkin(10));
t111 = cos(pkin(10));
t118 = t110 * t111;
t107 = t110 ^ 2;
t108 = t111 ^ 2;
t117 = t108 - t107;
t116 = pkin(8) * m(7) + mrSges(7,3);
t103 = t116 * pkin(5) + Ifges(6,4);
t115 = t103 * t118;
t104 = mrSges(6,2) - t116;
t106 = m(7) * pkin(5) + mrSges(6,1);
t114 = -t110 * t104 + t111 * t106;
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,3) + Ifges(5,2) + t107 * t102 + 0.2e1 * t115 + t108 * t105 + (2 * pkin(7) * mrSges(5,3)) + ((pkin(3) ^ 2 + pkin(7) ^ 2) * m(5)); mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,3); m(3); m(5) * pkin(3) + mrSges(4,1); -pkin(7) * m(5) + mrSges(4,2) - mrSges(5,3); m(4) + m(5); t120 * t117 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t115; t117 * t103 + t120 * t118 + Ifges(5,4); t111 * Ifges(6,5) - t110 * Ifges(6,6) + Ifges(5,5); t110 * Ifges(6,5) + t111 * Ifges(6,6) + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t112 + t113) * m(7)) + 0.2e1 * t114 * pkin(4) + t119; mrSges(5,1) + t114; t111 * t104 + t110 * t106 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
