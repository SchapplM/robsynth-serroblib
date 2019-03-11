% Return the minimum parameter vector for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t117 = m(5) + m(6);
t127 = (-Ifges(6,2) - Ifges(7,3));
t113 = (m(4) + t117);
t126 = 2 * pkin(9) * mrSges(6,3) - t127;
t125 = m(6) * pkin(9) + mrSges(6,3);
t124 = -pkin(8) * t117 - mrSges(5,3);
t123 = (mrSges(4,3) - t124);
t110 = -pkin(7) * t113 + mrSges(3,2) - t123;
t111 = pkin(2) * t113 + mrSges(3,1);
t115 = sin(pkin(10));
t116 = cos(pkin(10));
t122 = -t110 * t115 + t111 * t116;
t121 = pkin(3) ^ 2;
t120 = (pkin(4) ^ 2);
t119 = pkin(8) ^ 2;
t118 = pkin(9) ^ 2;
t112 = t119 + t121;
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + (t120 * m(6)) + Ifges(5,2) + (2 * pkin(8) * mrSges(5,3)) + t112 * t117 + (2 * pkin(7) * t123) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t113) + 0.2e1 * t122 * pkin(1); mrSges(2,1) + t122; t110 * t116 + t111 * t115 + mrSges(2,2); m(3) + t113; Ifges(4,1) - Ifges(4,2) + (-t112 + t119) * t117; Ifges(4,4); pkin(3) * t124 + Ifges(4,5); Ifges(4,6); t117 * t121 + Ifges(4,3); pkin(3) * t117 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + ((t118 - t120) * m(6)) + t126; pkin(4) * t125 + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t118 + t120) * m(6) + t126; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t125; Ifges(6,1) + Ifges(7,1) + t127; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
