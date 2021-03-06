% Return the minimum parameter vector for
% S6RPRRRP2
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
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t117 = m(5) + m(6);
t127 = (-Ifges(6,2) - Ifges(7,2));
t118 = (pkin(9) ^ 2);
t120 = (pkin(4) ^ 2);
t126 = (Ifges(5,2) + (t118 + t120) * m(6));
t113 = (m(4) + t117);
t125 = -pkin(9) * m(6) - mrSges(6,3);
t111 = (mrSges(5,3) - t125);
t124 = pkin(8) * t117 + t111;
t107 = -pkin(7) * t113 + mrSges(3,2) - mrSges(4,3);
t109 = pkin(2) * t113 + mrSges(3,1);
t115 = sin(pkin(10));
t116 = cos(pkin(10));
t123 = -t115 * t107 + t116 * t109;
t122 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t111 + t126 - t127;
t121 = pkin(3) ^ 2;
t119 = pkin(8) ^ 2;
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + t121 * t117 + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t113) + 0.2e1 * t123 * pkin(1); mrSges(2,1) + t123; t116 * t107 + t115 * t109 + mrSges(2,2); m(3) + t113; Ifges(4,1) - Ifges(4,2) + (t119 - t121) * t117 + t122; t124 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t119 + t121) * t117 + t122; pkin(3) * t117 + mrSges(4,1); mrSges(4,2) - t124; t118 * m(6) + Ifges(5,1) - t126; Ifges(5,4); t125 * pkin(4) + Ifges(5,5); Ifges(5,6); t120 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t127; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
