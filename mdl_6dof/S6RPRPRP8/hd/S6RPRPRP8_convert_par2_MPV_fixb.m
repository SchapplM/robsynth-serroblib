% Return the minimum parameter vector for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRP8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t114 = (pkin(8) ^ 2);
t122 = (-Ifges(6,2) - Ifges(7,3));
t119 = 2 * pkin(8) * mrSges(6,3) - t122;
t104 = t114 * m(6) + Ifges(5,1) + t119;
t115 = (pkin(4) ^ 2);
t107 = t115 * m(6) + Ifges(5,2);
t124 = t104 - t107;
t123 = (pkin(7) * m(4));
t112 = sin(pkin(9));
t113 = cos(pkin(9));
t121 = t112 * t113;
t109 = t112 ^ 2;
t110 = t113 ^ 2;
t120 = t110 - t109;
t118 = pkin(8) * m(6) + mrSges(6,3);
t105 = t118 * pkin(4) + Ifges(5,4);
t117 = t105 * t121;
t106 = mrSges(5,2) - t118;
t108 = m(6) * pkin(4) + mrSges(5,1);
t116 = -t112 * t106 + t113 * t108;
t1 = [0.2e1 * t117 + t109 * t104 + t110 * t107 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t123) * pkin(7)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t123; mrSges(3,3); m(3) + m(4); t124 * t120 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t117; t120 * t105 + t124 * t121 + Ifges(4,4); t113 * Ifges(5,5) - t112 * Ifges(5,6) + Ifges(4,5); t112 * Ifges(5,5) + t113 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t114 + t115) * m(6)) + 0.2e1 * t116 * pkin(3) + t119; mrSges(4,1) + t116; t113 * t106 + t112 * t108 + mrSges(4,2); mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t122; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
