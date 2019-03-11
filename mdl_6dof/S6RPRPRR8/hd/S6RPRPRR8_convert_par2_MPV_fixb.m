% Return the minimum parameter vector for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t119 = (m(6) + m(7));
t121 = (pkin(8) ^ 2);
t128 = -pkin(9) * m(7) - mrSges(7,3);
t112 = (mrSges(6,3) - t128);
t120 = (pkin(9) ^ 2);
t122 = (pkin(5) ^ 2);
t131 = (Ifges(6,2) + (t120 + t122) * m(7));
t125 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t112 + Ifges(7,2) + t131;
t105 = t121 * t119 + Ifges(5,1) + t125;
t123 = (pkin(4) ^ 2);
t110 = t123 * t119 + Ifges(5,2);
t133 = t105 - t110;
t132 = (pkin(7) * m(4));
t117 = sin(pkin(10));
t118 = cos(pkin(10));
t130 = t117 * t118;
t114 = t117 ^ 2;
t115 = t118 ^ 2;
t129 = t115 - t114;
t126 = pkin(8) * t119 + t112;
t106 = t126 * pkin(4) + Ifges(5,4);
t127 = t106 * t130;
t107 = mrSges(5,2) - t126;
t109 = pkin(4) * t119 + mrSges(5,1);
t124 = -t117 * t107 + t118 * t109;
t1 = [0.2e1 * t127 + t114 * t105 + t115 * t110 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + ((2 * mrSges(4,3) + t132) * pkin(7)); mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t132; mrSges(3,3); m(3) + m(4); t133 * t129 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t127; t129 * t106 + t133 * t130 + Ifges(4,4); t118 * Ifges(5,5) - t117 * Ifges(5,6) + Ifges(4,5); t117 * Ifges(5,5) + t118 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t121 + t123) * t119) + 0.2e1 * t124 * pkin(3) + t125; mrSges(4,1) + t124; t118 * t107 + t117 * t109 + mrSges(4,2); mrSges(5,3); m(5) + t119; m(7) * t120 + Ifges(6,1) - t131; Ifges(6,4); t128 * pkin(5) + Ifges(6,5); Ifges(6,6); t122 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
