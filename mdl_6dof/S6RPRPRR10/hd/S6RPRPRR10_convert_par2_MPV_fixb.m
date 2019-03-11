% Return the minimum parameter vector for
% S6RPRPRR10
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
% MPV [31x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t129 = -pkin(9) * m(7) - mrSges(7,3);
t108 = (mrSges(6,3) - t129);
t115 = (m(6) + m(7));
t128 = -pkin(8) * t115 - t108;
t127 = (pkin(7) * m(4));
t117 = (pkin(9) ^ 2);
t119 = (pkin(5) ^ 2);
t125 = (Ifges(6,2) + (t117 + t119) * m(7));
t113 = sin(pkin(10));
t114 = cos(pkin(10));
t124 = t113 * t114;
t123 = Ifges(5,4) * t124;
t120 = (pkin(4) ^ 2);
t122 = (t120 * t115 + Ifges(4,2) + Ifges(5,3));
t121 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t108 + Ifges(7,2) + t125;
t118 = pkin(8) ^ 2;
t111 = t114 ^ 2;
t110 = t113 ^ 2;
t105 = t128 * pkin(4) + Ifges(5,5);
t104 = t118 * t115 + Ifges(5,1) + t121;
t103 = Ifges(5,2) + (t118 + t120) * t115 + t121;
t1 = [Ifges(3,1) + Ifges(2,3) + (2 * mrSges(4,3) + t127) * pkin(7) + t122; mrSges(2,1); mrSges(2,2); mrSges(3,2) - mrSges(4,3) - t127; mrSges(3,3); m(3) + m(4); t110 * t103 + t111 * t104 + Ifges(4,1) - t122 - 0.2e1 * t123; t113 * Ifges(5,6) - t114 * t105 + Ifges(4,4); Ifges(4,5) + (t111 - t110) * Ifges(5,4) + (-t103 + t104) * t124; -t114 * Ifges(5,6) - t113 * t105 + Ifges(4,6); t111 * t103 + t110 * t104 + Ifges(4,3) + 0.2e1 * t123; mrSges(4,1); mrSges(4,2); pkin(4) * t115 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t128; m(5) + t115; m(7) * t117 + Ifges(6,1) - t125; Ifges(6,4); t129 * pkin(5) + Ifges(6,5); Ifges(6,6); t119 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
