% Return the minimum parameter vector for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% MPV [33x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRRP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t121 = (m(5) + m(6));
t134 = -pkin(9) * t121 - mrSges(5,3);
t114 = (mrSges(4,3) - t134);
t119 = (m(4) + t121);
t133 = -pkin(8) * t119 - t114;
t131 = (-Ifges(6,2) - Ifges(7,2));
t130 = 2 * pkin(10) * mrSges(6,3) - t131;
t129 = pkin(10) * m(6) + mrSges(6,3);
t128 = (mrSges(3,3) - t133);
t127 = (pkin(2) ^ 2);
t126 = (pkin(3) ^ 2);
t125 = (pkin(4) ^ 2);
t124 = (pkin(8) ^ 2);
t123 = (pkin(9) ^ 2);
t122 = pkin(10) ^ 2;
t118 = (t124 + t127);
t117 = (t123 + t126);
t115 = (m(3) + t119);
t1 = [Ifges(2,3) + Ifges(3,2) + Ifges(4,2) + t125 * m(6) + Ifges(5,2) + 2 * pkin(9) * mrSges(5,3) + t117 * t121 + 2 * pkin(8) * t114 + t118 * t119 + 2 * pkin(7) * t128 + (pkin(1) ^ 2 + pkin(7) ^ 2) * t115; pkin(1) * t115 + mrSges(2,1); -pkin(7) * t115 + mrSges(2,2) - t128; Ifges(3,1) - Ifges(3,2) + (-t118 + t124) * t119; Ifges(3,4); t133 * pkin(2) + Ifges(3,5); Ifges(3,6); t127 * t119 + Ifges(3,3); pkin(2) * t119 + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) + (-t117 + t123) * t121; Ifges(4,4); t134 * pkin(3) + Ifges(4,5); Ifges(4,6); t126 * t121 + Ifges(4,3); pkin(3) * t121 + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) + (t122 - t125) * m(6) + t130; t129 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t122 + t125) * m(6) + t130; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t129; Ifges(6,1) + Ifges(7,1) + t131; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); 2 * pkin(5) * mrSges(7,1) + Ifges(6,3) + Ifges(7,3); mrSges(6,1) + mrSges(7,1); mrSges(6,2) + mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
