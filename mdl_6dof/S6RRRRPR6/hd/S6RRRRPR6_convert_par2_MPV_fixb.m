% Return the minimum parameter vector for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t139 = (m(4) + m(5));
t141 = pkin(9) ^ 2;
t144 = pkin(3) ^ 2;
t155 = Ifges(4,2) + (t141 + t144) * m(5);
t154 = 2 * pkin(10) * mrSges(7,3) + Ifges(7,2);
t137 = sin(pkin(11));
t138 = cos(pkin(11));
t153 = t137 * t138;
t152 = Ifges(6,4) * t153;
t140 = pkin(10) ^ 2;
t143 = (pkin(5) ^ 2);
t123 = Ifges(6,2) + (t140 + t143) * m(7) + t154;
t125 = m(7) * t140 + Ifges(6,1) + t154;
t132 = t137 ^ 2;
t133 = t138 ^ 2;
t151 = t133 * t123 + t132 * t125 + Ifges(5,2);
t150 = -pkin(9) * m(5) - mrSges(5,3);
t149 = -pkin(10) * m(7) - mrSges(7,3);
t129 = (mrSges(4,3) - t150);
t148 = pkin(8) * t139 + t129;
t130 = m(7) * pkin(5) + mrSges(6,1);
t147 = -t137 * mrSges(6,2) + t138 * t130;
t146 = (2 * pkin(9) * mrSges(5,3)) + (2 * pkin(8) * t129) + t151 + 0.2e1 * t152 + t155;
t145 = (pkin(2) ^ 2);
t142 = pkin(8) ^ 2;
t134 = (m(3) + t139);
t126 = t149 * pkin(5) + Ifges(6,5);
t1 = [Ifges(2,3) + Ifges(3,2) + t145 * t139 + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t134; pkin(1) * t134 + mrSges(2,1); -pkin(7) * t134 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + ((t142 - t145) * t139) + t146; t148 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + ((t142 + t145) * t139) + t146; pkin(2) * t139 + mrSges(3,1); mrSges(3,2) - t148; t141 * m(5) + Ifges(4,1) - t155; Ifges(4,4); t150 * pkin(3) + Ifges(4,5); Ifges(4,6); t144 * m(5) + Ifges(4,3); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2); t132 * t123 + t133 * t125 + Ifges(5,1) - t151 - 0.4e1 * t152; Ifges(5,4) + (t133 - t132) * Ifges(6,4) + (-t123 + t125) * t153; -t137 * Ifges(6,6) + t138 * t126 + Ifges(5,5); t138 * Ifges(6,6) + t137 * t126 + Ifges(5,6); (t143 * m(7)) + 0.2e1 * pkin(4) * t147 + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t147; t138 * mrSges(6,2) + t137 * t130 + mrSges(5,2); mrSges(6,3) - t149; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
