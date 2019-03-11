% Return the minimum parameter vector for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPPR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t140 = (pkin(8) ^ 2);
t141 = (pkin(5) ^ 2);
t149 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t127 = Ifges(6,2) + (t140 + t141) * m(7) + t149;
t128 = m(7) * t140 + Ifges(6,1) + t149;
t136 = sin(pkin(10));
t131 = t136 ^ 2;
t138 = cos(pkin(10));
t133 = t138 ^ 2;
t148 = t136 * t138;
t145 = Ifges(6,4) * t148;
t124 = t133 * t127 + t131 * t128 + Ifges(4,2) + Ifges(5,3) + 0.2e1 * t145;
t130 = t141 * m(7) + Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t150 = -t124 + t130;
t137 = sin(pkin(9));
t139 = cos(pkin(9));
t147 = t137 * t139;
t132 = t137 ^ 2;
t134 = t139 ^ 2;
t146 = t134 - t132;
t144 = -pkin(8) * m(7) - mrSges(7,3);
t129 = t144 * pkin(5) + Ifges(6,5);
t126 = -t138 * Ifges(6,6) - t136 * t129 + Ifges(4,4) + Ifges(5,6);
t143 = t126 * t147;
t142 = t139 * mrSges(4,1) - t137 * mrSges(4,2);
t125 = -t136 * Ifges(6,6) + t138 * t129 - Ifges(5,4) + Ifges(4,5);
t123 = Ifges(4,6) - Ifges(5,5) - (t133 - t131) * Ifges(6,4) + (t127 - t128) * t148;
t1 = [Ifges(2,3) + Ifges(3,2) + t132 * t130 + 0.2e1 * t143 + t134 * t124 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t150 * t146 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t143; t146 * t126 + t150 * t147 + Ifges(3,4); -t137 * t123 + t139 * t125 + Ifges(3,5); t139 * t123 + t137 * t125 + Ifges(3,6); 0.2e1 * pkin(2) * t142 + t131 * t127 + t133 * t128 + Ifges(5,1) + Ifges(3,3) + Ifges(4,3) - 0.2e1 * t145; mrSges(3,1) + t142; t137 * mrSges(4,1) + t139 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t144; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
