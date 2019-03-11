% Return the minimum parameter vector for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRPRR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t152 = -pkin(9) * m(7) - mrSges(7,3);
t126 = (mrSges(6,3) - t152);
t137 = (m(6) + m(7));
t151 = -pkin(8) * t137 - t126;
t139 = (pkin(9) ^ 2);
t142 = (pkin(5) ^ 2);
t149 = (Ifges(6,2) + (t139 + t142) * m(7));
t133 = sin(pkin(11));
t135 = cos(pkin(11));
t148 = t133 * t135;
t143 = (pkin(4) ^ 2);
t147 = (t143 * t137 + Ifges(4,2) + Ifges(5,3));
t146 = Ifges(5,4) * t148;
t145 = 2 * pkin(7) * mrSges(4,3) + t147;
t144 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t126 + Ifges(7,2) + t149;
t141 = pkin(7) ^ 2;
t140 = pkin(8) ^ 2;
t136 = cos(pkin(10));
t134 = sin(pkin(10));
t130 = t135 ^ 2;
t129 = t133 ^ 2;
t123 = t151 * pkin(4) + Ifges(5,5);
t122 = t140 * t137 + Ifges(5,1) + t144;
t121 = Ifges(5,2) + (t140 + t143) * t137 + t144;
t1 = [Ifges(2,3) + t136 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t141) * m(4) + t145) + (0.2e1 * t136 * Ifges(3,4) + (t141 * m(4) + Ifges(3,1) + t145) * t134) * t134; mrSges(2,1); mrSges(2,2); m(4) * pkin(2) + mrSges(3,1); mrSges(3,2); pkin(7) * m(4) + mrSges(3,3) + mrSges(4,3); m(3) + m(4); t129 * t121 + t130 * t122 + Ifges(4,1) - 0.2e1 * t146 - t147; t133 * Ifges(5,6) - t135 * t123 + Ifges(4,4); Ifges(4,5) + (t130 - t129) * Ifges(5,4) + (-t121 + t122) * t148; -t135 * Ifges(5,6) - t133 * t123 + Ifges(4,6); t130 * t121 + t129 * t122 + Ifges(4,3) + 0.2e1 * t146; mrSges(4,1); mrSges(4,2); pkin(4) * t137 + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t151; m(5) + t137; m(7) * t139 + Ifges(6,1) - t149; Ifges(6,4); t152 * pkin(5) + Ifges(6,5); Ifges(6,6); t142 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
