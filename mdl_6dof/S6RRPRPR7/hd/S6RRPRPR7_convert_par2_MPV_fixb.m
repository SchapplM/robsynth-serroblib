% Return the minimum parameter vector for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% MPV [30x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t144 = (pkin(9) ^ 2);
t154 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t133 = m(7) * t144 + Ifges(6,1) + t154;
t146 = (pkin(5) ^ 2);
t136 = m(7) * t146 + Ifges(6,2);
t156 = t133 - t136;
t155 = -Ifges(3,2) - Ifges(4,3);
t142 = sin(pkin(10));
t143 = cos(pkin(10));
t153 = t142 * t143;
t139 = t142 ^ 2;
t140 = t143 ^ 2;
t152 = t140 - t139;
t151 = m(5) * pkin(8) + mrSges(5,3);
t150 = m(7) * pkin(9) + mrSges(7,3);
t134 = pkin(5) * t150 + Ifges(6,4);
t149 = t134 * t153;
t135 = mrSges(6,2) - t150;
t137 = m(7) * pkin(5) + mrSges(6,1);
t148 = -t135 * t142 + t137 * t143;
t147 = pkin(3) ^ 2;
t145 = pkin(8) ^ 2;
t138 = t145 + t147;
t1 = [Ifges(2,3) + Ifges(5,2) + t139 * t133 + 0.2e1 * t149 + t140 * t136 + (2 * pkin(8) * mrSges(5,3)) + t138 * m(5) + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)) - t155; m(3) * pkin(1) + mrSges(2,1); -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) + Ifges(4,1) + (-t138 + t145) * m(5) + t155; Ifges(3,4) - Ifges(4,5); pkin(3) * t151 + Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,6); m(5) * t147 + Ifges(4,2) + Ifges(3,3); mrSges(3,1); mrSges(3,2); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t151; mrSges(4,3); m(4) + m(5); t152 * t156 + Ifges(5,1) - Ifges(5,2) - 0.4e1 * t149; t134 * t152 + t153 * t156 + Ifges(5,4); Ifges(6,5) * t143 - Ifges(6,6) * t142 + Ifges(5,5); Ifges(6,5) * t142 + Ifges(6,6) * t143 + Ifges(5,6); Ifges(5,3) + Ifges(6,3) + ((t144 + t146) * m(7)) + 0.2e1 * t148 * pkin(4) + t154; mrSges(5,1) + t148; t135 * t143 + t137 * t142 + mrSges(5,2); mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
