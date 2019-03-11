% Return the minimum parameter vector for
% S6RRPRPR11
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRPR11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t141 = (pkin(9) ^ 2);
t143 = (pkin(5) ^ 2);
t152 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t129 = Ifges(6,2) + (t141 + t143) * m(7) + t152;
t131 = m(7) * t141 + Ifges(6,1) + t152;
t139 = sin(pkin(10));
t135 = t139 ^ 2;
t140 = cos(pkin(10));
t136 = t140 ^ 2;
t149 = t136 * t129 + t135 * t131 + Ifges(5,2);
t151 = t139 * t140;
t150 = Ifges(6,4) * t151;
t146 = (2 * pkin(8) * mrSges(5,3)) + t149 + 0.2e1 * t150;
t155 = -t146 - Ifges(3,2) - Ifges(4,3);
t148 = -pkin(8) * m(5) - mrSges(5,3);
t147 = -pkin(9) * m(7) - mrSges(7,3);
t133 = m(7) * pkin(5) + mrSges(6,1);
t145 = -t139 * mrSges(6,2) + t140 * t133;
t144 = (pkin(3) ^ 2);
t142 = pkin(8) ^ 2;
t134 = t142 + t144;
t132 = t147 * pkin(5) + Ifges(6,5);
t1 = [Ifges(2,3) + (t134 * m(5)) + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)) - t155; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) + Ifges(4,2) + ((-t134 + t144) * m(5)) + t155; Ifges(3,4) + Ifges(4,6); t148 * pkin(3) - Ifges(4,4) + Ifges(3,5); Ifges(3,6) - Ifges(4,5); (t142 * m(5)) + Ifges(4,1) + Ifges(3,3) + t146; mrSges(3,1); mrSges(3,2); m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) + t148; mrSges(4,3); m(4) + m(5); t135 * t129 + t136 * t131 + Ifges(5,1) - t149 - 0.4e1 * t150; Ifges(5,4) + (t136 - t135) * Ifges(6,4) + (-t129 + t131) * t151; -t139 * Ifges(6,6) + t140 * t132 + Ifges(5,5); t140 * Ifges(6,6) + t139 * t132 + Ifges(5,6); (t143 * m(7)) + 0.2e1 * pkin(4) * t145 + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t145; t140 * mrSges(6,2) + t139 * t133 + mrSges(5,2); mrSges(6,3) - t147; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
