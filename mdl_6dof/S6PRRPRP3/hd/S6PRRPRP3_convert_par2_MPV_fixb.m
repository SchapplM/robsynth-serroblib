% Return the minimum parameter vector for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRP3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t155 = (-Ifges(6,2) - Ifges(7,3));
t146 = sin(pkin(11));
t147 = cos(pkin(11));
t154 = t146 * t147;
t153 = 2 * pkin(9) * mrSges(6,3) - t155;
t152 = Ifges(5,4) * t154;
t151 = -pkin(9) * m(6) - mrSges(6,3);
t149 = (pkin(4) ^ 2);
t150 = t149 * m(6) + Ifges(4,2) + Ifges(5,3);
t148 = pkin(9) ^ 2;
t144 = t147 ^ 2;
t143 = t146 ^ 2;
t142 = t151 * pkin(4) + Ifges(5,5);
t141 = t148 * m(6) + Ifges(5,1) + t153;
t140 = Ifges(5,2) + (t148 + t149) * m(6) + t153;
t1 = [m(2) + m(3) + m(4); Ifges(3,3) + 2 * pkin(8) * mrSges(4,3) + (pkin(2) ^ 2 + pkin(8) ^ 2) * m(4) + t150; m(4) * pkin(2) + mrSges(3,1); -pkin(8) * m(4) + mrSges(3,2) - mrSges(4,3); t143 * t140 + t144 * t141 + Ifges(4,1) - t150 - 0.2e1 * t152; t146 * Ifges(5,6) - t147 * t142 + Ifges(4,4); Ifges(4,5) + (t144 - t143) * Ifges(5,4) + (-t140 + t141) * t154; -t147 * Ifges(5,6) - t146 * t142 + Ifges(4,6); t144 * t140 + t143 * t141 + Ifges(4,3) + 0.2e1 * t152; mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t151; m(5) + m(6); Ifges(6,1) + Ifges(7,1) + t155; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
