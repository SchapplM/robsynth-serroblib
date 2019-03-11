% Return the minimum parameter vector for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRRPRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t140 = (m(6) + m(7));
t142 = (pkin(10) ^ 2);
t145 = (pkin(5) ^ 2);
t153 = (Ifges(6,2) + (t142 + t145) * m(7));
t143 = (pkin(9) ^ 2);
t146 = (pkin(4) ^ 2);
t152 = (Ifges(4,2) + Ifges(5,3) + (t143 + t146) * t140);
t151 = pkin(8) * m(4) + mrSges(4,3);
t150 = -pkin(10) * m(7) - mrSges(7,3);
t135 = (mrSges(6,3) - t150);
t149 = pkin(9) * t140 + t135;
t148 = 2 * pkin(8) * mrSges(4,3) + 2 * pkin(10) * mrSges(7,3) + 2 * pkin(9) * t135 + Ifges(7,2) + t152 + t153;
t147 = (pkin(2) ^ 2);
t144 = pkin(8) ^ 2;
t141 = (m(3) + m(4));
t1 = [Ifges(2,3) + t147 * m(4) + Ifges(3,2) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * t141; pkin(1) * t141 + mrSges(2,1); -pkin(7) * t141 + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + (t144 - t147) * m(4) + t148; t151 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t144 + t147) * m(4) + t148; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t151; t143 * t140 + Ifges(4,1) + Ifges(5,1) - t152; Ifges(4,4) - Ifges(5,5); t149 * pkin(4) + Ifges(5,4) + Ifges(4,5); Ifges(4,6) - Ifges(5,6); t146 * t140 + Ifges(5,2) + Ifges(4,3); mrSges(4,1); mrSges(4,2); pkin(4) * t140 + mrSges(5,1); mrSges(5,2) - t149; mrSges(5,3); m(5) + t140; m(7) * t142 + Ifges(6,1) - t153; Ifges(6,4); t150 * pkin(5) + Ifges(6,5); Ifges(6,6); t145 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
