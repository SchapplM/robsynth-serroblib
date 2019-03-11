% Return the minimum parameter vector for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPPR1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t150 = (pkin(8) ^ 2);
t151 = (pkin(5) ^ 2);
t159 = 2 * pkin(8) * mrSges(7,3) + Ifges(7,2);
t135 = Ifges(5,2) + Ifges(6,3) + (t150 + t151) * m(7) + t159;
t136 = m(7) * t150 + Ifges(5,1) + Ifges(6,1) + t159;
t144 = sin(pkin(10));
t139 = t144 ^ 2;
t146 = cos(pkin(10));
t141 = t146 ^ 2;
t149 = Ifges(5,4) - Ifges(6,5);
t158 = t144 * t146;
t154 = t149 * t158;
t132 = t139 * t135 + t141 * t136 + Ifges(4,1) - 0.2e1 * t154;
t138 = t151 * m(7) + Ifges(4,2) + Ifges(6,2) + Ifges(5,3);
t160 = t132 - t138;
t145 = sin(pkin(9));
t147 = cos(pkin(9));
t157 = t145 * t147;
t140 = t145 ^ 2;
t142 = t147 ^ 2;
t156 = t142 - t140;
t155 = pkin(8) * m(7) + mrSges(7,3);
t137 = t155 * pkin(5) + Ifges(6,4) + Ifges(5,5);
t148 = Ifges(5,6) - Ifges(6,6);
t134 = -t146 * t137 + t144 * t148 + Ifges(4,4);
t153 = t134 * t157;
t152 = t147 * mrSges(4,1) - t145 * mrSges(4,2);
t133 = -t144 * t137 - t146 * t148 + Ifges(4,6);
t131 = Ifges(4,5) + (t141 - t139) * t149 + (-t135 + t136) * t158;
t1 = [Ifges(2,3) + Ifges(3,2) + t140 * t132 + 0.2e1 * t153 + t142 * t138 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t160 * t156 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t153; t156 * t134 + t160 * t157 + Ifges(3,4); t147 * t131 - t145 * t133 + Ifges(3,5); t145 * t131 + t147 * t133 + Ifges(3,6); 0.2e1 * pkin(2) * t152 + t141 * t135 + t139 * t136 + Ifges(3,3) + Ifges(4,3) + 0.2e1 * t154; mrSges(3,1) + t152; t145 * mrSges(4,1) + t147 * mrSges(4,2) + mrSges(3,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t155; mrSges(6,3); m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
