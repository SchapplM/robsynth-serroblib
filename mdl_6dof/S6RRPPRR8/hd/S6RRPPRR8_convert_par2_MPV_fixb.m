% Return the minimum parameter vector for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPPRR8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t142 = (m(6) + m(7));
t143 = (pkin(9) ^ 2);
t145 = (pkin(5) ^ 2);
t153 = (Ifges(6,2) + (t143 + t145) * m(7));
t138 = sin(pkin(10));
t139 = cos(pkin(10));
t152 = t138 * t139;
t151 = -pkin(9) * m(7) - mrSges(7,3);
t141 = Ifges(4,4) - Ifges(5,5);
t150 = t141 * t152;
t133 = (mrSges(6,3) - t151);
t149 = pkin(8) * t142 + t133;
t148 = 2 * pkin(9) * mrSges(7,3) + 2 * pkin(8) * t133 + Ifges(7,2) + t153;
t146 = (pkin(4) ^ 2);
t147 = (t146 * t142 + Ifges(3,2) + Ifges(5,2) + Ifges(4,3));
t144 = pkin(8) ^ 2;
t140 = Ifges(4,6) - Ifges(5,6);
t136 = t139 ^ 2;
t135 = t138 ^ 2;
t130 = t149 * pkin(4) + Ifges(5,4) + Ifges(4,5);
t129 = t144 * t142 + Ifges(4,1) + Ifges(5,1) + t148;
t128 = Ifges(4,2) + Ifges(5,3) + (t144 + t146) * t142 + t148;
t1 = [Ifges(2,3) + 2 * pkin(7) * mrSges(3,3) + (pkin(1) ^ 2 + pkin(7) ^ 2) * m(3) + t147; m(3) * pkin(1) + mrSges(2,1); -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3); t135 * t128 + t136 * t129 + Ifges(3,1) - t147 - 0.2e1 * t150; -t139 * t130 + t138 * t140 + Ifges(3,4); Ifges(3,5) + (t136 - t135) * t141 + (-t128 + t129) * t152; -t138 * t130 - t139 * t140 + Ifges(3,6); t136 * t128 + t135 * t129 + Ifges(3,3) + 0.2e1 * t150; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); pkin(4) * t142 + mrSges(5,1); mrSges(5,2) - t149; mrSges(5,3); m(5) + t142; m(7) * t143 + Ifges(6,1) - t153; Ifges(6,4); t151 * pkin(5) + Ifges(6,5); Ifges(6,6); m(7) * t145 + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
