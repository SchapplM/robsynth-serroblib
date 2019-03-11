% Return the minimum parameter vector for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t133 = (m(6) + m(7));
t128 = (m(5) + t133);
t149 = -pkin(8) * t128 - mrSges(5,3);
t136 = (pkin(8) ^ 2);
t140 = (pkin(3) ^ 2);
t148 = (Ifges(4,2) + (t136 + t140) * t128);
t139 = (pkin(4) ^ 2);
t147 = (t139 * t133 + Ifges(5,2));
t134 = (pkin(10) ^ 2);
t138 = (pkin(5) ^ 2);
t146 = (Ifges(6,2) + (t134 + t138) * m(7));
t145 = (mrSges(4,3) - t149);
t124 = (m(4) + t128);
t144 = -pkin(10) * m(7) - mrSges(7,3);
t122 = (mrSges(6,3) - t144);
t143 = pkin(9) * t133 + t122;
t142 = 2 * pkin(10) * mrSges(7,3) + 2 * pkin(9) * t122 + Ifges(7,2) + t146;
t141 = 2 * pkin(8) * mrSges(5,3) + 2 * pkin(7) * t145 + t147 + t148;
t137 = pkin(7) ^ 2;
t135 = pkin(9) ^ 2;
t132 = cos(pkin(11));
t131 = sin(pkin(11));
t1 = [Ifges(2,3) + t132 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t137) * t124 + t141) + (0.2e1 * t132 * Ifges(3,4) + (t137 * t124 + Ifges(3,1) + t141) * t131) * t131; mrSges(2,1); mrSges(2,2); pkin(2) * t124 + mrSges(3,1); mrSges(3,2); pkin(7) * t124 + mrSges(3,3) + t145; m(3) + t124; t136 * t128 + Ifges(4,1) - t148; Ifges(4,4); t149 * pkin(3) + Ifges(4,5); Ifges(4,6); t140 * t128 + Ifges(4,3); pkin(3) * t128 + mrSges(4,1); mrSges(4,2); t135 * t133 + Ifges(5,1) + t142 - t147; t143 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t135 + t139) * t133 + t142; pkin(4) * t133 + mrSges(5,1); mrSges(5,2) - t143; m(7) * t134 + Ifges(6,1) - t146; Ifges(6,4); t144 * pkin(5) + Ifges(6,5); Ifges(6,6); t138 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
