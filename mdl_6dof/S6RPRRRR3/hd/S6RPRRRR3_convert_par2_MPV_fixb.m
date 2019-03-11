% Return the minimum parameter vector for
% S6RPRRRR3
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t143 = -pkin(10) * m(7) - mrSges(7,3);
t120 = (mrSges(6,3) - t143);
t128 = (m(6) + m(7));
t142 = -pkin(9) * t128 - t120;
t131 = (pkin(9) ^ 2);
t134 = (pkin(4) ^ 2);
t140 = (Ifges(5,2) + (t131 + t134) * t128);
t130 = (pkin(10) ^ 2);
t133 = (pkin(5) ^ 2);
t139 = (Ifges(6,2) + (t130 + t133) * m(7));
t124 = m(5) + t128;
t121 = (m(4) + t124);
t114 = (mrSges(5,3) - t142);
t138 = pkin(8) * t124 + t114;
t115 = -pkin(7) * t121 + mrSges(3,2) - mrSges(4,3);
t116 = pkin(2) * t121 + mrSges(3,1);
t126 = sin(pkin(11));
t127 = cos(pkin(11));
t137 = -t126 * t115 + t127 * t116;
t136 = 2 * pkin(10) * mrSges(7,3) + 2 * pkin(8) * t114 + 2 * pkin(9) * t120 + Ifges(7,2) + t139 + t140;
t135 = pkin(3) ^ 2;
t132 = pkin(8) ^ 2;
t1 = [Ifges(2,3) + Ifges(3,3) + Ifges(4,2) + t135 * t124 + (2 * pkin(7) * mrSges(4,3)) + ((pkin(2) ^ 2 + pkin(7) ^ 2) * t121) + 0.2e1 * t137 * pkin(1); mrSges(2,1) + t137; t127 * t115 + t126 * t116 + mrSges(2,2); m(3) + t121; Ifges(4,1) - Ifges(4,2) + (t132 - t135) * t124 + t136; t138 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t132 + t135) * t124 + t136; pkin(3) * t124 + mrSges(4,1); mrSges(4,2) - t138; t131 * t128 + Ifges(5,1) - t140; Ifges(5,4); pkin(4) * t142 + Ifges(5,5); Ifges(5,6); t134 * t128 + Ifges(5,3); pkin(4) * t128 + mrSges(5,1); mrSges(5,2); m(7) * t130 + Ifges(6,1) - t139; Ifges(6,4); pkin(5) * t143 + Ifges(6,5); Ifges(6,6); t133 * m(7) + Ifges(6,3); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
