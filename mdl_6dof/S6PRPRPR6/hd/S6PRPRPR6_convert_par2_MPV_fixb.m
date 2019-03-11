% Return the minimum parameter vector for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRPRPR6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t143 = (pkin(8) * m(5));
t142 = 2 * pkin(9) * mrSges(7,3) + Ifges(7,2);
t134 = sin(pkin(11));
t135 = cos(pkin(11));
t141 = t134 * t135;
t140 = Ifges(6,4) * t141;
t139 = -pkin(9) * m(7) - mrSges(7,3);
t137 = (pkin(5) ^ 2);
t138 = t137 * m(7) + Ifges(5,2) + Ifges(6,3);
t136 = pkin(9) ^ 2;
t132 = t135 ^ 2;
t131 = t134 ^ 2;
t130 = t139 * pkin(5) + Ifges(6,5);
t129 = m(7) * t136 + Ifges(6,1) + t142;
t128 = Ifges(6,2) + (t136 + t137) * m(7) + t142;
t1 = [m(2) + m(3); Ifges(4,1) + Ifges(3,3) + (2 * mrSges(5,3) + t143) * pkin(8) + t138; mrSges(3,1); mrSges(3,2); mrSges(4,2) - mrSges(5,3) - t143; mrSges(4,3); m(4) + m(5); t131 * t128 + t132 * t129 + Ifges(5,1) - t138 - 0.2e1 * t140; t134 * Ifges(6,6) - t135 * t130 + Ifges(5,4); Ifges(5,5) + (t132 - t131) * Ifges(6,4) + (-t128 + t129) * t141; -t135 * Ifges(6,6) - t134 * t130 + Ifges(5,6); t132 * t128 + t131 * t129 + Ifges(5,3) + 0.2e1 * t140; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t139; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
