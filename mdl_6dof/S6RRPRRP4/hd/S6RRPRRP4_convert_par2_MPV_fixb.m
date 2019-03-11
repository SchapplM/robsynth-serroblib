% Return the minimum parameter vector for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t139 = (m(5) + m(6));
t141 = (pkin(8) ^ 2);
t148 = -pkin(9) * m(6) - mrSges(6,3);
t132 = (mrSges(5,3) - t148);
t140 = (pkin(9) ^ 2);
t142 = (pkin(4) ^ 2);
t151 = (Ifges(5,2) + (t140 + t142) * m(6));
t152 = (-Ifges(6,2) - Ifges(7,3));
t144 = 2 * pkin(9) * mrSges(6,3) + 2 * pkin(8) * t132 + t151 - t152;
t125 = t141 * t139 + Ifges(4,1) + t144;
t143 = (pkin(3) ^ 2);
t130 = t143 * t139 + Ifges(4,2);
t153 = t125 - t130;
t137 = sin(pkin(10));
t138 = cos(pkin(10));
t150 = t137 * t138;
t134 = t137 ^ 2;
t135 = t138 ^ 2;
t149 = t135 - t134;
t146 = pkin(8) * t139 + t132;
t126 = t146 * pkin(3) + Ifges(4,4);
t147 = t126 * t150;
t127 = mrSges(4,2) - t146;
t129 = pkin(3) * t139 + mrSges(4,1);
t145 = -t137 * t127 + t138 * t129;
t1 = [Ifges(2,3) + Ifges(3,2) + t134 * t125 + 0.2e1 * t147 + t135 * t130 + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)); m(3) * pkin(1) + mrSges(2,1); -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3); t153 * t149 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t147; t149 * t126 + t153 * t150 + Ifges(3,4); t138 * Ifges(4,5) - t137 * Ifges(4,6) + Ifges(3,5); t137 * Ifges(4,5) + t138 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t141 + t143) * t139) + 0.2e1 * t145 * pkin(2) + t144; mrSges(3,1) + t145; t138 * t127 + t137 * t129 + mrSges(3,2); mrSges(4,3); m(4) + t139; t140 * m(6) + Ifges(5,1) - t151; Ifges(5,4); pkin(4) * t148 + Ifges(5,5); Ifges(5,6); m(6) * t142 + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) + Ifges(7,1) + t152; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
