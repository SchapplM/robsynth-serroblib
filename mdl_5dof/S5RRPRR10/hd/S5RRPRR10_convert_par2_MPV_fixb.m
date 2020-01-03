% Return the minimum parameter vector for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPRR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t138 = (m(5) + m(6));
t141 = (pkin(8) ^ 2);
t142 = (pkin(4) ^ 2);
t152 = (t142 * m(6) + Ifges(5,2));
t148 = 2 * pkin(8) * mrSges(5,3) + t152;
t125 = t141 * t138 + Ifges(4,1) + t148;
t143 = (pkin(3) ^ 2);
t129 = t143 * t138 + Ifges(4,2);
t154 = t125 - t129;
t153 = (m(3) * pkin(7));
t151 = 2 * pkin(9) * mrSges(6,3) + Ifges(6,2);
t135 = sin(pkin(10));
t137 = cos(pkin(10));
t150 = t135 * t137;
t130 = t135 ^ 2;
t132 = t137 ^ 2;
t149 = t132 - t130;
t147 = pkin(9) * m(6) + mrSges(6,3);
t145 = pkin(8) * t138 + mrSges(5,3);
t126 = t145 * pkin(3) + Ifges(4,4);
t146 = t126 * t150;
t127 = mrSges(4,2) - t145;
t128 = pkin(3) * t138 + mrSges(4,1);
t144 = -t135 * t127 + t137 * t128;
t140 = pkin(9) ^ 2;
t136 = sin(pkin(5));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (0.2e1 * t146 + t130 * t125 + t132 * t129 + Ifges(3,2) + ((2 * mrSges(3,3) + t153) * pkin(7))) * t136 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t153) * t136; t149 * t154 + Ifges(3,1) - Ifges(3,2) - 0.4e1 * t146; t149 * t126 + t150 * t154 + Ifges(3,4); t137 * Ifges(4,5) - t135 * Ifges(4,6) + Ifges(3,5); t135 * Ifges(4,5) + t137 * Ifges(4,6) + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + ((t141 + t143) * t138) + 0.2e1 * t144 * pkin(2) + t148; mrSges(3,1) + t144; t137 * t127 + t135 * t128 + mrSges(3,2); mrSges(4,3); m(4) + t138; m(6) * t140 + Ifges(5,1) + t151 - t152; t147 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t140 + t142) * m(6) + t151; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t147; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
