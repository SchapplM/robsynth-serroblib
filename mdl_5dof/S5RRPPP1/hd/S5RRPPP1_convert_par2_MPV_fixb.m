% Return the minimum parameter vector for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRPPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t126 = Ifges(4,6) - Ifges(5,5) - Ifges(6,4);
t128 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t135 = sin(pkin(8));
t137 = cos(pkin(8));
t140 = t126 * t137 + t128 * t135;
t136 = sin(pkin(5));
t138 = cos(pkin(5));
t146 = t136 * t138;
t152 = t140 * t146;
t125 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t132 = t136 ^ 2;
t134 = t138 ^ 2;
t127 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t130 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t131 = t135 ^ 2;
t133 = t137 ^ 2;
t129 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t143 = t137 * t135 * t129;
t139 = t133 * t127 + t131 * t130 + 0.2e1 * t143;
t151 = t132 * t125 + t139 * t134 + Ifges(3,2);
t150 = (t133 - t131) * t129;
t147 = t128 * t137;
t145 = 0.2e1 * t152;
t141 = (-t127 + t130) * t137;
t1 = [Ifges(2,3) + (2 * pkin(7) * mrSges(3,3)) + ((pkin(1) ^ 2 + pkin(7) ^ 2) * m(3)) - 0.2e1 * t152 + t151; m(3) * pkin(1) + mrSges(2,1); -pkin(7) * m(3) + mrSges(2,2) - mrSges(3,3); t131 * t127 + t133 * t130 + Ifges(3,1) - 0.2e1 * t143 + t145 - t151; t138 * t150 - t136 * t147 + Ifges(3,4) + (t126 * t136 + t138 * t141) * t135; t136 * t150 + t138 * t147 + Ifges(3,5) + (-t126 * t138 + t136 * t141) * t135; Ifges(3,6) + t140 * (t134 - t132) + (-t125 + t139) * t146; t134 * t125 + t139 * t132 + Ifges(3,3) + t145; mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
