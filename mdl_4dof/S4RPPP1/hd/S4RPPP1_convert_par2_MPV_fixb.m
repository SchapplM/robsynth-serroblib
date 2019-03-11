% Return the minimum parameter vector for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [15x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RPPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t81 = cos(pkin(4));
t80 = cos(pkin(6));
t79 = sin(pkin(4));
t78 = sin(pkin(6));
t1 = [Ifges(2,3) + ((Ifges(3,3) + Ifges(4,1) + Ifges(5,1)) * t81 + 0.2e1 * (t78 * (Ifges(3,5) - Ifges(4,4) + Ifges(5,5)) + t80 * (Ifges(3,6) - Ifges(4,5) - Ifges(5,4))) * t79) * t81 + (t80 ^ 2 * (Ifges(3,2) + Ifges(4,3) + Ifges(5,2)) + (0.2e1 * t80 * (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) + (Ifges(3,1) + Ifges(4,2) + Ifges(5,3)) * t78) * t78) * t79 ^ 2; mrSges(2,1); mrSges(2,2); mrSges(3,1); mrSges(3,2); mrSges(3,3); m(3); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV  = t1;
