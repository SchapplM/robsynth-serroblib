% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:24:12
% EndTime: 2019-02-22 09:24:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
t136 = sin(pkin(11));
t142 = cos(pkin(6));
t154 = t136 * t142;
t137 = sin(pkin(7));
t138 = sin(pkin(6));
t153 = t137 * t138;
t152 = t137 * t142;
t141 = cos(pkin(7));
t151 = t138 * t141;
t139 = cos(pkin(12));
t150 = t139 * t141;
t140 = cos(pkin(11));
t149 = t140 * t142;
t135 = sin(pkin(12));
t131 = -t136 * t135 + t139 * t149;
t148 = t131 * t141 - t140 * t153;
t133 = -t140 * t135 - t139 * t154;
t147 = t133 * t141 + t136 * t153;
t146 = cos(qJ(3));
t145 = cos(qJ(4));
t144 = sin(qJ(3));
t143 = sin(qJ(4));
t134 = -t135 * t154 + t140 * t139;
t132 = t135 * t149 + t136 * t139;
t130 = -t139 * t153 + t142 * t141;
t129 = -t133 * t137 + t136 * t151;
t128 = -t131 * t137 - t140 * t151;
t127 = t144 * t152 + (t135 * t146 + t144 * t150) * t138;
t126 = t146 * t152 + (-t135 * t144 + t146 * t150) * t138;
t125 = t134 * t146 + t147 * t144;
t124 = -t134 * t144 + t147 * t146;
t123 = t132 * t146 + t148 * t144;
t122 = -t132 * t144 + t148 * t146;
t1 = [0, 0, t125, 0, 0, 0; 0, 0, t123, 0, 0, 0; 0, 0, t127, 0, 0, 0; 0, 0, -t124 * t145, t125 * t143 - t129 * t145, 0, 0; 0, 0, -t122 * t145, t123 * t143 - t128 * t145, 0, 0; 0, 0, -t126 * t145, t127 * t143 - t130 * t145, 0, 0; 0, 0, t124 * t143, t125 * t145 + t129 * t143, 0, 0; 0, 0, t122 * t143, t123 * t145 + t128 * t143, 0, 0; 0, 0, t126 * t143, t127 * t145 + t130 * t143, 0, 0;];
JR_rot  = t1;
