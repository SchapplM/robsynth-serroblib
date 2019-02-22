% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:55:34
% EndTime: 2019-02-22 11:55:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (109->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->37)
t134 = cos(pkin(6));
t136 = sin(qJ(2));
t140 = cos(qJ(1));
t143 = t140 * t136;
t137 = sin(qJ(1));
t139 = cos(qJ(2));
t144 = t137 * t139;
t126 = t134 * t143 + t144;
t135 = sin(qJ(3));
t138 = cos(qJ(3));
t133 = sin(pkin(6));
t147 = t133 * t140;
t118 = t126 * t135 + t138 * t147;
t142 = t140 * t139;
t145 = t137 * t136;
t125 = -t134 * t142 + t145;
t132 = pkin(11) + qJ(6);
t130 = sin(t132);
t131 = cos(t132);
t156 = -t118 * t130 - t125 * t131;
t155 = t118 * t131 - t125 * t130;
t152 = t130 * t135;
t151 = t131 * t135;
t150 = t133 * t135;
t149 = t133 * t138;
t148 = t133 * t139;
t146 = t135 * t139;
t141 = -t126 * t138 + t135 * t147;
t128 = -t134 * t145 + t142;
t127 = t134 * t144 + t143;
t124 = t134 * t135 + t136 * t149;
t123 = -t134 * t138 + t136 * t150;
t122 = t128 * t138 + t137 * t150;
t121 = t128 * t135 - t137 * t149;
t117 = t121 * t130 + t127 * t131;
t116 = t121 * t131 - t127 * t130;
t1 = [t156, -t127 * t152 + t128 * t131, t122 * t130, 0, 0, t116; t117, -t125 * t152 + t126 * t131, -t141 * t130, 0, 0, t155; 0 (t130 * t146 + t131 * t136) * t133, t124 * t130, 0, 0, t123 * t131 + t130 * t148; -t155, -t127 * t151 - t128 * t130, t122 * t131, 0, 0, -t117; t116, -t125 * t151 - t126 * t130, -t141 * t131, 0, 0, t156; 0 (-t130 * t136 + t131 * t146) * t133, t124 * t131, 0, 0, -t123 * t130 + t131 * t148; t141, -t127 * t138, -t121, 0, 0, 0; t122, -t125 * t138, -t118, 0, 0, 0; 0, t138 * t148, -t123, 0, 0, 0;];
JR_rot  = t1;
