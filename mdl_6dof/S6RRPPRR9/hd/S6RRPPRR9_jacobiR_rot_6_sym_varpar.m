% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:17:52
% EndTime: 2019-02-22 11:17:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (77->30), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
t137 = sin(qJ(2));
t138 = sin(qJ(1));
t141 = cos(qJ(2));
t142 = cos(qJ(1));
t154 = cos(pkin(6));
t143 = t142 * t154;
t128 = t137 * t143 + t138 * t141;
t136 = sin(qJ(5));
t140 = cos(qJ(5));
t134 = sin(pkin(6));
t149 = t134 * t142;
t122 = -t128 * t136 + t140 * t149;
t127 = t138 * t137 - t141 * t143;
t135 = sin(qJ(6));
t139 = cos(qJ(6));
t156 = t122 * t135 - t127 * t139;
t155 = t122 * t139 + t127 * t135;
t151 = t134 * t136;
t150 = t134 * t140;
t148 = t135 * t136;
t147 = t135 * t141;
t146 = t136 * t139;
t145 = t139 * t141;
t144 = t138 * t154;
t121 = t128 * t140 + t136 * t149;
t130 = -t137 * t144 + t142 * t141;
t129 = -t142 * t137 - t141 * t144;
t126 = t137 * t151 + t154 * t140;
t125 = -t154 * t136 + t137 * t150;
t120 = t130 * t136 + t138 * t150;
t119 = -t130 * t140 + t138 * t151;
t118 = t120 * t139 + t129 * t135;
t117 = -t120 * t135 + t129 * t139;
t1 = [t155, t129 * t146 - t130 * t135, 0, 0, -t119 * t139, t117; t118, -t127 * t146 - t128 * t135, 0, 0, t121 * t139, t156; 0 (-t135 * t137 + t136 * t145) * t134, 0, 0, t125 * t139, -t126 * t135 + t134 * t145; -t156, -t129 * t148 - t130 * t139, 0, 0, t119 * t135, -t118; t117, t127 * t148 - t128 * t139, 0, 0, -t121 * t135, t155; 0 (-t136 * t147 - t137 * t139) * t134, 0, 0, -t125 * t135, -t126 * t139 - t134 * t147; t121, -t129 * t140, 0, 0, t120, 0; t119, t127 * t140, 0, 0, -t122, 0; 0, -t141 * t150, 0, 0, t126, 0;];
JR_rot  = t1;
