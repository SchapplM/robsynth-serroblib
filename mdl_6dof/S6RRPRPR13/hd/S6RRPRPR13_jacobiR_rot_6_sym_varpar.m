% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:30
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR13_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:29:57
% EndTime: 2019-02-22 11:29:57
% DurationCPUTime: 0.10s
% Computational Cost: add. (115->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t137 = cos(pkin(6));
t142 = cos(qJ(2));
t143 = cos(qJ(1));
t144 = t143 * t142;
t139 = sin(qJ(2));
t140 = sin(qJ(1));
t147 = t140 * t139;
t127 = -t137 * t144 + t147;
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t136 = sin(pkin(6));
t149 = t136 * t143;
t122 = -t127 * t138 + t141 * t149;
t145 = t143 * t139;
t146 = t140 * t142;
t128 = t137 * t145 + t146;
t135 = pkin(11) + qJ(6);
t133 = sin(t135);
t134 = cos(t135);
t158 = t122 * t133 + t128 * t134;
t157 = t122 * t134 - t128 * t133;
t154 = t133 * t138;
t153 = t134 * t138;
t152 = t136 * t139;
t151 = t136 * t141;
t150 = t136 * t142;
t148 = t138 * t139;
t121 = t127 * t141 + t138 * t149;
t130 = -t137 * t147 + t144;
t129 = t137 * t146 + t145;
t126 = t137 * t141 - t138 * t150;
t125 = -t137 * t138 - t141 * t150;
t120 = t129 * t138 + t140 * t151;
t119 = t140 * t136 * t138 - t129 * t141;
t118 = t120 * t134 + t130 * t133;
t117 = -t120 * t133 + t130 * t134;
t1 = [t157, -t129 * t133 + t130 * t153, 0, -t119 * t134, 0, t117; t118, -t127 * t133 + t128 * t153, 0, t121 * t134, 0, t158; 0 (t133 * t142 + t134 * t148) * t136, 0, t125 * t134, 0, -t126 * t133 + t134 * t152; -t158, -t129 * t134 - t130 * t154, 0, t119 * t133, 0, -t118; t117, -t127 * t134 - t128 * t154, 0, -t121 * t133, 0, t157; 0 (-t133 * t148 + t134 * t142) * t136, 0, -t125 * t133, 0, -t126 * t134 - t133 * t152; t121, -t130 * t141, 0, t120, 0, 0; t119, -t128 * t141, 0, -t122, 0, 0; 0, -t139 * t151, 0, t126, 0, 0;];
JR_rot  = t1;
