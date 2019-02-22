% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:28
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:28:06
% EndTime: 2019-02-22 11:28:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (122->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t137 = cos(pkin(6));
t139 = sin(qJ(2));
t143 = cos(qJ(1));
t146 = t143 * t139;
t140 = sin(qJ(1));
t142 = cos(qJ(2));
t148 = t140 * t142;
t129 = t137 * t146 + t148;
t135 = pkin(11) + qJ(4);
t133 = sin(t135);
t134 = cos(t135);
t136 = sin(pkin(6));
t151 = t136 * t143;
t121 = t129 * t133 + t134 * t151;
t145 = t143 * t142;
t149 = t140 * t139;
t128 = -t137 * t145 + t149;
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t159 = -t121 * t138 - t128 * t141;
t158 = t121 * t141 - t128 * t138;
t155 = t133 * t138;
t154 = t133 * t141;
t153 = t136 * t139;
t152 = t136 * t140;
t150 = t138 * t142;
t147 = t141 * t142;
t144 = -t129 * t134 + t133 * t151;
t131 = -t137 * t149 + t145;
t130 = t137 * t148 + t146;
t127 = t137 * t133 + t134 * t153;
t126 = t133 * t153 - t137 * t134;
t125 = t131 * t134 + t133 * t152;
t124 = t131 * t133 - t134 * t152;
t120 = t124 * t138 + t130 * t141;
t119 = t124 * t141 - t130 * t138;
t1 = [t159, -t130 * t155 + t131 * t141, 0, t125 * t138, 0, t119; t120, -t128 * t155 + t129 * t141, 0, -t144 * t138, 0, t158; 0 (t133 * t150 + t139 * t141) * t136, 0, t127 * t138, 0, t126 * t141 + t136 * t150; -t158, -t130 * t154 - t131 * t138, 0, t125 * t141, 0, -t120; t119, -t128 * t154 - t129 * t138, 0, -t144 * t141, 0, t159; 0 (t133 * t147 - t138 * t139) * t136, 0, t127 * t141, 0, -t126 * t138 + t136 * t147; t144, -t130 * t134, 0, -t124, 0, 0; t125, -t128 * t134, 0, -t121, 0, 0; 0, t136 * t142 * t134, 0, -t126, 0, 0;];
JR_rot  = t1;
