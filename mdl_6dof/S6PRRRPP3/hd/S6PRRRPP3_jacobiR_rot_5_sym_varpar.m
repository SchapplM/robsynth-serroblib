% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:53
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPP3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:53:26
% EndTime: 2019-02-22 09:53:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t132 = sin(pkin(6));
t136 = sin(qJ(3));
t148 = t132 * t136;
t139 = cos(qJ(3));
t147 = t132 * t139;
t140 = cos(qJ(2));
t146 = t132 * t140;
t134 = cos(pkin(6));
t137 = sin(qJ(2));
t145 = t134 * t137;
t144 = t134 * t140;
t135 = sin(qJ(4));
t143 = t135 * t139;
t138 = cos(qJ(4));
t142 = t138 * t139;
t141 = t139 * t140;
t133 = cos(pkin(10));
t131 = sin(pkin(10));
t129 = t134 * t136 + t137 * t147;
t128 = t134 * t139 - t137 * t148;
t127 = -t131 * t145 + t133 * t140;
t126 = t131 * t144 + t133 * t137;
t125 = t131 * t140 + t133 * t145;
t124 = t131 * t137 - t133 * t144;
t123 = t127 * t139 + t131 * t148;
t122 = -t127 * t136 + t131 * t147;
t121 = t125 * t139 - t133 * t148;
t120 = -t125 * t136 - t133 * t147;
t1 = [0, -t126 * t136, t123, 0, 0, 0; 0, -t124 * t136, t121, 0, 0, 0; 0, t136 * t146, t129, 0, 0, 0; 0, t126 * t142 - t127 * t135, -t122 * t138, t123 * t135 - t126 * t138, 0, 0; 0, t124 * t142 - t125 * t135, -t120 * t138, t121 * t135 - t124 * t138, 0, 0; 0 (-t135 * t137 - t138 * t141) * t132, -t128 * t138, t129 * t135 + t138 * t146, 0, 0; 0, -t126 * t143 - t127 * t138, t122 * t135, t123 * t138 + t126 * t135, 0, 0; 0, -t124 * t143 - t125 * t138, t120 * t135, t121 * t138 + t124 * t135, 0, 0; 0 (t135 * t141 - t137 * t138) * t132, t128 * t135, t129 * t138 - t135 * t146, 0, 0;];
JR_rot  = t1;
